//! Parabolic line search along a step direction.
//!
//! Performs one-dimensional parabolic interpolation to choose a step length
//! along a search direction. Returns the best `(lambda, f(lambda))` pair found.

use crate::minimum::parameters::MinimumParameters;
use crate::mn_fcn::MnFcn;
use crate::parabola::{MnParabolaPoint, from_3_points};
use crate::precision::MnMachinePrecision;
use nalgebra::DVector;

/// Perform a parabolic line search along `step` from `params`.
///
/// `gdel` is the directional derivative: step · gradient (should be negative
/// for a descent direction). Returns `MnParabolaPoint { x: λ_opt, y: f_opt }`.
pub fn mn_linesearch(
    fcn: &MnFcn,
    params: &MinimumParameters,
    step: &DVector<f64>,
    gdel: f64,
    prec: &MnMachinePrecision,
) -> MnParabolaPoint {
    const LOWER_STEP_LIMIT: f64 = -100.0;
    const NEAR_EXISTING_POINT_TOLERANCE: f64 = 0.05;
    const INITIAL_STEP_LIMIT: f64 = 5.0;
    const BEST_STEP_EXPANSION: f64 = 2.0;
    // Longest forward step of the open search, as in MnLineSearch.cxx:62.
    const OVERALL_STEP_LIMIT: f64 = 1000.0;
    const MAX_PARABOLA_REFINEMENTS: usize = 12;

    let start_value = params.fval();
    let start = params.vec();
    let mut scratch = start.clone();

    // Smallest step that still changes a parameter. MnLineSearch.cxx:74-85.
    let mut slamin = 0.0;
    for (x_i, step_i) in start.iter().zip(step.iter()) {
        // exact zero only: no finite ratio
        if *step_i == 0.0 {
            continue;
        }
        let ratio = (x_i / step_i).abs();
        if slamin == 0.0 || ratio < slamin {
            slamin = ratio;
        }
    }
    if slamin.abs() < prec.eps() {
        slamin = prec.eps();
    }
    slamin *= prec.eps2();

    let toler = NEAR_EXISTING_POINT_TOLERANCE;
    let mut toler8 = toler;
    let mut slamax = INITIAL_STEP_LIMIT;
    let mut overall = OVERALL_STEP_LIMIT;
    let mut undral = LOWER_STEP_LIMIT;

    // ROOT starts at 1; shortening, failed halvings and refinements share 12.
    let mut niter = 1usize;
    let unit_value = evaluate_step(fcn, start, step, 1.0, &mut scratch);
    niter += 1;

    let mut best_lambda = 0.0;
    let mut best_value = start_value;
    if unit_value < start_value {
        best_lambda = 1.0;
        best_value = unit_value;
    }

    // Shorten the trial while it does not lower the function
    // (MnLineSearch.cxx:107-210). A full Newton step can raise the value.
    // The next trials must get shorter until one step improves it.
    let mut flast = unit_value;
    let mut slam = 1.0;
    let mut previous_lambda = 1.0;
    let mut previous_value = unit_value;
    let second_value = loop {
        let denom = 2.0 * (flast - start_value - gdel * slam) / (slam * slam);
        // exact zero only: no finite ratio
        if denom != 0.0 {
            slam = -gdel / denom;
        } else {
            slam = 1.0;
        }
        if slam < 0.0 {
            slam = slamax;
        }
        if slam > slamax {
            slam = slamax;
        }
        if slam < toler8 {
            slam = toler8;
        }
        if slam < slamin {
            return MnParabolaPoint::new(best_lambda, best_value);
        }
        if (slam - 1.0).abs() < toler8 && previous_value < start_value {
            return MnParabolaPoint::new(best_lambda, best_value);
        }
        if (slam - 1.0).abs() < toler8 {
            slam = 1.0 + toler8;
        }

        let value = evaluate_step(fcn, start, step, slam, &mut scratch);
        niter += 1;
        if value < best_value {
            best_value = value;
            best_lambda = slam;
        }
        let no_improvement = (start_value - best_value).abs() < best_value.abs() * prec.eps();
        if no_improvement && niter < MAX_PARABOLA_REFINEMENTS {
            flast = value;
            toler8 = toler * slam;
            overall = slam - toler8;
            slamax = overall;
            previous_lambda = slam;
            previous_value = flast;
            continue;
        }
        if niter >= MAX_PARABOLA_REFINEMENTS {
            return MnParabolaPoint::new(best_lambda, best_value);
        }
        break value;
    };

    let mut points = sorted_points([
        MnParabolaPoint::new(0.0, start_value),
        MnParabolaPoint::new(previous_lambda, previous_value),
        MnParabolaPoint::new(slam, second_value),
    ]);

    // Forward cap of a vertex. The two-point search may already have lowered
    // it. It grows with the best step on every pass (MnLineSearch.cxx:219).
    let mut upper_step_limit = slamax;
    // Length of the step taken from the best point when the parabola has no
    // minimum. It starts from the same cap and grows with the best step.
    let mut open_step_limit = slamax;

    while niter < MAX_PARABOLA_REFINEMENTS {
        let parabola = from_3_points(points[0], points[1], points[2]);
        // The cap never passes the overall step limit (MnLineSearch.cxx:239-245).
        upper_step_limit = upper_step_limit
            .max(BEST_STEP_EXPANSION * best_lambda.abs())
            .min(OVERALL_STEP_LIMIT);
        let is_open_step = parabola.a() < prec.eps2();
        let mut trial_lambda;

        if is_open_step {
            // A parabola without a minimum gives no vertex to test. Instead of
            // stopping, MnLineSearch.cxx:224-229 steps from the best point
            // along the slope, so a search that starts on a limit leaves it.
            open_step_limit = open_step_limit.max(BEST_STEP_EXPANSION * best_lambda.abs());
            let slope = 2.0 * parabola.a() * best_lambda + parabola.b();
            trial_lambda = if slope < 0.0 {
                best_lambda + open_step_limit
            } else {
                best_lambda - open_step_limit
            };
            if trial_lambda > 0.0 {
                if trial_lambda > overall {
                    trial_lambda = overall;
                }
            } else if trial_lambda < undral {
                trial_lambda = undral;
            }
        } else {
            trial_lambda = parabola.min();

            // Keep the vertex inside the window around the best step.
            if trial_lambda > best_lambda + upper_step_limit {
                trial_lambda = best_lambda + upper_step_limit;
            }
            if trial_lambda < best_lambda - upper_step_limit {
                trial_lambda = best_lambda - upper_step_limit;
            }
            if trial_lambda > 0.0 {
                if trial_lambda > overall {
                    trial_lambda = overall;
                }
            } else if trial_lambda < undral {
                trial_lambda = undral;
            }
            if trial_lambda < 0.0 && parabola.y(0.0) < parabola.y(trial_lambda) {
                break;
            }
        }

        // Halve a trial that is worse than the three stored points.
        let trial_value = loop {
            let min_spacing = toler8.max((toler8 * trial_lambda).abs());
            if points
                .iter()
                .any(|point| (trial_lambda - point.x).abs() < min_spacing)
            {
                return MnParabolaPoint::new(best_lambda, best_value);
            }
            let value = evaluate_step(fcn, start, step, trial_lambda, &mut scratch);
            let worse_than_all = value > points[0].y && value > points[1].y && value > points[2].y;
            if !worse_than_all {
                break value;
            }
            if trial_lambda > best_lambda {
                overall = overall.min(trial_lambda - toler8);
            }
            if trial_lambda < best_lambda {
                undral = undral.max(trial_lambda + toler8);
            }
            trial_lambda = 0.5 * (trial_lambda + best_lambda);
            niter += 1;
            if niter >= MAX_PARABOLA_REFINEMENTS {
                return MnParabolaPoint::new(best_lambda, best_value);
            }
        };

        if trial_value < best_value {
            best_lambda = trial_lambda;
            best_value = trial_value;
        }

        replace_worst_point(&mut points, MnParabolaPoint::new(trial_lambda, trial_value));
        points = sorted_points(points);

        if (best_value - start_value).abs() < start_value.abs() * prec.eps() {
            break;
        }
        niter += 1;
    }

    MnParabolaPoint::new(best_lambda, best_value)
}

fn evaluate_step(
    fcn: &MnFcn,
    start: &DVector<f64>,
    direction: &DVector<f64>,
    lambda: f64,
    scratch: &mut DVector<f64>,
) -> f64 {
    scratch.copy_from(start);
    scratch.axpy(lambda, direction, 1.0);
    fcn.call(scratch.as_slice())
}

fn sorted_points(mut points: [MnParabolaPoint; 3]) -> [MnParabolaPoint; 3] {
    points.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal));
    points
}

fn replace_worst_point(points: &mut [MnParabolaPoint; 3], new_point: MnParabolaPoint) {
    let mut worst_index = 0;
    for index in 1..points.len() {
        let order = points[index]
            .y
            .partial_cmp(&points[worst_index].y)
            .unwrap_or(std::cmp::Ordering::Equal);
        if order != std::cmp::Ordering::Less {
            worst_index = index;
        }
    }
    points[worst_index] = new_point;
}

pub fn do_eval(fcn: &MnFcn, x: &[f64]) -> f64 {
    fcn.call(x)
}

/// Compatibility wrapper for callers using cubic-search terminology.
pub fn cubic_search(
    fcn: &MnFcn,
    params: &MinimumParameters,
    step: &DVector<f64>,
    gdel: f64,
    prec: &MnMachinePrecision,
) -> MnParabolaPoint {
    mn_linesearch(fcn, params, step, gdel, prec)
}

/// Compatibility wrapper for callers using Brent-search terminology.
pub fn brent_search(
    fcn: &MnFcn,
    params: &MinimumParameters,
    step: &DVector<f64>,
    gdel: f64,
    prec: &MnMachinePrecision,
) -> MnParabolaPoint {
    mn_linesearch(fcn, params, step, gdel, prec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fcn::FCN;
    use crate::parameter::MinuitParameter;
    use crate::user_transformation::MnUserTransformation;

    struct Quadratic;
    impl FCN for Quadratic {
        fn value(&self, p: &[f64]) -> f64 {
            p[0] * p[0]
        }
    }

    #[test]
    fn linesearch_quadratic() {
        let params = vec![MinuitParameter::new(0, "x", 2.0, 0.1)];
        let trafo = MnUserTransformation::new(params);
        let fcn = MnFcn::new(&Quadratic, &trafo);

        // Start at x=2, step direction = -1 (downhill)
        let start = MinimumParameters::new(DVector::from_vec(vec![2.0]), 4.0);
        let step = DVector::from_vec(vec![-1.0]);
        let gdel = step.dot(&DVector::from_vec(vec![4.0])); // grad = 2x = 4 at x=2
        let prec = MnMachinePrecision::new();

        let result = mn_linesearch(&fcn, &start, &step, gdel, &prec);

        // Optimal step for f(2-λ) = (2-λ)² is λ=2, giving f=0
        assert!(
            (result.x - 2.0).abs() < 0.1,
            "line search should find lambda near 2, got {}",
            result.x
        );
        assert!(result.y < 4.0, "line search should improve: f={}", result.y);
        assert!(
            result.y < 0.1,
            "line search should find near-minimum: f={}",
            result.y
        );
    }

    #[test]
    fn linesearch_non_descent_direction_does_not_make_result_worse() {
        let params = vec![MinuitParameter::new(0, "x", 1.0, 0.1)];
        let trafo = MnUserTransformation::new(params);
        let fcn = MnFcn::new(&Quadratic, &trafo);

        let start = MinimumParameters::new(DVector::from_vec(vec![1.0]), 1.0);
        let step = DVector::from_vec(vec![1.0]);
        let gdel = 2.0;
        let prec = MnMachinePrecision::new();

        let result = mn_linesearch(&fcn, &start, &step, gdel, &prec);

        assert!(result.x.is_finite());
        assert!(
            result.y <= 1.0,
            "non-descent search should keep the best known value, got {}",
            result.y
        );
    }

    #[test]
    fn linesearch_clamps_large_interpolated_step() {
        let params = vec![MinuitParameter::new(0, "x", 10.0, 0.1)];
        let trafo = MnUserTransformation::new(params);
        let fcn = MnFcn::new(&Quadratic, &trafo);

        let start = MinimumParameters::new(DVector::from_vec(vec![10.0]), 100.0);
        let step = DVector::from_vec(vec![-1.0]);
        let gdel = -20.0;
        let prec = MnMachinePrecision::new();

        let result = mn_linesearch(&fcn, &start, &step, gdel, &prec);

        assert!(result.x.is_finite());
        assert!(
            result.x.abs() <= 20.0,
            "line search ran away to {}",
            result.x
        );
        assert!(result.y < 100.0);
    }

    struct Constant;
    impl FCN for Constant {
        fn value(&self, _p: &[f64]) -> f64 {
            7.0
        }
    }

    #[test]
    fn linesearch_degenerate_parabola_is_finite_and_no_worse() {
        let params = vec![MinuitParameter::new(0, "x", 0.0, 0.1)];
        let trafo = MnUserTransformation::new(params);
        let fcn = MnFcn::new(&Constant, &trafo);

        let start = MinimumParameters::new(DVector::from_vec(vec![0.0]), 7.0);
        let step = DVector::from_vec(vec![1.0]);
        let prec = MnMachinePrecision::new();

        let result = mn_linesearch(&fcn, &start, &step, 0.0, &prec);

        assert!(result.x.is_finite());
        assert!(result.y.is_finite());
        assert_eq!(result.y, 7.0);
    }

    /// f(t) = (t - 10000)^2 along the step. The minimum lies at 10000, far
    /// beyond the step limit of 1000 (MnLineSearch.cxx:239-245).
    struct FarQuadratic;
    impl FCN for FarQuadratic {
        fn value(&self, p: &[f64]) -> f64 {
            (p[0] - 10000.0).powi(2)
        }
    }

    #[test]
    fn linesearch_vertex_cap_stays_within_overall_step_limit() {
        let params = vec![MinuitParameter::new(0, "x", 0.0, 0.1)];
        let trafo = MnUserTransformation::new(params);
        let fcn = MnFcn::new(&FarQuadratic, &trafo);

        let start = MinimumParameters::new(DVector::from_vec(vec![0.0]), 1.0e8);
        let step = DVector::from_vec(vec![1.0]);
        // gdel = step . grad = 1 * 2 * (0 - 10000)
        let gdel = -20000.0;
        let prec = MnMachinePrecision::new();

        let result = mn_linesearch(&fcn, &start, &step, gdel, &prec);

        assert!(
            result.x <= 1000.0,
            "accepted step {} is above the overall limit 1000",
            result.x
        );
        assert!(result.y < 1.0e8, "line search should improve: {}", result.y);
    }
}
