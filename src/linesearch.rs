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
    const MAX_PARABOLA_REFINEMENTS: usize = 12;

    let start_value = params.fval();
    let start = params.vec();
    let mut scratch = start.clone();

    let unit_value = evaluate_step(fcn, start, step, 1.0, &mut scratch);

    let mut best_lambda = 0.0;
    let mut best_value = start_value;
    if unit_value < start_value {
        best_lambda = 1.0;
        best_value = unit_value;
    }

    let mut trial_lambda = 1.0;
    let curvature =
        2.0 * (unit_value - start_value - gdel * trial_lambda) / (trial_lambda * trial_lambda);
    if curvature.abs() > prec.eps2() {
        trial_lambda = -gdel / curvature;
    } else {
        trial_lambda = INITIAL_STEP_LIMIT;
    }

    trial_lambda = trial_lambda.clamp(NEAR_EXISTING_POINT_TOLERANCE, INITIAL_STEP_LIMIT);

    if (trial_lambda - 1.0).abs() < NEAR_EXISTING_POINT_TOLERANCE && unit_value < start_value {
        return MnParabolaPoint::new(best_lambda, best_value);
    }

    let trial_value = evaluate_step(fcn, start, step, trial_lambda, &mut scratch);

    if trial_value < best_value {
        best_lambda = trial_lambda;
        best_value = trial_value;
    }

    let mut points = sorted_points([
        MnParabolaPoint::new(0.0, start_value),
        MnParabolaPoint::new(1.0, unit_value),
        MnParabolaPoint::new(trial_lambda, trial_value),
    ]);

    let upper_step_limit = (BEST_STEP_EXPANSION * best_lambda.abs()).max(INITIAL_STEP_LIMIT);

    for _ in 0..MAX_PARABOLA_REFINEMENTS {
        let parabola = from_3_points(points[0], points[1], points[2]);

        if parabola.a() < prec.eps2() {
            break;
        }

        trial_lambda = parabola.min();

        if trial_lambda > upper_step_limit {
            trial_lambda = upper_step_limit;
        }
        if trial_lambda < LOWER_STEP_LIMIT {
            trial_lambda = LOWER_STEP_LIMIT;
        }
        if trial_lambda < 0.0 && parabola.y(0.0) < parabola.y(trial_lambda) {
            break;
        }

        let min_spacing = NEAR_EXISTING_POINT_TOLERANCE * trial_lambda.abs().max(1.0);
        if points
            .iter()
            .any(|point| (trial_lambda - point.x).abs() < min_spacing)
        {
            break;
        }

        let trial_value = evaluate_step(fcn, start, step, trial_lambda, &mut scratch);

        if trial_value < best_value {
            best_lambda = trial_lambda;
            best_value = trial_value;
        }

        replace_worst_point(&mut points, MnParabolaPoint::new(trial_lambda, trial_value));
        points = sorted_points(points);

        if (best_value - start_value).abs() < start_value.abs() * prec.eps() {
            break;
        }
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
}
