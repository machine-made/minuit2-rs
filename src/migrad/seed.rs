//! Migrad seed.
//!
//! Evaluates the function and the gradient, then sets the covariance
//! diagonal to 1/g2. Strategy 2 replaces that covariance with the
//! inverse Hessian.

use nalgebra::{DMatrix, DVector};

use crate::fcn::FCNGradient;
use crate::gradient::{
    AnalyticalGradientCalculator, InitialGradientCalculator, Numerical2PGradientCalculator,
};
use crate::hesse::calculator as hesse_calculator;
use crate::linesearch::mn_linesearch;
use crate::minimum::error::MinimumError;
use crate::minimum::gradient::FunctionGradient;
use crate::minimum::parameters::MinimumParameters;
use crate::minimum::seed::MinimumSeed;
use crate::minimum::state::MinimumState;
use crate::mn_fcn::MnFcn;
use crate::strategy::MnStrategy;
use crate::user_transformation::MnUserTransformation;

pub struct MigradSeedGenerator;

impl MigradSeedGenerator {
    /// Generate seed using numerical gradients (central differences).
    pub fn generate(
        fcn: &MnFcn,
        trafo: &MnUserTransformation,
        strategy: &MnStrategy,
    ) -> MinimumSeed {
        let n = trafo.variable_parameters();
        let eps = trafo.precision().eps();

        // 1. Get initial internal parameter values
        let int_values = trafo.initial_internal_values();
        let int_vec = DVector::from_vec(int_values.clone());

        // 2. Evaluate FCN at starting point
        let fval = fcn.call(&int_values);
        let params = MinimumParameters::new(int_vec, fval);

        // 3. Compute heuristic gradient (no FCN calls — just from step sizes)
        let heuristic_calc = InitialGradientCalculator::new(*strategy);
        let heuristic_grad = heuristic_calc.compute(fcn, &params, trafo);

        // 4. Compute numerical gradient (2-point central differences)
        let numerical_calc = Numerical2PGradientCalculator::new(*strategy);
        let gradient = numerical_calc.compute(fcn, &params, trafo, &heuristic_grad);

        // ROOT NegativeG2LineSearch is a seed-only repair: when a diagonal
        // second derivative is non-positive, line-search along that coordinate
        // and recompute all gradients before building the initial covariance.
        let had_negative_g2 = has_negative_g2(&gradient);
        let (params, gradient) = if had_negative_g2 {
            escape_negative_curvature(fcn, params, gradient, trafo, strategy)
        } else {
            (params, gradient)
        };

        // 5. Build V₀. ROOT's ordinary MnSeedGenerator site uses positive
        // fallback for non-positive G2; the post-NegativeG2LineSearch rebuild
        // keeps signed 1/G2 when |G2| is significant.
        let (error, edm) = build_seed_error_and_edm(&gradient, n, eps, had_negative_g2);

        let state = MinimumState::new(params, error, gradient, edm, fcn.num_of_calls());
        // Strategy 2 replaces the diagonal seed with the inverse Hessian.
        let state = hesse_seed(fcn, state, trafo, strategy);

        MinimumSeed::new(state, trafo.clone())
    }

    /// Generate seed using analytical gradients from user.
    pub fn generate_with_gradient(
        fcn: &dyn FCNGradient,
        trafo: &MnUserTransformation,
        strategy: &MnStrategy,
    ) -> MinimumSeed {
        let n = trafo.variable_parameters();
        let eps = trafo.precision().eps();

        // 1. Get initial internal parameter values
        let int_values = trafo.initial_internal_values();
        let int_vec = DVector::from_vec(int_values.clone());

        // 2. Evaluate FCN at starting point (needed for EDM calculation)
        let fval = fcn.value(&trafo.transform(&int_values));
        let params = MinimumParameters::new(int_vec, fval);

        // The calculator stores a user g2 only when no Hessian exists.
        // The seed stores the Hessian diagonal when a Hessian exists.
        // The covariance diagonal is 1/g2.
        let gradient = hessian_diagonal_g2(
            fcn,
            trafo,
            &params,
            AnalyticalGradientCalculator::compute(fcn, trafo, &params),
        );

        // A user g2 or a Hessian diagonal can be non-positive.
        // The seed then runs the negative-curvature escape.
        let had_negative_g2 = has_negative_g2(&gradient);
        let (params, gradient, nfcn) = if had_negative_g2 {
            escape_negative_curvature_analytical(fcn, params, gradient, trafo)
        } else {
            (params, gradient, 1)
        };

        // 4. Build V₀ using the same ordinary-vs-post-NG2LS split as the
        // numerical seed path.
        let (error, edm) = build_seed_error_and_edm(&gradient, n, eps, had_negative_g2);

        // The analytical seed counts the initial evaluation. A negative-curvature
        // escape adds its own evaluations to `nfcn`.
        let state = MinimumState::new(params, error, gradient, edm, nfcn);
        // Strategy 2 replaces the diagonal seed with the inverse Hessian.
        let state = if strategy.strategy() == 2 {
            let mn_fcn = MnFcn::with_num_calls(fcn, trafo, state.nfcn());
            hesse_seed(&mn_fcn, state, trafo, strategy)
        } else {
            state
        };

        MinimumSeed::new(state, trafo.clone())
    }

    pub fn call_with_analytical_gradient_calculator(
        fcn: &dyn FCNGradient,
        trafo: &MnUserTransformation,
        strategy: &MnStrategy,
    ) -> MinimumSeed {
        Self::generate_with_gradient(fcn, trafo, strategy)
    }
}

/// Store the Hessian diagonal when the function has a Hessian.
///
/// The calculator leaves g2 unset in that case.
/// The seed then stores dxdi squared times that diagonal.
/// A user g2 stays only when the function has no Hessian.
/// A missing or short Hessian leaves the gradient unchanged.
fn hessian_diagonal_g2(
    fcn: &dyn FCNGradient,
    trafo: &MnUserTransformation,
    params: &MinimumParameters,
    gradient: FunctionGradient,
) -> FunctionGradient {
    if !fcn.has_hessian() {
        return gradient;
    }
    let external = trafo.transform(params.vec().as_slice());
    let Some(ext_g2) = AnalyticalGradientCalculator::hessian_diagonal(fcn, &external) else {
        return gradient;
    };
    store_external_diagonal(trafo, params, gradient, &ext_g2)
}

/// Store dxdi squared times an external diagonal.
///
/// A bounded parameter adds no second derivative of the transform.
/// A short diagonal leaves the gradient unchanged.
fn store_external_diagonal(
    trafo: &MnUserTransformation,
    params: &MinimumParameters,
    gradient: FunctionGradient,
    ext_g2: &[f64],
) -> FunctionGradient {
    let n = gradient.grad().len();
    let mut g2 = gradient.g2().clone();
    for i in 0..n {
        let ext_idx = trafo.ext_of_int(i);
        if ext_idx >= ext_g2.len() {
            return gradient;
        }
        let dxdi = trafo.dint2ext(ext_idx, params.vec()[i]);
        g2[i] = dxdi * dxdi * ext_g2[ext_idx];
    }
    let mut out = FunctionGradient::new(gradient.grad().clone(), g2, gradient.gstep().clone());
    out.set_analytical(true);
    out
}

fn has_negative_g2(gradient: &FunctionGradient) -> bool {
    gradient.g2().iter().any(|&g2| g2 <= 0.0)
}

/// Replace a strategy-2 seed with the Hesse covariance.
///
/// Other strategies keep the diagonal seed. Hesse uses the call budget
/// `200 + 100 n + 5 n^2` when the caller passes no limit.
fn hesse_seed(
    fcn: &MnFcn,
    state: MinimumState,
    trafo: &MnUserTransformation,
    strategy: &MnStrategy,
) -> MinimumState {
    if strategy.strategy() != 2 {
        return state;
    }
    let n = state.parameters().vec().len();
    let maxcalls = 200 + 100 * n + 5 * n * n;
    hesse_calculator::calculate(fcn, &state, trafo, strategy, maxcalls).state
}

fn build_seed_error_and_edm(
    gradient: &FunctionGradient,
    n: usize,
    eps: f64,
    signed_negative_g2_rebuild: bool,
) -> (MinimumError, f64) {
    let mut v0 = DMatrix::zeros(n, n);
    for i in 0..n {
        let g2i = gradient.g2()[i];
        v0[(i, i)] = if signed_negative_g2_rebuild {
            if g2i.abs() > eps { 1.0 / g2i } else { 1.0 }
        } else if g2i > eps {
            1.0 / g2i
        } else {
            1.0
        };
    }

    let mut error = MinimumError::new(v0, 1.0);
    let edm = {
        let g = gradient.grad();
        let e = error.matrix();
        0.5 * g.dot(&(e * g))
    };

    if signed_negative_g2_rebuild && edm < 0.0 {
        error.set_invert_failed(true);
    }

    (error, edm)
}

fn is_at_zero_jacobian_limit(
    trafo: &MnUserTransformation,
    params: &MinimumParameters,
    int_idx: usize,
) -> bool {
    let ext_idx = trafo.ext_of_int(int_idx);
    let parameter = trafo.parameter(ext_idx);
    (parameter.has_lower_limit() || parameter.has_upper_limit())
        && trafo.dint2ext(ext_idx, params.vec()[int_idx]).abs() < trafo.precision().eps2()
}

fn escape_negative_curvature(
    fcn: &MnFcn,
    params: MinimumParameters,
    gradient: FunctionGradient,
    trafo: &MnUserTransformation,
    strategy: &MnStrategy,
) -> (MinimumParameters, FunctionGradient) {
    let mut recompute_gradient = |params: &MinimumParameters, previous: &FunctionGradient| {
        Numerical2PGradientCalculator::new(*strategy).compute(fcn, params, trafo, previous)
    };
    escape_negative_curvature_with(fcn, params, gradient, trafo, true, &mut recompute_gradient)
}

fn escape_negative_curvature_analytical(
    fcn: &dyn FCNGradient,
    params: MinimumParameters,
    gradient: FunctionGradient,
    trafo: &MnUserTransformation,
) -> (MinimumParameters, FunctionGradient, usize) {
    let mn_fcn = MnFcn::new(fcn, trafo);
    // After a step the escape prefers a finite user g2.
    // The Hessian diagonal is used when that g2 is absent.
    let mut recompute_gradient = |params: &MinimumParameters, _previous: &FunctionGradient| {
        let grad = AnalyticalGradientCalculator::compute(fcn, trafo, params);
        let external = trafo.transform(params.vec().as_slice());
        let Some(ext_g2) =
            AnalyticalGradientCalculator::external_second_derivative(fcn, trafo, &external)
        else {
            return grad;
        };
        store_external_diagonal(trafo, params, grad, &ext_g2)
    };
    let (params, gradient) = escape_negative_curvature_with(
        &mn_fcn,
        params,
        gradient,
        trafo,
        false,
        &mut recompute_gradient,
    );
    (params, gradient, 1 + mn_fcn.num_of_calls())
}

fn escape_negative_curvature_with(
    fcn: &MnFcn,
    mut params: MinimumParameters,
    mut gradient: FunctionGradient,
    trafo: &MnUserTransformation,
    has_gstep: bool,
    recompute_gradient: &mut impl FnMut(&MinimumParameters, &FunctionGradient) -> FunctionGradient,
) -> (MinimumParameters, FunctionGradient) {
    let n = gradient.g2().len();
    let mut iter = 0usize;

    // Mirrors ROOT NegativeG2LineSearch.cxx:62-123 (v6-36-08): repair one
    // offending coordinate per pass, then recompute the full gradient before
    // scanning again. Rust permits one extra pass to preserve the established
    // start-at-limit escape regression at zero-Jacobian transform singularities.
    loop {
        let mut iterate = false;

        for i in 0..n {
            if gradient.g2()[i] <= 0.0 {
                let mut step = DVector::zeros(n);
                let gstep = if has_gstep { gradient.gstep()[i] } else { 1.0 };
                let zero_grad = gradient.grad()[i].abs() < trafo.precision().eps();
                let zero_g2 = gradient.g2()[i].abs() < trafo.precision().eps();
                let zero_gstep = gstep.abs() < trafo.precision().eps();
                if zero_grad && zero_g2 && zero_gstep {
                    continue;
                }

                let at_limit = is_at_zero_jacobian_limit(trafo, &params, i);
                step[i] = if zero_grad && at_limit {
                    // Waived branch: ROOT skips the exact zero-gradient,
                    // zero-Jacobian limit singularity and therefore defines no
                    // direction there. Keep Rust's documented start-at-limit
                    // escape direction; use a finite internal probe so the seed
                    // has enough curvature information not to converge at the
                    // transform singularity.
                    gstep.abs().max(1.0)
                } else if zero_grad && zero_g2 {
                    // Same ROOT-skipped branch away from a detected transform
                    // singularity: preserve the established positive probe.
                    gstep
                } else if gradient.grad()[i] < 0.0 {
                    gstep
                } else {
                    -gstep
                };

                let gdel = step[i] * gradient.grad()[i];
                let ls = mn_linesearch(fcn, &params, &step, gdel, trafo.precision());

                let actual_step = ls.x * &step;
                params =
                    MinimumParameters::with_step(params.vec() + &actual_step, actual_step, ls.y);
                gradient = recompute_gradient(&params, &gradient);

                iterate = true;
                break;
            }
        }

        if !(iter < 2 * n + 1 && iterate) {
            break;
        }
        iter += 1;
    }

    (params, gradient)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fcn::FCN;
    use crate::parameter::MinuitParameter;

    /// f(x, y) = (x^2 - 1)^2 + 4 y^2 with a user Hessian and no g2.
    /// The Hessian diagonal is 12 x^2 - 4 and 8. It is negative at x = 0.1.
    struct DoubleWell;

    impl FCN for DoubleWell {
        fn value(&self, p: &[f64]) -> f64 {
            (p[0] * p[0] - 1.0).powi(2) + 4.0 * p[1] * p[1]
        }

        fn has_hessian(&self) -> bool {
            true
        }

        /// Packed upper triangle: h00, h01, h11.
        fn hessian(&self, p: &[f64]) -> Vec<f64> {
            vec![12.0 * p[0] * p[0] - 4.0, 0.0, 8.0]
        }
    }

    impl FCNGradient for DoubleWell {
        fn gradient(&self, p: &[f64]) -> Vec<f64> {
            vec![4.0 * p[0] * (p[0] * p[0] - 1.0), 8.0 * p[1]]
        }
    }

    /// A positive Hessian diagonal gives the seed g2 `dxdi^2 * H_ii`.
    /// For the sine transform of a parameter in (lo, up), `dxdi^2 = (x - lo) * (up - x)`.
    #[test]
    fn seed_g2_is_hessian_diagonal_times_jacobian_squared() {
        let (lo, up, x0) = (0.5, 3.0, 1.5);
        let trafo = MnUserTransformation::new(vec![
            MinuitParameter::with_limits(0, "x", x0, 0.1, lo, up),
            MinuitParameter::new(1, "y", 0.5, 0.1),
        ]);

        let seed =
            MigradSeedGenerator::generate_with_gradient(&DoubleWell, &trafo, &MnStrategy::new(1));

        let g2 = seed.state().gradient().g2();
        let expected_x = (x0 - lo) * (up - x0) * (12.0 * x0 * x0 - 4.0);
        // The round trip x -> asin -> cos loses a few ulp of 34.5 (ulp 7.1e-15).
        // The bound 1e-10 is about 14000 ulp, far above that loss.
        assert!(
            (g2[0] - expected_x).abs() < 1e-10,
            "g2[0] = {}, expected {}",
            g2[0],
            expected_x
        );
        // y is unbounded, so dxdi is exactly 1 and g2[1] is exactly 8. The bound
        // 1e-12 is about 560 ulp of 8 (ulp 1.8e-15).
        assert!((g2[1] - 8.0).abs() < 1e-12, "g2[1] = {}", g2[1]);
    }

    /// A non-positive Hessian diagonal at the start sends the seed through the
    /// negative-curvature escape. The recomputed gradient keeps the Hessian g2.
    #[test]
    fn seed_escape_keeps_hessian_g2() {
        let trafo = MnUserTransformation::new(vec![
            MinuitParameter::new(0, "x", 0.1, 0.1),
            MinuitParameter::new(1, "y", 0.5, 0.1),
        ]);

        let seed =
            MigradSeedGenerator::generate_with_gradient(&DoubleWell, &trafo, &MnStrategy::new(1));

        let x_end = seed.state().parameters().vec()[0];
        assert!(
            (x_end - 0.1).abs() > 0.5,
            "the escape did not move x from 0.1: x = {x_end}"
        );
        let g2 = seed.state().gradient().g2();
        let expected_x = 12.0 * x_end * x_end - 4.0;
        // Both sides evaluate the same expression of x_end, so they agree to a few
        // ulp of 7.7 (ulp 8.9e-16). The bound 1e-10 is far above that.
        assert!(
            (g2[0] - expected_x).abs() < 1e-10,
            "g2[0] = {}, Hessian diagonal {} (the heuristic is 2/0.1^2 = 200)",
            g2[0],
            expected_x
        );
        // y is unbounded, so dxdi is exactly 1 and g2[1] is exactly 8. The bound
        // 1e-12 is about 560 ulp of 8 (ulp 1.8e-15).
        assert!((g2[1] - 8.0).abs() < 1e-12, "g2[1] = {}", g2[1]);
    }
}
