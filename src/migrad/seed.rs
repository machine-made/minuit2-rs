//! Migrad seed generator.
//!
//! Creates the initial `MinimumSeed` by evaluating the FCN, computing a
//! numerical gradient, and building `V0 = diag(1/g2_i)`.

use nalgebra::{DMatrix, DVector};

use crate::fcn::FCNGradient;
use crate::gradient::{
    AnalyticalGradientCalculator, InitialGradientCalculator, Numerical2PGradientCalculator,
};
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

        MinimumSeed::new(state, trafo.clone())
    }

    /// Generate seed using analytical gradients from user.
    pub fn generate_with_gradient(
        fcn: &dyn FCNGradient,
        trafo: &MnUserTransformation,
        _strategy: &MnStrategy,
    ) -> MinimumSeed {
        let n = trafo.variable_parameters();
        let eps = trafo.precision().eps();

        // 1. Get initial internal parameter values
        let int_values = trafo.initial_internal_values();
        let int_vec = DVector::from_vec(int_values.clone());

        // 2. Evaluate FCN at starting point (needed for EDM calculation)
        let fval = fcn.value(&trafo.transform(&int_values));
        let params = MinimumParameters::new(int_vec, fval);

        // 3. Compute analytical gradient (user-provided, with g2/gstep heuristics).
        // A user Hessian replaces that heuristic. MnSeedGenerator.cxx:168-180
        // stores the Hessian diagonal as g2 and sets V to 1/g2.
        let gradient = hessian_diagonal_g2(
            fcn,
            trafo,
            &params,
            AnalyticalGradientCalculator::compute(fcn, trafo, &params),
        );

        // ROOT also routes analytical-gradient seeds through NegativeG2LineSearch
        // when the gradient calculator supplies non-positive G2 values. The
        // heuristic g2 of AnalyticalGradientCalculator is positive. A user
        // Hessian can give a non-positive diagonal, so the guard is reachable
        // with a user Hessian.
        let had_negative_g2 = has_negative_g2(&gradient);
        let (params, gradient, nfcn) = if had_negative_g2 {
            escape_negative_curvature_analytical(fcn, params, gradient, trafo)
        } else {
            (params, gradient, 1)
        };

        // 4. Build V₀ using the same ordinary-vs-post-NG2LS split as the
        // numerical seed path.
        let (error, edm) = build_seed_error_and_edm(&gradient, n, eps, had_negative_g2);

        // Note: analytical-gradient seeding counts the initial FCN evaluation;
        // any NegativeG2 line-search evaluations are included in `nfcn`.
        let state = MinimumState::new(params, error, gradient, edm, nfcn);

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

/// Replace the heuristic g2 of `gradient` by the user Hessian diagonal.
///
/// The function keeps the heuristic g2 when it has no Hessian or returns a
/// Hessian that does not cover every external index. It gives no error and no
/// log in the last case. ROOT indexes the Hessian by `iext * next + jext` with
/// no check (`AnalyticalGradientCalculator.cxx:89`).
///
/// It also keeps the heuristic g2 when the function has a g2 of its own. That
/// is an open gap: ROOT reads the g2 of the function, and
/// `AnalyticalGradientCalculator::compute` does not.
fn hessian_diagonal_g2(
    fcn: &dyn FCNGradient,
    trafo: &MnUserTransformation,
    params: &MinimumParameters,
    gradient: FunctionGradient,
) -> FunctionGradient {
    if fcn.has_g2() || !fcn.has_hessian() {
        return gradient;
    }
    let external = trafo.transform(params.vec().as_slice());
    let Some(ext_g2) = AnalyticalGradientCalculator::g2(fcn, &external) else {
        return gradient;
    };
    let n = gradient.grad().len();
    let mut g2 = gradient.g2().clone();
    for i in 0..n {
        let ext_idx = trafo.ext_of_int(i);
        if ext_idx >= ext_g2.len() {
            return gradient;
        }
        let dxdi = trafo.dint2ext(ext_idx, params.vec()[i]);
        // ROOT computes dxdi * H_ext * dxdj with no curvature term
        // d2x/dint2 * grad_ext for a bounded parameter
        // (AnalyticalGradientCalculator.cxx:83-91). The port does the same.
        g2[i] = dxdi * dxdi * ext_g2[ext_idx];
    }
    let mut out = FunctionGradient::new(gradient.grad().clone(), g2, gradient.gstep().clone());
    out.set_analytical(true);
    out
}

fn has_negative_g2(gradient: &FunctionGradient) -> bool {
    gradient.g2().iter().any(|&g2| g2 <= 0.0)
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
    // The recompute keeps the user Hessian diagonal as g2, as the first gradient does.
    let mut recompute_gradient = |params: &MinimumParameters, _previous: &FunctionGradient| {
        hessian_diagonal_g2(
            fcn,
            trafo,
            params,
            AnalyticalGradientCalculator::compute(fcn, trafo, params),
        )
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
            "the escape did not move x from 0.1: x = {}",
            x_end
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
