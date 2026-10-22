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

        // 3. Compute analytical gradient (user-provided, with g2/gstep heuristics)
        let gradient = AnalyticalGradientCalculator::compute(fcn, trafo, &params);

        // ROOT also routes analytical-gradient seeds through NegativeG2LineSearch
        // when the gradient calculator supplies non-positive G2 values.  In this
        // Rust implementation AnalyticalGradientCalculator currently synthesizes
        // positive G2 heuristics, so this is normally a no-op; keep the guard for
        // traceability if custom G2 support is added later.
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
    let mut recompute_gradient = |params: &MinimumParameters, _previous: &FunctionGradient| {
        AnalyticalGradientCalculator::compute(fcn, trafo, params)
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
