//! Variable-metric (Migrad) iteration loop with DFP update.
//!
//! Main quasi-Newton iteration: compute step, line search, update gradient,
//! and apply the DFP rank-2 inverse Hessian update.

use nalgebra::DMatrix;

#[cfg(feature = "trace")]
use std::io::Write;

use crate::fcn::FCNGradient;
use crate::gradient::{AnalyticalGradientCalculator, Numerical2PGradientCalculator};
use crate::hesse::calculator as hesse_calculator;
use crate::linesearch::mn_linesearch;
use crate::minimum::error::{ErrorMatrixStatus, MinimumError};
use crate::minimum::gradient::FunctionGradient;
use crate::minimum::parameters::MinimumParameters;
use crate::minimum::seed::MinimumSeed;
use crate::minimum::state::MinimumState;
use crate::mn_fcn::MnFcn;
use crate::posdef::make_pos_def;
use crate::strategy::MnStrategy;

pub struct VariableMetricBuilder;

#[cfg(feature = "trace")]
fn trace_iteration(
    iter: usize,
    nfcn: usize,
    fval: f64,
    edm: f64,
    lambda: f64,
    grad_norm: f64,
    dcovar: f64,
) {
    let Ok(path) = std::env::var("MINUIT2_RS_TRACE_JSONL") else {
        return;
    };
    if path.is_empty() {
        return;
    }
    let Ok(mut file) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
    else {
        return;
    };
    let _ = writeln!(
        file,
        "{{\"runner\":\"minuit2-rs\",\"iter\":{},\"nfcn\":{},\"fval\":{:.17},\"edm\":{:.17},\"lambda\":{:.17},\"grad_norm\":{:.17},\"dcovar\":{:.17}}}",
        iter, nfcn, fval, edm, lambda, grad_norm, dcovar
    );
}

#[cfg(not(feature = "trace"))]
#[inline]
fn trace_iteration(
    _iter: usize,
    _nfcn: usize,
    _fval: f64,
    _edm: f64,
    _lambda: f64,
    _grad_norm: f64,
    _dcovar: f64,
) {
}

impl VariableMetricBuilder {
    /// Top-level Migrad minimization: run iterations, optionally re-seed on failure.
    ///
    /// Returns the iteration history as a `Vec<MinimumState>`.
    pub fn minimum(
        fcn: &MnFcn,
        seed: &MinimumSeed,
        strategy: &MnStrategy,
        maxfcn: usize,
        edmval: f64,
    ) -> Vec<MinimumState> {
        let grad_calc = Numerical2PGradientCalculator::new(*strategy);
        let next_grad = |p: &MinimumParameters, prev: &FunctionGradient| {
            grad_calc.compute_with_previous(fcn, p, seed.trafo(), prev)
        };
        Self::minimize_with_reseed(fcn, seed, strategy, maxfcn, edmval, next_grad)
    }

    /// Top-level Migrad minimization with analytical gradients.
    pub fn minimum_with_gradient(
        fcn: &MnFcn,
        gradient_fcn: &dyn FCNGradient,
        seed: &MinimumSeed,
        _strategy: &MnStrategy,
        maxfcn: usize,
        edmval: f64,
    ) -> Vec<MinimumState> {
        let next_grad = |p: &MinimumParameters, _prev: &FunctionGradient| {
            AnalyticalGradientCalculator::compute(gradient_fcn, seed.trafo(), p)
        };
        Self::minimize_with_reseed(fcn, seed, _strategy, maxfcn, edmval, next_grad)
    }

    /// Run variable-metric passes, re-seeding from the last state when ROOT's
    /// Hesse verification or the EDM test says more work is needed.
    ///
    /// `next_grad(new_params, prev_grad)` computes the gradient at `new_params`;
    /// the numerical strategy uses `prev_grad` for step-size warm-starting while
    /// the analytical strategy ignores it.
    fn minimize_with_reseed(
        fcn: &MnFcn,
        seed: &MinimumSeed,
        strategy: &MnStrategy,
        maxfcn: usize,
        edmval: f64,
        mut next_grad: impl FnMut(&MinimumParameters, &FunctionGradient) -> FunctionGradient,
    ) -> Vec<MinimumState> {
        let mut states = Self::iterate(fcn, seed, maxfcn, edmval, &mut next_grad);

        // ROOT Minuit2 verifies a nominally converged variable-metric result
        // with MnHesse for strategy >= 2, and for strategy 1 when the updated
        // covariance is still uncertain (`Dcovar() > 0.05`). If Hesse raises
        // EDM above the requested tolerance and the value is not below the
        // machine-accuracy floor, MIGRAD continues from the Hesse state with a
        // 30% larger call budget (ROOT v6-36-08,
        // math/minuit2/src/VariableMetricBuilder.cxx:112-151). Hahn1's first
        // pass stops in a local basin unless this Hesse-verified continuation
        // is performed.
        let should_hesse = |state: &MinimumState| {
            strategy.strategy() >= 2 || (strategy.strategy() == 1 && state.error().dcovar() > 0.05)
        };

        let maxfcn2 = (maxfcn as f64 * 1.3) as usize;
        for _pass in 0..5 {
            let Some(last) = states.last() else {
                return states;
            };

            let mut must_continue = last.edm() >= edmval;
            if should_hesse(last) {
                let mut hesse_strategy = *strategy;
                hesse_strategy.set_hessian_force_pos_def(1);
                let hesse =
                    hesse_calculator::calculate(fcn, last, seed.trafo(), &hesse_strategy, maxfcn);
                let hesse_state = hesse.state;
                let hesse_is_valid = hesse_state.is_valid() && hesse_state.error().is_valid();
                let machine_limit = (seed.precision().eps2() * hesse_state.fval()).abs();
                must_continue = hesse_is_valid
                    && hesse_state.edm() > edmval
                    && hesse_state.edm() >= machine_limit;
                if !hesse_is_valid {
                    // ROOT keeps a failed Hesse result diagnostic, but the
                    // returned minimum is not advanced to an unusable Hesse
                    // point unless the error is available/call-limit/above-EDM
                    // (VariableMetricBuilder.cxx:183-198). Preserve the last
                    // usable parameters while carrying the failed error state
                    // so public validity does not silently fall back to
                    // parameters-only success.
                    states.push(MinimumState::new(
                        last.parameters().clone(),
                        hesse_state.error().clone(),
                        last.gradient().clone(),
                        last.edm(),
                        hesse_state.nfcn(),
                    ));
                    return states;
                }
                states.push(hesse_state);
            }

            if !must_continue {
                return states;
            }

            let remaining = maxfcn2.saturating_sub(fcn.num_of_calls());
            if remaining < 10 {
                return states;
            }

            let last = states.last().unwrap_or_else(|| seed.state());
            let seed2 = MinimumSeed::new(
                MinimumState::new(
                    last.parameters().clone(),
                    last.error().clone(),
                    last.gradient().clone(),
                    last.edm(),
                    last.nfcn(),
                ),
                seed.trafo().clone(),
            );

            let states2 = Self::iterate(fcn, &seed2, maxfcn2, edmval, &mut next_grad);
            if states2.is_empty() {
                return states;
            }
            states.extend(states2);
        }
        states
    }

    /// Core quasi-Newton iteration: Newton step → pos-def fallback → line search
    /// → gradient update (via `next_grad`) → DFP update → EDM check.
    fn iterate(
        fcn: &MnFcn,
        seed: &MinimumSeed,
        maxfcn: usize,
        edmval: f64,
        next_grad: &mut impl FnMut(&MinimumParameters, &FunctionGradient) -> FunctionGradient,
    ) -> Vec<MinimumState> {
        let n = seed.n_variable_params();
        let prec = seed.precision();

        let mut params = seed.parameters().clone();
        let mut error = seed.error().clone();
        let mut gradient = seed.gradient().clone();
        let mut edm = seed.edm();

        let mut states = Vec::new();
        let mut iter = 0usize;

        loop {
            // 1. Newton step: step = -V * grad
            let v = error.matrix();
            let g = gradient.grad();
            let step = -(v * g);

            // 2. Ensure descent direction: gdel = step · grad must be negative.
            //    If not, make V positive-definite; if still bad, fall back to
            //    steepest descent with the identity matrix.
            let mut gdel = step.dot(g);

            let (current_step, mut current_error) = if gdel > 0.0 {
                let (v_fixed, _was_modified) = make_pos_def(v, prec);
                let mut err_fixed = MinimumError::new(v_fixed.clone(), error.dcovar());
                err_fixed.set_made_pos_def(true);
                let step_fixed = -(&v_fixed * g);
                gdel = step_fixed.dot(g);

                if gdel > 0.0 {
                    let step_sd = -g.clone();
                    gdel = step_sd.dot(g);
                    let err_sd = MinimumError::new(DMatrix::identity(n, n), 1.0);
                    (step_sd, err_sd)
                } else {
                    (step_fixed, err_fixed)
                }
            } else {
                (step, error.clone())
            };

            // 3. Line search: parabolic interpolation along current_step
            let ls_result = mn_linesearch(fcn, &params, &current_step, gdel, prec);
            let lambda = ls_result.x;
            let f_new = ls_result.y;

            // 4. No-improvement guard: stop if we are at machine-precision plateau
            if (f_new - params.fval()).abs() <= params.fval().abs() * prec.eps() {
                let new_params = MinimumParameters::with_step(
                    params.vec() + lambda * &current_step,
                    lambda * &current_step,
                    f_new,
                );
                trace_iteration(
                    iter,
                    fcn.num_of_calls(),
                    f_new,
                    edm,
                    lambda,
                    gradient.grad().norm(),
                    current_error.dcovar(),
                );
                states.push(MinimumState::new(
                    new_params,
                    current_error.clone(),
                    gradient.clone(),
                    edm,
                    fcn.num_of_calls(),
                ));
                break;
            }

            // 5. Advance parameters
            let p_new = params.vec() + lambda * &current_step;
            let new_params = MinimumParameters::with_step(p_new, lambda * &current_step, f_new);

            // 6. Compute gradient at the new point
            let new_gradient = next_grad(&new_params, &gradient);

            // 7. Estimate EDM with the pre-update covariance, matching ROOT
            // Minuit2's VariableMetricBuilder.cxx: after `FunctionGradient g = gc(...)`,
            // it calls `edm = Estimator().Estimate(g, s0.Error())` before
            // `MinimumError e = ErrorUpdator().Update(s0, p, g)` (ROOT v6-36-08,
            // math/minuit2/src/VariableMetricBuilder.cxx:291-312).  The DFP-updated
            // covariance is stored in the state, but it is not used for this
            // iteration's EDM estimate; only the loop's local convergence check is
            // corrected later by `(1 + 3 * e.Dcovar())`.
            let new_g = new_gradient.grad();
            edm = 0.5 * new_g.dot(&(current_error.matrix() * new_g));
            if edm.is_nan() {
                states.push(MinimumState::new(
                    params.clone(),
                    current_error.clone(),
                    gradient.clone(),
                    edm,
                    fcn.num_of_calls(),
                ));
                break;
            }
            if edm < 0.0 {
                let (v_fixed, _was_modified) = make_pos_def(current_error.matrix(), prec);
                let mut err_fixed = MinimumError::new(v_fixed, current_error.dcovar());
                err_fixed.set_made_pos_def(true);
                current_error = err_fixed;
                edm = 0.5 * new_g.dot(&(current_error.matrix() * new_g));
                if edm < 0.0 {
                    states.push(MinimumState::new(
                        params.clone(),
                        current_error.clone(),
                        gradient.clone(),
                        edm,
                        fcn.num_of_calls(),
                    ));
                    break;
                }
            }

            // 8. DFP rank-2 update of the inverse Hessian approximation
            let (v_updated, new_dcovar) = Self::dfp_update(
                &current_error,
                &new_params,
                &params,
                &new_gradient,
                &gradient,
            );

            let mut new_error = MinimumError::new(v_updated, new_dcovar);
            if current_error.status() == ErrorMatrixStatus::MadePositiveDefinite {
                new_error.set_made_pos_def(true);
            }

            trace_iteration(
                iter,
                fcn.num_of_calls(),
                new_params.fval(),
                edm,
                lambda,
                new_gradient.grad().norm(),
                new_dcovar,
            );
            states.push(MinimumState::new(
                new_params.clone(),
                new_error.clone(),
                new_gradient.clone(),
                edm,
                fcn.num_of_calls(),
            ));

            let corrected_edm = edm * (1.0 + 3.0 * new_dcovar);
            if corrected_edm < edmval {
                break;
            }

            if fcn.num_of_calls() >= maxfcn {
                break;
            }

            iter += 1;
            params = new_params;
            error = new_error;
            gradient = new_gradient;
        }

        states
    }

    /// Rank-2 DFP update of the inverse Hessian approximation.
    ///
    /// Returns `(V_new, dcovar)` where `dcovar` measures how much the matrix changed.
    pub fn update(
        error: &MinimumError,
        p_new: &MinimumParameters,
        p_old: &MinimumParameters,
        g_new: &FunctionGradient,
        g_old: &FunctionGradient,
    ) -> (DMatrix<f64>, f64) {
        Self::dfp_update(error, p_new, p_old, g_new, g_old)
    }

    pub fn nrow(error: &MinimumError) -> usize {
        error.matrix().nrows()
    }

    pub fn has_negative_g2(gradient: &FunctionGradient) -> bool {
        gradient.g2().iter().any(|g2| *g2 <= 0.0)
    }

    pub fn add_result(states: &mut Vec<MinimumState>, state: MinimumState) {
        states.push(state);
    }

    fn dfp_update(
        error: &MinimumError,
        p_new: &MinimumParameters,
        p_old: &MinimumParameters,
        g_new: &FunctionGradient,
        g_old: &FunctionGradient,
    ) -> (DMatrix<f64>, f64) {
        let v = error.matrix();

        let dx = p_new.vec() - p_old.vec();
        let dg = g_new.grad() - g_old.grad();

        // delgam = dx · dg  (curvature along the step)
        let delgam = dx.dot(&dg);

        // vg = V * dg;  gvg = dg^T * V * dg
        let vg = v * &dg;
        let gvg = dg.dot(&vg);

        // ROOT Minuit2's DavidonErrorUpdator only skips the update for
        // `delgam == 0` or `gvg <= 0`; for `delgam < 0` it warns but still
        // applies the Davidon formula (ROOT v6-36-08,
        // math/minuit2/src/DavidonErrorUpdator.cxx:42-78).  Preserving that
        // negative-curvature update is important for difficult trajectories
        // such as Hahn1, where dropping it changes the quasi-Newton path.
        if delgam == 0.0 || gvg <= 0.0 {
            return (v.clone(), error.dcovar());
        }

        // Rank-2 update matrix: outer(dx, dx)/delgam - outer(vg, vg)/gvg
        let mut v_upd = &dx * dx.transpose() / delgam - &vg * vg.transpose() / gvg;

        // Optional rank-1 correction when curvature condition is strong
        if delgam > gvg {
            let flnu = &dx / delgam - &vg / gvg;
            // Scale flnu by gvg before the outer product so each element is
            // (gvg * flnu[i]) * flnu[j], matching the original left-to-right
            // associativity exactly (f64 multiplication is not associative).
            v_upd += (&flnu * gvg) * flnu.transpose();
        }

        let v_new = v + &v_upd;

        let sum_upd: f64 = v_upd.iter().map(|x| x.abs()).sum();
        let sum_new: f64 = v_new.iter().map(|x| x.abs()).sum();
        let dcovar = if sum_new > 0.0 {
            0.5 * (error.dcovar() + sum_upd / sum_new)
        } else {
            error.dcovar()
        };

        (v_new, dcovar)
    }
}
