//! MnFunctionCross: iterative crossing-point finder.
//!
//! Finds the point where the function value equals `fmin + Up`, by iteratively
//! running Migrad with the scanned parameter fixed and using parabolic
//! interpolation to converge on the crossing.

use crate::fcn::FCN;
use crate::migrad::MnMigrad;
use crate::minimum::FunctionMinimum;
use crate::parabola::{MnParabolaPoint, from_3_points};
use crate::parameter::MinuitParameter;
use crate::strategy::MnStrategy;

use super::cross::MnCross;

/// Find where F(par) = Fmin + Up along one parameter direction.
///
/// # Parameters
/// - `fcn`: the user function
/// - `minimum`: the current minimum
/// - `par`: external parameter index being scanned
/// - `pmid`: midpoint parameter value (starting point for scan)
/// - `pdir`: scan direction magnitude
/// - `tlr`: tolerance for convergence (default 0.1)
/// - `maxcalls`: maximum function calls
/// - `strategy`: minimization strategy
#[allow(clippy::too_many_arguments)]
pub fn find_crossing(
    fcn: &dyn FCN,
    minimum: &FunctionMinimum,
    par: usize,
    pmid: f64,
    pdir: f64,
    tlr: f64,
    maxcalls: usize,
    strategy: &MnStrategy,
) -> MnCross {
    let up = minimum.up();
    let fmin = minimum.fval();
    let _nvar = minimum.n_variable_params();

    // Use lower strategy for internal Migrad calls
    let mgr_strategy = if strategy.strategy() > 0 {
        MnStrategy::new(strategy.strategy() - 1)
    } else {
        MnStrategy::new(0)
    };
    let mgr_tlr = 0.5 * tlr;

    let npar = minimum.user_state().len();

    // Tolerances
    let tlf = tlr * up; // function tolerance
    let tla = tlr; // parameter tolerance
    let fmin_delta = 0.01 * up;

    // --- Phase 1: Check limits ---
    let p = minimum.user_state().parameter(par);
    let limset = p.has_lower_limit() || p.has_upper_limit() || p.has_limits();
    if limset && npar == 1 {
        // Single parameter at limit — can't cross
        if pdir > 0.0 && p.has_upper_limit() && pmid >= p.upper_limit() {
            return MnCross::limit_reached(0);
        }
        if pdir < 0.0 && p.has_lower_limit() && pmid <= p.lower_limit() {
            return MnCross::limit_reached(0);
        }
    }

    // --- Phase 2: First Migrad at pmid ---
    let migrad_result = run_migrad_fixed(fcn, minimum, par, pmid, &mgr_strategy, mgr_tlr, maxcalls);

    let mut nfcn_total = migrad_result.nfcn();
    if !migrad_result.is_valid() {
        return MnCross::invalid(nfcn_total);
    }

    // Check if we found a new minimum
    if is_new_minimum(migrad_result.fval(), fmin, fmin_delta) {
        let state = migrad_result.user_state().clone();
        return MnCross::new_minimum_found(state, nfcn_total);
    }

    let f0 = migrad_result.fval();
    let a0 = 0.0_f64; // relative position along scan direction

    // --- Phase 3: Heuristic step ---
    let aopt = if (f0 - fmin).abs() < up * 0.01 {
        // Very close to fmin — take a unit step
        1.0
    } else {
        let ratio = up / (f0 - fmin);
        if ratio > 0.0 {
            (ratio.sqrt() - 1.0).clamp(-0.5, 1.0)
        } else {
            1.0
        }
    };

    // --- Phase 4: Second Migrad ---
    let p1 = pmid + aopt * pdir;
    let migrad2 = run_migrad_fixed(fcn, minimum, par, p1, &mgr_strategy, mgr_tlr, maxcalls);
    nfcn_total += migrad2.nfcn();

    if !migrad2.is_valid() {
        return MnCross::invalid(nfcn_total);
    }

    if is_new_minimum(migrad2.fval(), fmin, fmin_delta) {
        let state = migrad2.user_state().clone();
        return MnCross::new_minimum_found(state, nfcn_total);
    }

    let f1 = migrad2.fval();
    let a1 = aopt;

    // --- Phase 5: Ensure positive slope ---
    let f_left = f0;
    let a_left = a0;
    let mut f_right = f1;
    let mut a_right = a1;
    // dfda = (f1 - f0) / (a1 - a0)
    let mut dfda = if (a1 - a0).abs() > 1e-15 {
        (f1 - f0) / (a1 - a0)
    } else {
        0.0
    };

    // If slope is negative, we need to go further
    let mut maxiter_slope = 15;
    while dfda < 0.0 && maxiter_slope > 0 {
        maxiter_slope -= 1;
        a_right += 0.2;
        let p_try = pmid + a_right * pdir;

        // Check limits
        if limset && at_limit(pdir, p, p_try) {
            return MnCross::limit_reached(nfcn_total);
        }

        let mgr = run_migrad_fixed(fcn, minimum, par, p_try, &mgr_strategy, mgr_tlr, maxcalls);
        nfcn_total += mgr.nfcn();

        if !mgr.is_valid() {
            return MnCross::invalid(nfcn_total);
        }
        if is_new_minimum(mgr.fval(), fmin, fmin_delta) {
            let state = mgr.user_state().clone();
            return MnCross::new_minimum_found(state, nfcn_total);
        }

        f_right = mgr.fval();
        dfda = (f_right - f_left) / (a_right - a_left);
    }

    if dfda < 0.0 {
        return MnCross::invalid(nfcn_total);
    }

    // --- Phase 6: Linear extrapolation to crossing ---
    // We want f(a) = fmin + up
    // Linear: a_cross = a_left + (fmin + up - f_left) / dfda
    let mut a_cross = a_left + (fmin + up - f_left) / dfda;

    // Evaluate
    let p_cross = pmid + a_cross * pdir;
    let mgr_cross = run_migrad_fixed(fcn, minimum, par, p_cross, &mgr_strategy, mgr_tlr, maxcalls);
    nfcn_total += mgr_cross.nfcn();

    if !mgr_cross.is_valid() {
        return MnCross::invalid(nfcn_total);
    }
    if is_new_minimum(mgr_cross.fval(), fmin, fmin_delta) {
        let state = mgr_cross.user_state().clone();
        return MnCross::new_minimum_found(state, nfcn_total);
    }

    let f_cross = mgr_cross.fval();

    // Check convergence
    let adist = (a_cross - a_right).abs();
    let fdist = (f_cross - fmin - up).abs();
    let tla_scaled = if aopt.abs() > 1.0 {
        tla * aopt.abs()
    } else {
        tla
    };

    if adist < tla_scaled && fdist < tlf {
        let state = mgr_cross.user_state().clone();
        return MnCross::valid(a_cross, state, nfcn_total);
    }

    // --- Phase 7: Parabolic convergence ---
    // We have 3 points: (a_left, f_left), (a_right, f_right), (a_cross, f_cross)
    let mut pts = Vec::from([(a_left, f_left), (a_right, f_right), (a_cross, f_cross)]);
    pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let maxitr = 15;
    for _itr in 0..maxitr {
        if nfcn_total >= maxcalls {
            return MnCross::call_limit_reached(nfcn_total);
        }

        // Sort points by parameter value
        pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        // Fit parabola through 3 points (function value vs parameter value)
        let p1 = MnParabolaPoint::new(pts[0].0, pts[0].1);
        let p2 = MnParabolaPoint::new(pts[1].0, pts[1].1);
        let p3 = MnParabolaPoint::new(pts[2].0, pts[2].1);

        let parab = from_3_points(p1, p2, p3);

        // Where does the parabola equal fmin + up?
        // a*x^2 + b*x + c = fmin + up
        // a*x^2 + b*x + (c - fmin - up) = 0
        let target = fmin + up;
        let disc = parab.b() * parab.b() - 4.0 * parab.a() * (parab.c() - target);

        if disc < 0.0 || parab.a().abs() < 1e-15 {
            // Parabola doesn't cross target — fall back to linear
            let slope = (pts[2].1 - pts[0].1) / (pts[2].0 - pts[0].0);
            if slope.abs() < 1e-15 {
                return MnCross::invalid(nfcn_total);
            }
            a_cross = pts[0].0 + (target - pts[0].1) / slope;
        } else {
            let sqrt_disc = disc.sqrt();
            let root1 = (-parab.b() + sqrt_disc) / (2.0 * parab.a());
            let root2 = (-parab.b() - sqrt_disc) / (2.0 * parab.a());

            // Choose root closest to the bracket
            let mid_a = 0.5 * (pts[0].0 + pts[2].0);
            a_cross = if (root1 - mid_a).abs() < (root2 - mid_a).abs() {
                root1
            } else {
                root2
            };
        }

        // Clamp to reasonable range (slightly beyond bracket)
        let smalla = 0.01 * (pts[2].0 - pts[0].0).abs().max(1e-10);
        let a_lo = pts[0].0 - smalla;
        let a_hi = pts[2].0 + smalla;
        a_cross = a_cross.clamp(a_lo, a_hi);

        // Evaluate at new point
        let p_try = pmid + a_cross * pdir;

        // Check limits
        if limset && at_limit(pdir, p, p_try) {
            return MnCross::limit_reached(nfcn_total);
        }

        let mgr = run_migrad_fixed(fcn, minimum, par, p_try, &mgr_strategy, mgr_tlr, maxcalls);
        nfcn_total += mgr.nfcn();

        if !mgr.is_valid() {
            return MnCross::invalid(nfcn_total);
        }
        if is_new_minimum(mgr.fval(), fmin, fmin_delta) {
            let state = mgr.user_state().clone();
            return MnCross::new_minimum_found(state, nfcn_total);
        }

        let f_new = mgr.fval();

        // Check convergence
        let adist = (a_cross - pts[1].0).abs();
        let fdist = (f_new - target).abs();
        let tla_scaled = if aopt.abs() > 1.0 {
            tla * aopt.abs()
        } else {
            tla
        };

        if adist < tla_scaled && fdist < tlf {
            let state = mgr.user_state().clone();
            return MnCross::valid(a_cross, state, nfcn_total);
        }

        // Replace the farthest-from-target point
        let new_pt = (a_cross, f_new);
        // Find which existing point to replace: the one whose f is farthest from target
        let mut worst_idx = 0;
        let mut worst_dist = (pts[0].1 - target).abs();
        for (idx, pt) in pts.iter().enumerate().skip(1) {
            let d = (pt.1 - target).abs();
            if d > worst_dist {
                worst_dist = d;
                worst_idx = idx;
            }
        }
        pts[worst_idx] = new_pt;
    }

    // Didn't converge after maxitr — return best estimate
    MnCross::invalid(nfcn_total)
}

/// Run Migrad with one parameter fixed at a given value.
fn run_migrad_fixed(
    fcn: &dyn FCN,
    minimum: &FunctionMinimum,
    fix_par: usize,
    fix_val: f64,
    strategy: &MnStrategy,
    tolerance: f64,
    maxcalls: usize,
) -> FunctionMinimum {
    let user_state = minimum.user_state();
    let nparams = user_state.len();

    let mut builder = MnMigrad::new()
        .with_strategy(strategy.strategy())
        .tolerance(tolerance)
        .max_fcn(maxcalls);

    // Add all parameters from the minimum, with the scan parameter fixed
    for i in 0..nparams {
        let p = user_state.parameter(i);
        let val = if i == fix_par { fix_val } else { p.value() };
        builder = add_parameter_to_builder(builder, p, val);
        if i != fix_par && p.is_fixed() && !p.is_const() {
            builder = builder.fix(i);
        }
    }

    // Fix the scan parameter.
    builder = builder.fix(fix_par);

    builder.minimize(fcn)
}

fn at_limit(pdir: f64, p: &MinuitParameter, val: f64) -> bool {
    (pdir > 0.0 && p.has_upper_limit() && val > p.upper_limit())
        || (pdir < 0.0 && p.has_lower_limit() && val < p.lower_limit())
}

fn is_new_minimum(fval: f64, fmin: f64, fmin_delta: f64) -> bool {
    fval < fmin - fmin_delta
}

fn add_parameter_to_builder(mut builder: MnMigrad, p: &MinuitParameter, val: f64) -> MnMigrad {
    let err = p.error();
    if p.has_limits() {
        builder = builder.add_limited(p.name(), val, err, p.lower_limit(), p.upper_limit());
    } else if p.has_lower_limit() {
        builder = builder.add_lower_limited(p.name(), val, err, p.lower_limit());
    } else if p.has_upper_limit() {
        builder = builder.add_upper_limited(p.name(), val, err, p.upper_limit());
    } else if p.is_const() {
        builder = builder.add_const(p.name(), val);
    } else {
        builder = builder.add(p.name(), val, err.max(1e-10));
    }
    builder
}
