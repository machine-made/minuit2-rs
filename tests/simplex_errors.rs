//! Simplex errors and covariance against ROOT Minuit2 `SimplexBuilder.cxx` and
//! `MnUserParameterState.cxx`: the parameter error is the final simplex step
//! scaled by `sqrt(up / edm)`, and the result has no covariance.

use std::cell::RefCell;
use std::path::Path;
use std::process::Command;

use minuit2::minimum::error::ErrorMatrixStatus;
use minuit2::{FunctionMinimum, MnSimplex};
use serde_json::Value;

/// Relative tolerance for the comparison with the stored ROOT output. The two
/// runners execute the same algorithm in binary64 and the stored values carry
/// 17 digits. Measured worst relative difference: 0, both errors agree bit for
/// bit. A wrong scaling moves an error by a factor of 2 or more, so 1e-6 still
/// detects the defect and leaves room for a different compiler or math library.
const ROOT_TOLERANCE: f64 = 1e-6;

/// Relative tolerance for the errors computed in these tests from the recorded
/// simplex vertices. These are not a comparison with ROOT. The test repeats the
/// library's floating-point operations (spread of the vertices, times
/// `sqrt(up / edm)`), so the two values agree bit for bit unless the compiler
/// changes the order of the operations. Measured worst relative difference over
/// all unlimited errors below: 0. The tolerance 1e-12 is a few thousand units in
/// the last place and still far below any wrong mapping, which changes an error
/// by more than 1 percent in the tests below.
const RECORDED_TOLERANCE: f64 = 1e-12;

/// Relative tolerance for a limited parameter. The test converts the reported
/// external value back to the internal coordinate before it calls
/// `int2ext_error`, and this round trip changes the last digits. Measured worst
/// relative difference: 5.8e-15, so 1e-10 leaves a margin of about 1.7e4.
const LIMITED_TOLERANCE: f64 = 1e-10;

fn assert_relative(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let difference = (actual - expected).abs();
    assert!(
        difference <= tolerance * expected.abs(),
        "{what}: {actual:?} against {expected:?}, relative difference {}",
        difference / expected.abs()
    );
}

/// The `errors` of the Rust runner follow the stored ROOT output of the
/// `quadratic2_simplex` workload, and neither result has a covariance.
#[test]
fn simplex_errors_match_root() {
    let output = Command::new(env!("CARGO_BIN_EXE_ref_compare_runner"))
        .args(["--workload", "quadratic2_simplex"])
        .output()
        .expect("the runner should start");
    assert!(
        output.status.success(),
        "the runner failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let rust: Value = serde_json::from_slice(&output.stdout).expect("the runner should print JSON");

    let root_path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("reports/verification/raw/ref/quadratic2_simplex.json");
    let root_text = std::fs::read_to_string(&root_path).expect("the ROOT output should exist");
    let root: Value = serde_json::from_str(&root_text).expect("the ROOT output should be JSON");

    let root_errors = root["errors"].as_array().expect("ROOT errors");
    let rust_errors = rust["errors"].as_array().expect("Rust errors");
    assert_eq!(rust_errors.len(), root_errors.len());
    for (i, (rust_error, root_error)) in rust_errors.iter().zip(root_errors).enumerate() {
        assert_relative(
            rust_error.as_f64().expect("Rust error"),
            root_error.as_f64().expect("ROOT error"),
            ROOT_TOLERANCE,
            &format!("error of parameter {i}"),
        );
    }

    assert_eq!(root["has_covariance"], Value::Bool(false));
    assert_eq!(
        rust["has_covariance"],
        Value::Bool(false),
        "Rust reports a covariance: {}",
        rust["covariance"]
    );
    assert!(
        rust["covariance"].is_null(),
        "Rust reports a covariance: {}",
        rust["covariance"]
    );
}

/// `f(x) = (x - 1)^2` from 0 with step 0.5: the error is the spread of the
/// final two simplex vertices times `sqrt(up / edm)`. The function records its
/// calls, so the test finds the vertices without the library's help.
#[test]
fn simplex_error_is_scaled_final_step() {
    let calls = RefCell::new(Vec::<(f64, f64)>::new());
    let f = |p: &[f64]| {
        let value = (p[0] - 1.0).powi(2);
        calls.borrow_mut().push((p[0], value));
        value
    };
    let min = MnSimplex::new().add("x", 0.0, 0.5).minimize(&f);

    let up = min.up();
    let edm = min.edm();
    let best_x = min.params()[0];
    let best_value = min.fval();

    // With one parameter the final simplex is the best vertex and the worst
    // vertex. The edm is the difference of their function values, so the worst
    // vertex is the call whose difference equals the edm bit for bit.
    let calls = calls.borrow();
    let (worst_x, _) = calls
        .iter()
        .rev()
        .find(|(x, value)| {
            x.total_cmp(&best_x).is_ne() && (value - best_value).total_cmp(&edm).is_eq()
        })
        .expect("the worst vertex should be among the calls");
    let expected = (best_x - worst_x).abs() * (up / edm).sqrt();
    assert!(expected > 0.0);

    let errors = min.user_state().errors();
    assert_relative(errors[0], expected, RECORDED_TOLERANCE, "error of x");

    let step = min.state().parameters().step()[0];
    assert_relative(
        errors[0],
        step,
        RECORDED_TOLERANCE,
        "error of x against the final step",
    );

    assert!(!min.user_state().has_covariance());
    assert!(min.user_state().covariance().is_none());
}

/// ROOT's Simplex result holds a default `MinimumError`: no matrix is available
/// and no failure flag is set. No flag may report a failed inversion or Hesse.
#[test]
fn simplex_error_matrix_is_unset_without_failure_flags() {
    let min = MnSimplex::new()
        .add("x", 0.0, 0.5)
        .minimize(&|p: &[f64]| (p[0] - 1.0).powi(2));

    let error = min.state().error();
    assert_eq!(error.status(), ErrorMatrixStatus::NotAvailable);
    assert!(!error.is_valid());
    assert!(!error.is_available());
    assert!(!min.state().has_covariance());
    assert!(!error.invert_failed());
    assert!(!error.is_not_pos_def());
    assert!(!error.hesse_failed());
    assert!(!error.is_made_pos_def());
    assert!(!error.has_reached_call_limit());
}

// The three tests below have no ROOT file to compare with. Their expected errors
// come from the vertices of the final simplex, which the test finds in the
// calls its function recorded. A Simplex run keeps the `n + 1` lowest distinct
// function values it has seen, so the final simplex is the `n + 1` lowest
// recorded points. `final_simplex` checks this against the reported `fval` and
// `edm`, bit for bit, and fails if the identification is wrong.

type Calls = RefCell<Vec<(Vec<f64>, f64)>>;

fn final_simplex(calls: &Calls, vertices: usize, fval: f64, edm: f64) -> Vec<Vec<f64>> {
    let calls = calls.borrow();
    let mut distinct: Vec<&(Vec<f64>, f64)> = Vec::new();
    for call in calls.iter() {
        if !distinct.iter().any(|seen| seen.0 == call.0) {
            distinct.push(call);
        }
    }
    distinct.sort_by(|a, b| a.1.total_cmp(&b.1));
    let simplex = &distinct[..vertices];
    let best = simplex[0].1;
    let worst = simplex[vertices - 1].1;
    assert!(
        best.total_cmp(&fval).is_eq(),
        "best recorded value {best:?} is not the reported fval {fval:?}"
    );
    assert!(
        (worst - best).total_cmp(&edm).is_eq(),
        "recorded spread {:?} is not the reported edm {edm:?}",
        worst - best
    );
    simplex.iter().map(|call| call.0.clone()).collect()
}

/// The error of external parameter `ext`: the spread of the final vertices in
/// the internal coordinate, times `sqrt(up / edm)`, mapped to the external
/// coordinate by `int2ext_error` at the reported value.
fn expected_error(min: &FunctionMinimum, vertices: &[Vec<f64>], ext: usize) -> f64 {
    let trafo = min.seed().trafo();
    let internal: Vec<f64> = vertices
        .iter()
        .map(|vertex| trafo.ext2int(ext, vertex[ext]))
        .collect();
    let high = internal.iter().copied().fold(f64::MIN, f64::max);
    let low = internal.iter().copied().fold(f64::MAX, f64::min);
    let step = (high - low) * (min.up() / min.edm()).sqrt();
    trafo.int2ext_error(ext, trafo.ext2int(ext, min.params()[ext]), step)
}

/// `x` with the lower limit 0 and a free `y`, in `(x - 1)^2 + 4 (y + 2)^2`. The
/// error of `x` is the scaled step mapped by `int2ext_error`, which differs
/// from the step itself. The error of `y` is its scaled step.
#[test]
fn simplex_limited_error_goes_through_the_external_transform() {
    let calls: Calls = RefCell::new(Vec::new());
    let f = |p: &[f64]| {
        let value = (p[0] - 1.0).powi(2) + 4.0 * (p[1] + 2.0).powi(2);
        calls.borrow_mut().push((p.to_vec(), value));
        value
    };
    let min = MnSimplex::new()
        .add_lower_limited("x", 0.5, 0.3, 0.0)
        .add("y", 0.0, 0.3)
        .tolerance(0.001)
        .minimize(&f);

    let vertices = final_simplex(&calls, 3, min.fval(), min.edm());
    let errors = min.user_state().errors();
    let expected_x = expected_error(&min, &vertices, 0);
    let expected_y = expected_error(&min, &vertices, 1);
    assert_relative(errors[0], expected_x, LIMITED_TOLERANCE, "error of x");
    assert_relative(errors[1], expected_y, RECORDED_TOLERANCE, "error of y");

    // The transform matters: the reported error is not the internal step.
    let step = min.state().parameters().step();
    assert_relative(errors[1], step[1], RECORDED_TOLERANCE, "y against its step");
    assert!(
        (errors[0] - step[0]).abs() > 0.01 * step[0],
        "the error of x {} equals its internal step {}",
        errors[0],
        step[0]
    );
    assert!(!min.user_state().has_covariance());
}

/// `a` is fixed and comes before the variable `b` and `c`, whose curvatures
/// differ by a factor of 100. A wrong index from the external to the internal
/// list would give `b` the step of `c`. The fixed parameter keeps its error.
#[test]
fn simplex_errors_follow_the_variable_parameters_after_a_fixed_one() {
    let calls: Calls = RefCell::new(Vec::new());
    let f = |p: &[f64]| {
        let value = (p[0] - 3.0).powi(2) + (p[1] - 1.0).powi(2) + 100.0 * (p[2] + 2.0).powi(2);
        calls.borrow_mut().push((p.to_vec(), value));
        value
    };
    let min = MnSimplex::new()
        .add("a", 0.7, 0.3)
        .add("b", 0.0, 0.5)
        .add("c", 0.0, 0.1)
        .fix(0)
        .tolerance(0.01)
        .minimize(&f);

    let vertices = final_simplex(&calls, 3, min.fval(), min.edm());
    let errors = min.user_state().errors();
    assert_eq!(errors.len(), 3);
    assert!(
        errors[0].total_cmp(&0.3).is_eq(),
        "the fixed parameter lost its error: {}",
        errors[0]
    );
    let expected_b = expected_error(&min, &vertices, 1);
    let expected_c = expected_error(&min, &vertices, 2);
    assert!(
        (expected_b / expected_c - 1.0).abs() > 0.1,
        "the steps of b and c are too close to tell them apart: {expected_b} and {expected_c}"
    );
    assert_relative(errors[1], expected_b, RECORDED_TOLERANCE, "error of b");
    assert_relative(errors[2], expected_c, RECORDED_TOLERANCE, "error of c");
    assert!(!min.user_state().has_covariance());
}

/// A run stopped by a small call limit. ROOT's behaviour at the call limit is
/// not measured here, because no ROOT file exists for this case, so the test
/// checks internal consistency only: the result is flagged, has no covariance,
/// and every error is finite and positive.
#[test]
fn simplex_at_the_call_limit_reports_flagged_finite_errors_without_covariance() {
    let min = MnSimplex::new()
        .add("x", 0.0, 0.5)
        .add("y", 0.0, 0.5)
        .max_fcn(12)
        .minimize(&|p: &[f64]| (p[0] - 1.0).powi(2) + 4.0 * (p[1] + 2.0).powi(2));

    assert!(min.reached_call_limit());
    assert!(!min.is_valid());
    assert!(min.nfcn() >= 12);
    assert!(!min.user_state().has_covariance());
    assert!(min.user_state().covariance().is_none());
    let errors = min.user_state().errors();
    assert_eq!(errors.len(), 2);
    for (i, error) in errors.iter().enumerate() {
        assert!(
            error.is_finite() && *error > 0.0,
            "error of parameter {i} is {error}"
        );
    }
}
