//! Simplex errors and covariance against ROOT Minuit2 `SimplexBuilder.cxx` and
//! `MnUserParameterState.cxx`: the parameter error is the final simplex step
//! scaled by `sqrt(up / edm)`, and the result has no covariance.

use std::cell::RefCell;
use std::path::Path;
use std::process::Command;

use minuit2::MnSimplex;
use minuit2::minimum::error::ErrorMatrixStatus;
use serde_json::Value;

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
            1e-6,
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
    assert_relative(errors[0], expected, 1e-12, "error of x");

    let step = min.state().parameters().step()[0];
    assert_relative(errors[0], step, 1e-12, "error of x against the final step");

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
