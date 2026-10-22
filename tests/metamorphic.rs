use minuit2::{MnHesse, MnMigrad};

fn assert_close(got: f64, want: f64, tol: f64, label: &str) {
    assert!(
        (got - want).abs() <= tol,
        "{label}: expected {want}, got {got}, diff={}",
        (got - want).abs()
    );
}

fn assert_params_close(got: &[f64], want: &[f64], tol: f64, label: &str) {
    assert_eq!(got.len(), want.len(), "{label}: parameter length mismatch");
    for (i, (&g, &w)) in got.iter().zip(want.iter()).enumerate() {
        assert_close(g, w, tol, &format!("{label}[{i}]"));
    }
}

/// Adding a constant to an objective must not move the argmin; it only shifts fval.
#[test]
fn migrad_argmin_invariant_under_objective_translation() {
    let base = |p: &[f64]| (p[0] - 2.5).powi(2) + 3.0 * (p[1] + 1.25).powi(2);
    let shifted = |p: &[f64]| base(p) + 123.75;

    let base_min = MnMigrad::new()
        .add("x", -4.0, 0.8)
        .add("y", 3.0, 0.8)
        .minimize(&base);
    let shifted_min = MnMigrad::new()
        .add("x", -4.0, 0.8)
        .add("y", 3.0, 0.8)
        .minimize(&shifted);

    assert!(base_min.is_valid());
    assert!(shifted_min.is_valid());
    let base_params = base_min.params();
    let shifted_params = shifted_min.params();
    assert_params_close(&shifted_params, &base_params, 1e-7, "translated argmin");
    assert_close(
        shifted_min.fval() - base_min.fval(),
        123.75,
        1e-8,
        "translated fval shift",
    );
}

/// Multiplying a positive objective by a positive scalar must preserve the argmin
/// and scale the reported fval by the same scalar.
#[test]
fn migrad_argmin_invariant_under_positive_objective_scaling() {
    let base = |p: &[f64]| 2.0 * (p[0] - 1.5).powi(2) + 0.5 * (p[1] + 4.0).powi(2);
    let scaled = |p: &[f64]| 7.0 * base(p);

    let base_min = MnMigrad::new()
        .add("x", 8.0, 1.0)
        .add("y", -8.0, 1.0)
        .minimize(&base);
    let scaled_min = MnMigrad::new()
        .add("x", 8.0, 1.0)
        .add("y", -8.0, 1.0)
        .minimize(&scaled);

    assert!(base_min.is_valid());
    assert!(scaled_min.is_valid());
    let base_params = base_min.params();
    let scaled_params = scaled_min.params();
    assert_params_close(&scaled_params, &base_params, 1e-7, "scaled argmin");
    assert_close(
        scaled_min.fval(),
        7.0 * base_min.fval(),
        1e-8,
        "scaled fval",
    );
}

/// For a separable quadratic, permuting parameter order should permute the solution
/// in the same way. This catches index/name mixups in builders and result mapping.
#[test]
fn migrad_separable_quadratic_commutes_with_parameter_permutation() {
    let original = |p: &[f64]| {
        1.0 * (p[0] - 1.0).powi(2) + 2.0 * (p[1] + 2.0).powi(2) + 4.0 * (p[2] - 3.0).powi(2)
    };
    let permuted = |p: &[f64]| {
        4.0 * (p[0] - 3.0).powi(2) + 1.0 * (p[1] - 1.0).powi(2) + 2.0 * (p[2] + 2.0).powi(2)
    };

    let original_min = MnMigrad::new()
        .add("a", -5.0, 1.0)
        .add("b", 6.0, 1.0)
        .add("c", -7.0, 1.0)
        .minimize(&original);
    let permuted_min = MnMigrad::new()
        .add("c", -7.0, 1.0)
        .add("a", -5.0, 1.0)
        .add("b", 6.0, 1.0)
        .minimize(&permuted);

    assert!(original_min.is_valid());
    assert!(permuted_min.is_valid());

    let p = permuted_min.params();
    let permuted_back = [p[1], p[2], p[0]];
    let original_params = original_min.params();
    assert_params_close(&permuted_back, &original_params, 1e-7, "permuted solution");
}

/// Hesse refines covariance/error estimates; it must not move the fitted parameters.
#[test]
fn hesse_preserves_migrad_parameter_values() {
    let objective =
        |p: &[f64]| (p[0] - 2.0).powi(2) + 5.0 * (p[1] + 0.5).powi(2) + 0.25 * p[0] * p[1];

    let migrad = MnMigrad::new()
        .add("x", -3.0, 0.5)
        .add("y", 4.0, 0.5)
        .minimize(&objective);
    let hesse = MnHesse::new().calculate(&objective, &migrad);

    assert!(migrad.is_valid());
    assert!(hesse.is_valid());
    let migrad_params = migrad.params();
    let hesse_params = hesse.params();
    assert_params_close(
        &hesse_params,
        &migrad_params,
        1e-12,
        "hesse parameter preservation",
    );
    assert!(hesse.user_state().has_covariance());
}
