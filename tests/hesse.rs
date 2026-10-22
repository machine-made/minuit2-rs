use minuit2::{MnHesse, MnMigrad};

/// Quadratic: f(x,y) = a*x^2 + b*y^2
/// ROOT Minuit2 user covariance convention: V = 2 * up * H^-1.
/// Therefore sigma_x = sqrt(2*up / (2*a)) = sqrt(up/a), and similarly for y.
#[test]
fn hesse_quadratic_errors() {
    let a = 2.0;
    let b = 8.0;

    let result = MnMigrad::new()
        .add("x", 5.0, 1.0)
        .add("y", -3.0, 1.0)
        .minimize(&|p: &[f64]| a * p[0] * p[0] + b * p[1] * p[1]);

    assert!(result.is_valid());

    // Run Hesse
    let hesse_result =
        MnHesse::new().calculate(&|p: &[f64]| a * p[0] * p[0] + b * p[1] * p[1], &result);

    assert!(hesse_result.is_valid());

    let state = hesse_result.user_state();
    assert!(state.has_covariance(), "Hesse should produce covariance");

    // For chi-square: up = 1.0
    // H_xx = 2*a = 4, H_yy = 2*b = 16
    // V_xx = 2/H_xx = 1/2, V_yy = 2/H_yy = 1/8
    // sigma_x = sqrt(1/2) ≈ 0.7071, sigma_y = sqrt(1/8) ≈ 0.3536
    let err_x = state.error("x").unwrap();
    let err_y = state.error("y").unwrap();

    assert!(
        (err_x - std::f64::consts::FRAC_1_SQRT_2).abs() < 0.05,
        "sigma_x should be ~0.707, got {err_x}"
    );
    assert!(
        (err_y - 0.3535533905932738).abs() < 0.03,
        "sigma_y should be ~0.354, got {err_y}"
    );
}

/// Rosenbrock: Hesse after Migrad should give valid covariance.
#[test]
fn hesse_rosenbrock_valid() {
    let rosenbrock = |p: &[f64]| (1.0 - p[0]).powi(2) + 100.0 * (p[1] - p[0] * p[0]).powi(2);

    let result = MnMigrad::new()
        .add("x", 0.0, 0.1)
        .add("y", 0.0, 0.1)
        .minimize(&rosenbrock);

    assert!(result.is_valid());

    let hesse_result = MnHesse::new().calculate(&rosenbrock, &result);
    assert!(hesse_result.is_valid());

    let state = hesse_result.user_state();
    assert!(state.has_covariance());

    // At minimum (1,1): H = [[802, -400], [-400, 200]]
    // V = H^-1 = [[1/200, 1/200], [1/200, 401/20000]]
    // err_x = sqrt(1/200) ≈ 0.0707
    // Check errors are reasonable (not exact due to numerical Hessian)
    let err_x = state.error("x").unwrap();
    let err_y = state.error("y").unwrap();
    assert!(err_x > 0.01, "err_x should be positive, got {err_x}");
    assert!(err_y > 0.01, "err_y should be positive, got {err_y}");
}

#[test]
fn strategy2_migrad_returns_hesse_verified_covariance() {
    let rosenbrock = |p: &[f64]| (1.0 - p[0]).powi(2) + 100.0 * (p[1] - p[0] * p[0]).powi(2);

    let result = MnMigrad::new()
        .add("x", -1.2, 0.2)
        .add("y", 1.0, 0.2)
        .with_strategy(2)
        .max_fcn(100_000)
        .minimize(&rosenbrock);

    assert!(
        result.is_valid(),
        "strategy-2 Migrad should converge: fval={} edm={} nfcn={}",
        result.fval(),
        result.edm(),
        result.nfcn()
    );

    let explicit = MnHesse::new()
        .with_strategy(2)
        .with_max_calls(100_000)
        .calculate(&rosenbrock, &result);
    assert!(explicit.is_valid());

    let returned_cov = result.state().error().matrix();
    let explicit_cov = explicit.state().error().matrix();
    assert_eq!(returned_cov.shape(), explicit_cov.shape());
    for i in 0..returned_cov.nrows() {
        for j in 0..returned_cov.ncols() {
            let got = returned_cov[(i, j)];
            let want = explicit_cov[(i, j)];
            let scale = want.abs().max(1.0);
            assert!(
                (got - want).abs() <= 2e-3 * scale,
                "cov[{i},{j}] returned {got:.12e} vs explicit Hesse {want:.12e}"
            );
        }
    }
}

/// Global correlations on correlated quadratic.
#[test]
fn hesse_global_correlations() {
    // f(x,y) = x^2 + y^2 + x*y  (correlated)
    let result = MnMigrad::new()
        .add("x", 5.0, 1.0)
        .add("y", -3.0, 1.0)
        .minimize(&|p: &[f64]| p[0] * p[0] + p[1] * p[1] + p[0] * p[1]);

    assert!(result.is_valid());

    let hesse_result = MnHesse::new().calculate(
        &|p: &[f64]| p[0] * p[0] + p[1] * p[1] + p[0] * p[1],
        &result,
    );

    let state = hesse_result.user_state();
    if let Some(gcc) = state.global_cc() {
        // With cross-term, gcc should be non-zero
        assert!(
            gcc[0] > 0.1,
            "gcc[0] should show correlation, got {}",
            gcc[0]
        );
        assert!(
            gcc[1] > 0.1,
            "gcc[1] should show correlation, got {}",
            gcc[1]
        );
    }
}

/// Hesse with calculate_errors (doesn't modify minimum).
#[test]
fn hesse_calculate_errors() {
    let result = MnMigrad::new()
        .add("x", 5.0, 1.0)
        .add("y", -3.0, 1.0)
        .minimize(&|p: &[f64]| 2.0 * p[0] * p[0] + 8.0 * p[1] * p[1]);

    let state = MnHesse::new()
        .calculate_errors(&|p: &[f64]| 2.0 * p[0] * p[0] + 8.0 * p[1] * p[1], &result);

    assert!(state.has_covariance());
    let err_x = state.error("x").unwrap();
    assert!(
        (err_x - std::f64::consts::FRAC_1_SQRT_2).abs() < 0.05,
        "sigma_x should be ~0.707, got {err_x}"
    );
}

#[test]
fn hesse_flat_parameter_returns_failed_state_without_covariance() {
    let minimum = MnMigrad::new()
        .add("x", 2.0, 0.2)
        .add("flat", 1.0, 0.2)
        .minimize(&|p: &[f64]| p[0] * p[0]);

    let hesse = MnHesse::new().calculate(&|p: &[f64]| p[0] * p[0], &minimum);

    assert!(
        !hesse.is_valid(),
        "flat direction must not produce a valid Hesse covariance"
    );
    assert!(
        !hesse.user_state().has_covariance(),
        "failed Hesse state must not expose covariance"
    );
}
