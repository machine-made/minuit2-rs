use minuit2::{MnHesse, MnMigrad};

/// Quadratic: f(x,y) = a*x^2 + b*y^2
/// ROOT Minuit2 user covariance convention: V = 2 * up * H^-1.
/// Therefore `sigma_x` = sqrt(2*up / (2*a)) = sqrt(up/a), and similarly for y.
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
        (err_y - 0.353_553_390_593_273_8).abs() < 0.03,
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

// ROOT Minuit2 v6-36-08 on the same workload: strategy-2 Migrad on Rosenbrock, then a
// second explicit strategy-2 Hesse pass. ROOT keeps the refined gradient steps of the
// first pass and the second pass starts from them, so the error matrix moves.
// The constants below are the `migrad_inverse_hessian` and `hesse_inverse_hessian`
// entries of reports/verification/raw/trace/repeated_hesse_strategy2.root.json.
const ROOT_AFTER_MIGRAD: [[f64; 2]; 2] = [
    [0.509_400_261_528_040_9, 1.017_172_600_092_039],
    [1.017_172_600_092_039, 2.036_094_556_713_421],
];
const ROOT_AFTER_HESSE: [[f64; 2]; 2] = [
    [0.505_347_499_645_794_7, 1.009_070_067_042_981_7],
    [1.009_070_067_042_981_7, 2.019_895_494_517_758],
];

/// Every entry of `matrix` must be within a relative 1e-5 of the ROOT entry.
/// The scale is the absolute value of the ROOT entry or 1.0, whichever is larger.
/// Rust and ROOT run the same algorithm in floating point and differ in the seventh
/// digit. For entry [0,0], whose scale is 1.0, the measured differences are 6.6e-8 and
/// 5.2e-8; the largest over all entries is 1.2e-7. The tolerance leaves a margin of
/// about 80. It still fails if a refinement step is wrong. A first probe step of 0.5
/// times the stored step instead of 0.2 moves entry [0,0] by 3.8e-4. Ignoring the
/// stored steps in the second pass moves it by 0.004.
fn assert_matches_root(matrix: &nalgebra::DMatrix<f64>, root: &[[f64; 2]; 2], label: &str) {
    assert_eq!(matrix.shape(), (2, 2));
    for i in 0..2 {
        for j in 0..2 {
            let got = matrix[(i, j)];
            let want = root[i][j];
            let scale = want.abs().max(1.0);
            assert!(
                (got - want).abs() <= 1e-5 * scale,
                "{label}[{i},{j}] is {got:.16e}, ROOT has {want:.16e}"
            );
        }
    }
}

#[test]
fn strategy2_repeated_hesse_matches_root() {
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

    let returned = result.state().error().matrix();
    let repeated = explicit.state().error().matrix();
    assert_matches_root(returned, &ROOT_AFTER_MIGRAD, "Migrad error matrix");
    assert_matches_root(repeated, &ROOT_AFTER_HESSE, "second Hesse error matrix");

    // ROOT's second pass lowers entry [0,0] by 0.0040528.
    let shift = returned[(0, 0)] - repeated[(0, 0)];
    assert!(
        (0.004..0.0041).contains(&shift),
        "second Hesse should lower entry [0,0] by about 0.0040528, got {shift:.7e}"
    );
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

/// Hesse with `calculate_errors` (doesn't modify minimum).
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
