mod common;

use minuit2::{FCN, MnSimplex};

/// Rosenbrock function: f(x,y) = (1-x)^2 + 100(y-x^2)^2
/// Minimum at (1, 1) with f = 0.
#[test]
fn rosenbrock_2d() {
    let result = MnSimplex::new()
        .add("x", 0.0, 0.1)
        .add("y", 0.0, 0.1)
        // With the Simplex EDM criterion, use a tighter tolerance for
        // Rosenbrock to ensure full convergence.
        .tolerance(1.0e-5)
        .minimize(&|p: &[f64]| (1.0 - p[0]).powi(2) + 100.0 * (p[1] - p[0] * p[0]).powi(2));

    assert!(result.is_valid(), "minimization should converge");

    let params = result.params();
    assert!(
        (params[0] - 1.0).abs() < 0.1,
        "x should be near 1.0, got {}",
        params[0]
    );
    assert!(
        (params[1] - 1.0).abs() < 0.1,
        "y should be near 1.0, got {}",
        params[1]
    );
    assert!(
        result.fval() < 0.01,
        "fval should be near 0, got {}",
        result.fval()
    );
}

/// Simple quadratic bowl: f(x,y) = x^2 + y^2
/// Minimum at (0, 0) with f = 0.
#[test]
fn quadratic_bowl() {
    let result = MnSimplex::new()
        .add("x", 5.0, 1.0)
        .add("y", -3.0, 1.0)
        .minimize(&|p: &[f64]| p[0] * p[0] + p[1] * p[1]);

    assert!(result.is_valid(), "minimization should converge");

    let params = result.params();
    assert!(
        params[0].abs() < 0.01,
        "x should be near 0, got {}",
        params[0]
    );
    assert!(
        params[1].abs() < 0.01,
        "y should be near 0, got {}",
        params[1]
    );
    assert!(
        result.fval() < 1e-4,
        "fval should be near 0, got {}",
        result.fval()
    );
}

/// Gaussian fit to synthetic data.
/// Three parameters: amplitude, mean, sigma.
#[test]
fn gaussian_fit() {
    // Generate synthetic data: Gaussian with amp=10, mean=5, sigma=2
    let true_amp = 10.0;
    let true_mean = 5.0;
    let true_sigma = 2.0;

    let n_data = 50;
    let x_data: Vec<f64> = (0..n_data).map(|i| i as f64 * 0.2).collect();
    let y_data: Vec<f64> = x_data
        .iter()
        .map(|&x| true_amp * (-(x - true_mean).powi(2) / (2.0 * true_sigma * true_sigma)).exp())
        .collect();

    // Chi-square: sum of (model - data)^2
    let chi2 = move |p: &[f64]| -> f64 {
        let amp = p[0];
        let mean = p[1];
        let sigma = p[2];
        x_data
            .iter()
            .zip(y_data.iter())
            .map(|(&x, &y)| {
                let model = amp * (-(x - mean).powi(2) / (2.0 * sigma * sigma)).exp();
                (model - y).powi(2)
            })
            .sum()
    };

    let result = MnSimplex::new()
        .add("amp", 8.0, 1.0)
        .add("mean", 4.0, 0.5)
        .add_lower_limited("sigma", 1.5, 0.5, 0.01)
        .tolerance(0.001)
        .minimize(&chi2);

    assert!(result.is_valid(), "gaussian fit should converge");

    let state = result.user_state();
    let amp = state.value("amp").unwrap();
    let mean = state.value("mean").unwrap();
    let sigma = state.value("sigma").unwrap();

    assert!(
        (amp - true_amp).abs() < 0.5,
        "amp should be near {true_amp}, got {amp}"
    );
    assert!(
        (mean - true_mean).abs() < 0.5,
        "mean should be near {true_mean}, got {mean}"
    );
    assert!(
        (sigma - true_sigma).abs() < 0.5,
        "sigma should be near {true_sigma}, got {sigma}"
    );
}

/// Test bounded parameters with the sin transform.
#[test]
fn bounded_parameters() {
    // Minimize (x-3)^2 with x in [0, 5]
    let result = MnSimplex::new()
        .add_limited("x", 1.0, 0.5, 0.0, 5.0)
        .tolerance(0.001)
        .minimize(&|p: &[f64]| (p[0] - 3.0).powi(2));

    assert!(result.is_valid(), "bounded min should converge");
    let x = result.params()[0];
    assert!((x - 3.0).abs() < 0.1, "x should be near 3.0, got {x}");
    assert!((0.0..=5.0).contains(&x), "x should be within bounds");
}

/// Test with a fixed parameter.
#[test]
fn fixed_parameter() {
    struct QuadWithFixed;
    impl FCN for QuadWithFixed {
        fn value(&self, p: &[f64]) -> f64 {
            (p[0] - 2.0).powi(2) + (p[1] - 3.0).powi(2)
        }
    }

    let result = MnSimplex::new()
        .add("x", 0.0, 0.5)
        .add("y", 0.0, 0.5)
        .fix(1) // fix y at 0
        .minimize(&QuadWithFixed);

    assert!(result.is_valid());
    let params = result.params();
    // x should minimize to 2.0
    assert!(
        (params[0] - 2.0).abs() < 0.1,
        "x should be near 2.0, got {}",
        params[0]
    );
    // y should stay at 0.0 (fixed)
    assert!(
        (params[1] - 0.0).abs() < 1e-15,
        "y should be 0.0 (fixed), got {}",
        params[1]
    );
}

/// Display output should not panic.
#[test]
fn display_output() {
    let result = MnSimplex::new()
        .add("x", 1.0, 0.1)
        .minimize(&|p: &[f64]| p[0] * p[0]);

    common::assert_function_minimum_display(&result, &["x"]);
}

#[test]
fn respects_call_limit() {
    let result = MnSimplex::new()
        .add("x", 5.0, 1.0)
        .add("y", -3.0, 1.0)
        .max_fcn(4)
        .minimize(&|p: &[f64]| p[0] * p[0] + p[1] * p[1]);

    assert!(result.reached_call_limit());
    assert!(result.nfcn() >= 4);
}
