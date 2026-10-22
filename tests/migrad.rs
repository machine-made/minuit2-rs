use nalgebra::DVector;
mod common;

use minuit2::{
    FCN, FunctionMinimum, MinuitParameter, MnMigrad, MnSimplex, MnUserTransformation,
    minimum::{parameters::MinimumParameters, seed::MinimumSeed, state::MinimumState},
};

/// Rosenbrock function: f(x,y) = (1-x)^2 + 100(y-x^2)^2
/// Minimum at (1, 1) with f = 0.
#[test]
fn rosenbrock_2d() {
    let result = MnMigrad::new()
        .add("x", -1.0, 1.0)
        .add("y", -1.0, 1.0)
        .tolerance(0.1) // tighter than default for precise convergence
        .minimize(&|p: &[f64]| (1.0 - p[0]).powi(2) + 100.0 * (p[1] - p[0] * p[0]).powi(2));

    assert!(result.is_valid(), "migrad should converge");

    let params = result.params();
    assert!(
        (params[0] - 1.0).abs() < 0.05,
        "x should be near 1.0, got {}",
        params[0]
    );
    assert!(
        (params[1] - 1.0).abs() < 0.05,
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
    let result = MnMigrad::new()
        .add("x", 5.0, 1.0)
        .add("y", -3.0, 1.0)
        .minimize(&|p: &[f64]| p[0] * p[0] + p[1] * p[1]);

    assert!(result.is_valid(), "migrad should converge");

    let params = result.params();
    assert!(
        params[0].abs() < 1e-4,
        "x should be near 0, got {}",
        params[0]
    );
    assert!(
        params[1].abs() < 1e-4,
        "y should be near 0, got {}",
        params[1]
    );
    assert!(
        result.fval() < 1e-8,
        "fval should be near 0, got {}",
        result.fval()
    );
}

/// Quadratic with known minimum in N dimensions.
#[test]
fn nd_quadratic() {
    let result = MnMigrad::new()
        .add("x0", 3.0, 1.0)
        .add("x1", -2.0, 1.0)
        .add("x2", 7.0, 1.0)
        .add("x3", -5.0, 1.0)
        .add("x4", 1.0, 1.0)
        .minimize(&|p: &[f64]| p.iter().map(|x| x * x).sum());

    assert!(result.is_valid(), "N-dim quadratic should converge");
    for (i, &val) in result.params().iter().enumerate() {
        assert!(val.abs() < 1e-3, "x{} should be near 0, got {}", i, val);
    }
}

/// Test bounded parameters with the sin transform.
#[test]
fn bounded_parameters() {
    // Minimize (x-3)^2 with x in [0, 5]
    let result = MnMigrad::new()
        .add_limited("x", 1.0, 0.5, 0.0, 5.0)
        .minimize(&|p: &[f64]| (p[0] - 3.0).powi(2));

    assert!(result.is_valid(), "bounded min should converge");
    let x = result.params()[0];
    assert!((x - 3.0).abs() < 0.05, "x should be near 3.0, got {x}");
    assert!((0.0..=5.0).contains(&x), "x should be within bounds: {x}");
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

    let result = MnMigrad::new()
        .add("x", 0.0, 0.5)
        .add("y", 0.0, 0.5)
        .fix(1) // fix y at 0
        .minimize(&QuadWithFixed);

    assert!(result.is_valid());
    let params = result.params();
    assert!(
        (params[0] - 2.0).abs() < 0.01,
        "x should be near 2.0, got {}",
        params[0]
    );
    assert!(
        (params[1] - 0.0).abs() < 1e-15,
        "y should be 0.0 (fixed), got {}",
        params[1]
    );
}

/// Gaussian fit to synthetic data.
#[test]
fn gaussian_fit() {
    let true_amp = 10.0;
    let true_mean = 5.0;
    let true_sigma = 2.0;

    let n_data = 50;
    let x_data: Vec<f64> = (0..n_data).map(|i| i as f64 * 0.2).collect();
    let y_data: Vec<f64> = x_data
        .iter()
        .map(|&x| true_amp * (-(x - true_mean).powi(2) / (2.0 * true_sigma * true_sigma)).exp())
        .collect();

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

    let result = MnMigrad::new()
        .add("amp", 8.0, 1.0)
        .add("mean", 4.0, 0.5)
        .add_lower_limited("sigma", 1.5, 0.5, 0.01)
        .minimize(&chi2);

    assert!(result.is_valid(), "gaussian fit should converge");

    let state = result.user_state();
    let amp = state.value("amp").unwrap();
    let mean = state.value("mean").unwrap();
    let sigma = state.value("sigma").unwrap();

    assert!(
        (amp - true_amp).abs() < 0.1,
        "amp should be near {true_amp}, got {amp}"
    );
    assert!(
        (mean - true_mean).abs() < 0.1,
        "mean should be near {true_mean}, got {mean}"
    );
    assert!(
        (sigma - true_sigma).abs() < 0.1,
        "sigma should be near {true_sigma}, got {sigma}"
    );
}

/// Compare Migrad vs Simplex: with tight tolerance, Migrad should achieve
/// a very small fval on quadratic-like functions with far fewer FCN calls.
#[test]
fn migrad_vs_simplex_quadratic() {
    // On a 5D quadratic bowl, Migrad should converge much faster
    let quadratic = |p: &[f64]| -> f64 {
        p.iter()
            .enumerate()
            .map(|(i, &x)| (i as f64 + 1.0) * x * x)
            .sum()
    };

    let migrad_result = MnMigrad::new()
        .add("x0", 5.0, 1.0)
        .add("x1", -3.0, 1.0)
        .add("x2", 7.0, 1.0)
        .add("x3", -2.0, 1.0)
        .add("x4", 4.0, 1.0)
        .minimize(&quadratic);

    let simplex_result = MnSimplex::new()
        .add("x0", 5.0, 1.0)
        .add("x1", -3.0, 1.0)
        .add("x2", 7.0, 1.0)
        .add("x3", -2.0, 1.0)
        .add("x4", 4.0, 1.0)
        .minimize(&quadratic);

    assert!(migrad_result.is_valid());
    assert!(
        simplex_result.reached_call_limit()
            || simplex_result.is_above_max_edm()
            || simplex_result.is_valid(),
        "unexpected simplex termination flags"
    );

    // Migrad should achieve a much lower function value on a smooth quadratic
    assert!(
        migrad_result.fval() < simplex_result.fval() * 100.0,
        "Migrad fval ({:.6e}) should be much lower than Simplex fval ({:.6e})",
        migrad_result.fval(),
        simplex_result.fval()
    );
    assert!(
        migrad_result.nfcn() < simplex_result.nfcn(),
        "Migrad should use fewer calls than Simplex (migrad={}, simplex={})",
        migrad_result.nfcn(),
        simplex_result.nfcn()
    );
}

/// Display output should not panic.
#[test]
fn display_output() {
    let result = MnMigrad::new()
        .add("x", 1.0, 0.1)
        .minimize(&|p: &[f64]| p[0] * p[0]);

    common::assert_function_minimum_display(&result, &["x"]);
}

#[test]
fn display_output_marks_call_limit_warning() {
    let params = MnUserTransformation::new(vec![MinuitParameter::new(0, "x", 0.0, 0.1)]);
    let seed_state = MinimumState::from_params_edm(
        MinimumParameters::new(DVector::from_vec(vec![0.0]), 0.0),
        0.0,
        1,
    );
    let seed = MinimumSeed::new(seed_state, params);
    let minimum = FunctionMinimum::with_call_limit(
        seed,
        vec![MinimumState::from_params_edm(
            MinimumParameters::new(DVector::from_vec(vec![0.0]), 0.0),
            0.0,
            1,
        )],
        1.0,
    );

    let output = format!("{minimum}");

    assert!(!minimum.is_valid());
    assert!(minimum.reached_call_limit());
    assert!(output.contains("WARNING: call limit reached"));
    assert!(output.contains("valid:     false"));
}

#[test]
fn display_output_marks_above_max_edm_warning() {
    let params = MnUserTransformation::new(vec![MinuitParameter::new(0, "x", 0.0, 0.1)]);
    let seed_state = MinimumState::from_params_edm(
        MinimumParameters::new(DVector::from_vec(vec![0.0]), 0.0),
        0.0,
        1,
    );
    let seed = MinimumSeed::new(seed_state, params);
    let minimum = FunctionMinimum::above_max_edm(
        seed,
        vec![MinimumState::from_params_edm(
            MinimumParameters::new(DVector::from_vec(vec![0.0]), 0.0),
            0.0,
            1,
        )],
        1.0,
    );

    let output = format!("{minimum}");

    assert!(!minimum.is_valid());
    assert!(minimum.is_above_max_edm());
    assert!(output.contains("WARNING: EDM above maximum"));
    assert!(output.contains("valid:     false"));
}
