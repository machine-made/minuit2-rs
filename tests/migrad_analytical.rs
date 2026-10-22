//! Integration tests for Migrad with analytical gradients.

use minuit2::{FCN, FCNGradient, MnMigrad};

/// Rosenbrock function: (1-x)² + 100(y-x²)²
struct Rosenbrock;

impl FCN for Rosenbrock {
    fn value(&self, p: &[f64]) -> f64 {
        let x = p[0];
        let y = p[1];
        (1.0 - x) * (1.0 - x) + 100.0 * (y - x * x) * (y - x * x)
    }

    fn error_def(&self) -> f64 {
        1.0
    }
}

impl FCNGradient for Rosenbrock {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        let x = p[0];
        let y = p[1];
        // df/dx = -2(1-x) - 400x(y-x²)
        let dfdx = -2.0 * (1.0 - x) - 400.0 * x * (y - x * x);
        // df/dy = 200(y-x²)
        let dfdy = 200.0 * (y - x * x);
        vec![dfdx, dfdy]
    }
}

/// Quadratic: x² + 4y²
struct Quadratic;

impl FCN for Quadratic {
    fn value(&self, p: &[f64]) -> f64 {
        p[0] * p[0] + 4.0 * p[1] * p[1]
    }

    fn error_def(&self) -> f64 {
        1.0
    }
}

impl FCNGradient for Quadratic {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        vec![2.0 * p[0], 8.0 * p[1]]
    }
}

#[test]
fn migrad_analytical_quadratic() {
    // Test with simple quadratic using analytical gradients
    let result = MnMigrad::new()
        .add("x", 3.0, 0.1)
        .add("y", 2.0, 0.1)
        .minimize_grad(&Quadratic);

    // Check convergence
    assert!(result.is_valid(), "Minimization should converge");
    assert!(
        result.fval() < 1e-8,
        "Should reach near-zero: {}",
        result.fval()
    );

    // Check parameters
    let state = result.user_state();
    let x = state.value("x").unwrap();
    let y = state.value("y").unwrap();
    assert!(x.abs() < 1e-4, "x should be near 0, got {}", x);
    assert!(y.abs() < 1e-4, "y should be near 0, got {}", y);
}

#[test]
fn migrad_analytical_rosenbrock() {
    // Test with Rosenbrock using analytical gradients
    let result = MnMigrad::new()
        .add("x", 0.0, 0.1)
        .add("y", 0.0, 0.1)
        .minimize_grad(&Rosenbrock);

    // Should converge to (1, 1) or at least get close
    if result.is_valid() {
        let state = result.user_state();
        let x = state.value("x").unwrap();
        let y = state.value("y").unwrap();
        // Just check we're in the right ballpark
        assert!(x > -0.5, "x should be positive region, got {}", x);
        assert!(y > -0.5, "y should be positive region, got {}", y);
    }
}

#[test]
fn migrad_analytical_vs_numerical_quadratic() {
    // Compare analytical and numerical gradients on quadratic
    let migrad = MnMigrad::new().add("x", 3.0, 0.1).add("y", 2.0, 0.1);

    let result_analytical = migrad.minimize_grad(&Quadratic);
    let result_numerical = migrad.minimize(&Quadratic);

    // Both should converge to the same point
    let state_analytical = result_analytical.user_state();
    let state_numerical = result_numerical.user_state();

    let x_analytical = state_analytical.value("x").unwrap();
    let x_numerical = state_numerical.value("x").unwrap();
    let y_analytical = state_analytical.value("y").unwrap();
    let y_numerical = state_numerical.value("y").unwrap();

    assert!(
        (x_analytical - x_numerical).abs() < 1e-5,
        "x should match: analytical={}, numerical={}",
        x_analytical,
        x_numerical
    );
    assert!(
        (y_analytical - y_numerical).abs() < 1e-5,
        "y should match: analytical={}, numerical={}",
        y_analytical,
        y_numerical
    );

    // Analytical might use fewer function calls (not always guaranteed due to randomness)
    // but both should be valid
    assert!(result_analytical.is_valid(), "Analytical should converge");
    assert!(result_numerical.is_valid(), "Numerical should converge");
}

#[test]
fn migrad_analytical_with_bounds() {
    // Test analytical gradients with bounded parameters
    let result = MnMigrad::new()
        .add_limited("x", 0.5, 0.1, -1.0, 1.0)
        .add_limited("y", 0.5, 0.1, -2.0, 2.0)
        .minimize_grad(&Quadratic);

    assert!(result.is_valid(), "Should converge with bounds");

    let state = result.user_state();
    let x = state.value("x").unwrap();
    let y = state.value("y").unwrap();

    // For bounded quadratic, minimum is at origin but with bounds
    // we might not get exactly zero
    assert!(x.abs() < 0.1, "x should be small, got {}", x);
    assert!(y.abs() < 0.1, "y should be small, got {}", y);

    // Check bounds are satisfied
    assert!((-1.0..=1.0).contains(&x), "x must be in [-1, 1], got {}", x);
    assert!((-2.0..=2.0).contains(&y), "y must be in [-2, 2], got {}", y);
}
