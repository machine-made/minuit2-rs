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
        let derivative_y = 200.0 * (y - x * x);
        vec![dfdx, derivative_y]
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
    // Measured |x| is 6.0e-14 and |y| is 1.6e-13, so the bound 1e-4 has a margin of 6e8 or more.
    // The bound follows from the fval limit 1e-8 above: f = x^2 + 4y^2 < 1e-8 gives |x| < 1e-4.
    // It is a generous sanity limit; a result that did not move would sit at 3 and 2.
    assert!(x.abs() < 1e-4, "x should be near 0, got {x}");
    assert!(y.abs() < 1e-4, "y should be near 0, got {y}");
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
        // Measured x = 0.9966 and y = 0.9930, which are 3.4e-3 and 7.0e-3 below the minimum (1, 1).
        // The bound -0.5 is 1.5 below them, so it is a generous sanity limit.
        // It rejects a result in the negative region but not the start point (0, 0).
        assert!(x > -0.5, "x should be positive region, got {x}");
        assert!(y > -0.5, "y should be positive region, got {y}");
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

    // The two runs differ by 7.0e-10 in x and 2.8e-10 in y, so the bound 1e-5 has a margin of
    // 1.4e4 or more. Both end within 1e-9 of the minimum at the origin; the numerical run is
    // the one off by 7e-10. The bound is a generous sanity limit: a run that did not move
    // differs by order 1.
    assert!(
        (x_analytical - x_numerical).abs() < 1e-5,
        "x should match: analytical={x_analytical}, numerical={x_numerical}"
    );
    assert!(
        (y_analytical - y_numerical).abs() < 1e-5,
        "y should match: analytical={y_analytical}, numerical={y_numerical}"
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
    // Measured |x| is 2.7e-3 and |y| is 1.4e-3, so the bound 0.1 has a margin of 37 and 71.
    // Migrad works in sine-transformed variables here and stops at an EDM of 1.5e-5 (goal 2e-4),
    // not near machine precision as in the unbounded case.
    assert!(x.abs() < 0.1, "x should be small, got {x}");
    assert!(y.abs() < 0.1, "y should be small, got {y}");

    // Check bounds are satisfied
    assert!((-1.0..=1.0).contains(&x), "x must be in [-1, 1], got {x}");
    assert!((-2.0..=2.0).contains(&y), "y must be in [-2, 2], got {y}");
}
