use minuit2::{MnMigrad, MnSimplex};

/// A function that returns NaN sometimes.
/// The optimizer should handle this gracefully (e.g., by retreating or retrying).
#[test]
fn nan_resilience() {
    // f(x) = x^2, but returns NaN if x > 5.0
    let result = MnMigrad::new()
        .add("x", 4.0, 1.0) // Start close to the danger zone
        .minimize(&|p: &[f64]| {
            if p[0] > 5.0 { f64::NAN } else { p[0] * p[0] }
        });

    // It should not panic. Ideally it converges to 0.
    // If it steps into NaN territory, it should treat it as a high value.
    assert!(result.is_valid());
    assert!(result.fval() < 1e-4);
}

/// A function that returns Infinity sometimes.
#[test]
fn inf_resilience() {
    let result = MnMigrad::new().add("x", 4.0, 1.0).minimize(&|p: &[f64]| {
        if p[0] > 5.0 {
            f64::INFINITY
        } else {
            p[0] * p[0]
        }
    });

    assert!(result.is_valid());
    assert!(result.fval() < 1e-4);
}

/// Stress test with high dimensionality (50 parameters).
/// Ensures that the nalgebra matrix operations and generic logic hold up.
#[test]
fn high_dim_stress() {
    let n = 50;
    let mut builder = MnMigrad::new();

    // Add 50 parameters: x0, x1, ...
    for i in 0..n {
        builder = builder.add(format!("x{}", i), i as f64, 0.1);
    }

    // Minimize sum(x_i^2)
    let result = builder
        .max_fcn(10000) // Increase budget for high dim
        .minimize(&|p: &[f64]| p.iter().map(|x| x * x).sum());

    assert!(result.is_valid());
    assert!(result.fval() < 1e-4);

    for i in 0..n {
        assert!(result.user_state().value(&format!("x{}", i)).unwrap().abs() < 1e-2);
    }
}

/// Goldstein-Price function.
/// A classic test function for global optimization.
/// Usually minimized on [-2, 2]. Global minimum at (0, -1) is 3.
#[test]
fn goldstein_price() {
    let gp = |p: &[f64]| {
        let x = p[0];
        let y = p[1];

        let part1 = 1.0
            + (x + y + 1.0).powi(2)
                * (19.0 - 14.0 * x + 3.0 * x * x - 14.0 * y + 6.0 * x * y + 3.0 * y * y);
        let part2 = 30.0
            + (2.0 * x - 3.0 * y).powi(2)
                * (18.0 - 32.0 * x + 12.0 * x * x + 48.0 * y - 36.0 * x * y + 27.0 * y * y);

        part1 * part2
    };

    // Start somewhat near the minimum (0, -1).
    // Use MnSimplex because Goldstein-Price has very steep walls that confuse Migrad's gradient logic.
    let result = MnSimplex::new()
        .add("x", 0.5, 0.5)
        .add("y", -0.5, 0.5)
        .tolerance(0.0001) // Simplex needs tighter tolerance to pin down the exact minimum
        .max_fcn(5000)
        .minimize(&gp);

    assert!(result.is_valid());

    // Global minimum is 3.0
    assert!((result.fval() - 3.0).abs() < 1e-4);

    let params = result.params();
    assert!((params[0] - 0.0).abs() < 0.1);
    assert!((params[1] - (-1.0)).abs() < 0.1);
}

/// Test boundary edge cases.
/// Starting EXACTLY on a boundary should be fine.
/// Parameters pushed against a boundary should stick.
#[test]
fn boundary_edge_case() {
    // Minimize (x-5)^2 with x in [0, 5]. Minimum is at 5 (the boundary).

    // Case 1: Start at the boundary (5.0)
    let result_at_bound = MnMigrad::new()
        .add_limited("x", 5.0, 0.1, 0.0, 5.0)
        .minimize(&|p: &[f64]| (p[0] - 5.0).powi(2));

    assert!(result_at_bound.is_valid());
    assert!((result_at_bound.params()[0] - 5.0).abs() < 1e-4);

    // Case 2: Start outside the boundary? (Should panic or clamp?)
    // The API documentation says we must provide valid start values.
    // Let's test starting VERY close to the boundary.
    let result_near_bound = MnMigrad::new()
        .add_limited("x", 4.999999, 0.1, 0.0, 5.0)
        .minimize(&|p: &[f64]| (p[0] - 5.0).powi(2));

    assert!(result_near_bound.is_valid());
    assert!((result_near_bound.params()[0] - 5.0).abs() < 1e-4);
}

/// Rosenbrock with the classic difficult start point (-1.2, 1.0).
#[test]
fn rosenbrock_hard_start() {
    let result = MnMigrad::new()
        .add("x", -1.2, 0.1)
        .add("y", 1.0, 0.1)
        .tolerance(0.1)
        .minimize(&|p: &[f64]| (1.0 - p[0]).powi(2) + 100.0 * (p[1] - p[0] * p[0]).powi(2));

    assert!(result.is_valid());
    assert!(result.fval() < 1e-4);
    assert!((result.params()[0] - 1.0).abs() < 1e-2);
    assert!((result.params()[1] - 1.0).abs() < 1e-2);
}
