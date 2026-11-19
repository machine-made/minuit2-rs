//! Covariance of Rosenbrock after Migrad with a user gradient, against the closed form.

use minuit2::{FCN, FCNGradient, MnHesse, MnMigrad};
use nalgebra::Matrix2;

struct Rosenbrock;

impl FCN for Rosenbrock {
    fn value(&self, p: &[f64]) -> f64 {
        100.0 * (p[1] - p[0] * p[0]).powi(2) + (1.0 - p[0]).powi(2)
    }

    fn error_def(&self) -> f64 {
        1.0
    }
}

impl FCNGradient for Rosenbrock {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        vec![
            -400.0 * p[0] * (p[1] - p[0] * p[0]) + 2.0 * (p[0] - 1.0),
            200.0 * (p[1] - p[0] * p[0]),
        ]
    }
}

/// Covariance at the minimum (1, 1) for error definition 1, from the closed form.
///
/// The Hessian is [[802, -400], [-400, 200]] with determinant 802 * 200 - 400 * 400 = 400.
/// Its inverse is [[200, 400], [400, 802]] / 400 = [[1/2, 1], [1, 401/200]].
/// The covariance is twice the error definition times the inverse: [[1, 2], [2, 401/100]].
const REFERENCE: [[f64; 2]; 2] = [[1.0, 2.0], [2.0, 4.01]];

/// Largest entry error divided by the geometric mean of the matching reference diagonals.
///
/// This is the `cov_error` metric of the comparison report.
fn covariance_error(covariance: &[[f64; 2]; 2]) -> f64 {
    let mut worst = 0.0_f64;
    for (i, row) in covariance.iter().enumerate() {
        for (j, entry) in row.iter().enumerate() {
            let scale = (REFERENCE[i][i] * REFERENCE[j][j]).sqrt();
            worst = worst.max((entry - REFERENCE[i][j]).abs() / scale);
        }
    }
    worst
}

fn entries(minimum: &minuit2::FunctionMinimum) -> [[f64; 2]; 2] {
    let covariance = minimum
        .user_state()
        .covariance()
        .expect("the minimum has a covariance");
    assert_eq!(covariance.nrow(), 2);
    [
        [covariance.get(0, 0), covariance.get(0, 1)],
        [covariance.get(1, 0), covariance.get(1, 1)],
    ]
}

#[test]
fn gradient_covariance_matches_reference() {
    // The literals above must be the inverse of the Hessian, times twice the error definition.
    let hessian = Matrix2::new(802.0, -400.0, -400.0, 200.0);
    assert_eq!(hessian.determinant(), 400.0);
    let reference = Matrix2::new(
        REFERENCE[0][0],
        REFERENCE[0][1],
        REFERENCE[1][0],
        REFERENCE[1][1],
    );
    assert!((reference * hessian - Matrix2::identity() * 2.0).amax() < 1e-12);

    // Same case as the gradient row of the comparison report: start (0, 0), steps 0.1,
    // strategy 1, tolerance 1e-4, no limits.
    let minimum = MnMigrad::new()
        .with_strategy(1)
        .tolerance(1e-4)
        .add("p0", 0.0, 0.1)
        .add("p1", 0.0, 0.1)
        .minimize_grad(&Rosenbrock);
    assert!(minimum.is_valid(), "Migrad did not converge");

    // The report scores the covariance that Migrad returns, so this test does the same.
    // iminuit 2.33.0 gives an error of 0.0382 on this case, measured by the comparison harness.
    // This crate gave 0.662 before the Hesse gradient steps and refinement followed ROOT.
    // It gives 0.0399 now. The bound 0.0765 is twice the iminuit value.
    let error = covariance_error(&entries(&minimum));
    assert!(
        error <= 0.0765,
        "covariance error {error} is above 0.0765, twice the iminuit value 0.0382"
    );
    // A value below 1e-6 would mean the reference came from the run and not from the closed form.
    assert!(
        error > 1e-6,
        "covariance error {error} is implausibly small"
    );

    // Hesse on the same minimum must still return a valid result.
    let hesse = MnHesse::new()
        .with_strategy(1)
        .calculate(&Rosenbrock, &minimum);
    assert!(hesse.is_valid(), "Hesse did not return a valid minimum");
}
