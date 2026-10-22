//! Global correlation coefficients from a covariance matrix.
//!
//! Computes `gcc(i) = sqrt(1 - 1/(V_inv(i,i) * V(i,i)))`, which measures how
//! correlated each parameter is with all others combined.

use nalgebra::DMatrix;

/// Compute global correlation coefficients from a covariance matrix.
///
/// Returns `(coefficients, is_valid)`. If the matrix cannot be inverted,
/// returns zeros and `is_valid = false`.
pub fn global_correlation_coefficients(cov: &DMatrix<f64>) -> (Vec<f64>, bool) {
    let n = cov.nrows();
    assert_eq!(n, cov.ncols(), "covariance matrix must be square");

    let Some(inv) = cov.clone().try_inverse() else {
        return (vec![0.0; n], false);
    };

    let mut gcc = Vec::with_capacity(n);
    let mut valid = true;

    for i in 0..n {
        let denom = inv[(i, i)] * cov[(i, i)];
        if denom < 1.0 {
            // Rounding: correlation should be 0 if denom < 1
            gcc.push(0.0);
        } else {
            gcc.push((1.0 - 1.0 / denom).sqrt());
        }
    }

    // Check for NaN (shouldn't happen with the guard, but be safe)
    if gcc.iter().any(|&g| g.is_nan()) {
        valid = false;
    }

    (gcc, valid)
}

pub fn is_valid(cov: &DMatrix<f64>) -> bool {
    global_correlation_coefficients(cov).1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagonal_has_zero_correlation() {
        // Diagonal covariance → no correlations
        let cov = DMatrix::from_diagonal_element(3, 3, 1.0);
        let (gcc, valid) = global_correlation_coefficients(&cov);
        assert!(valid);
        for &g in &gcc {
            assert!(g.abs() < 1e-12, "diagonal should give gcc=0, got {g}");
        }
    }

    #[test]
    fn correlated_2x2() {
        // V = [[1, 0.8], [0.8, 1]] → known correlation
        let mut cov = DMatrix::zeros(2, 2);
        cov[(0, 0)] = 1.0;
        cov[(0, 1)] = 0.8;
        cov[(1, 0)] = 0.8;
        cov[(1, 1)] = 1.0;

        let (gcc, valid) = global_correlation_coefficients(&cov);
        assert!(valid);
        // For 2x2, gcc(i) = |rho| = 0.8
        assert!(
            (gcc[0] - 0.8).abs() < 1e-10,
            "gcc[0] should be 0.8, got {}",
            gcc[0]
        );
        assert!(
            (gcc[1] - 0.8).abs() < 1e-10,
            "gcc[1] should be 0.8, got {}",
            gcc[1]
        );
    }

    #[test]
    fn singular_matrix_returns_invalid() {
        let cov = DMatrix::zeros(2, 2);
        let (_, valid) = global_correlation_coefficients(&cov);
        assert!(!valid);
    }
}
