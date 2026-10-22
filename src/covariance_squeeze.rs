//! Remove a parameter from a covariance/Hessian matrix.
//!
//! Three covariance squeeze operations:
//! 1. Squeeze a DMatrix (remove row/col n)
//! 2. Squeeze an MnUserCovariance (invert → squeeze Hessian → invert back)
//! 3. Squeeze a MinimumError (extract Hessian → squeeze → invert → new MinimumError)

use nalgebra::DMatrix;

use crate::minimum::error::MinimumError;
use crate::user_covariance::MnUserCovariance;

/// Remove row and column `n` from a square matrix.
pub fn squeeze_matrix(mat: &DMatrix<f64>, n: usize) -> DMatrix<f64> {
    let dim = mat.nrows();
    assert!(n < dim, "index {n} out of range for {dim}x{dim} matrix");
    let new_dim = dim - 1;
    let mut result = DMatrix::zeros(new_dim, new_dim);

    let mut ri = 0;
    for i in 0..dim {
        if i == n {
            continue;
        }
        let mut rj = 0;
        for j in 0..dim {
            if j == n {
                continue;
            }
            result[(ri, rj)] = mat[(i, j)];
            rj += 1;
        }
        ri += 1;
    }
    result
}

/// Squeeze a user covariance matrix: invert → remove row/col → invert back.
///
/// If inversion fails, returns a diagonal matrix from the remaining elements.
pub fn squeeze_user_covariance(cov: &MnUserCovariance, n: usize) -> MnUserCovariance {
    let nrow = cov.nrow();
    assert!(
        n < nrow,
        "index {n} out of range for {nrow}x{nrow} covariance"
    );

    // Convert to DMatrix
    let mut mat = DMatrix::zeros(nrow, nrow);
    for i in 0..nrow {
        for j in 0..nrow {
            mat[(i, j)] = cov.get(i, j);
        }
    }

    // Invert to get Hessian
    let hessian = match mat.clone().try_inverse() {
        Some(h) => h,
        None => {
            // Fallback: return diagonal of remaining elements
            return make_diagonal_user_cov(cov, n);
        }
    };

    // Squeeze the Hessian
    let squeezed_h = squeeze_matrix(&hessian, n);

    // Invert back to covariance
    let squeezed_cov = match squeezed_h.try_inverse() {
        Some(c) => c,
        None => {
            return make_diagonal_user_cov(cov, n);
        }
    };

    // Convert back to MnUserCovariance
    let new_dim = nrow - 1;
    let mut result = MnUserCovariance::new(new_dim);
    for i in 0..new_dim {
        for j in i..new_dim {
            result.set(i, j, squeezed_cov[(i, j)]);
        }
    }
    result
}

/// Squeeze a MinimumError: extract Hessian → squeeze → invert → new error.
pub fn squeeze_error(err: &MinimumError, n: usize) -> MinimumError {
    let mat = err.matrix();
    let dim = mat.nrows();
    assert!(
        n < dim,
        "index {n} out of range for {dim}x{dim} error matrix"
    );

    // Get Hessian (inverse of error matrix)
    let hessian = match mat.clone().try_inverse() {
        Some(h) => h,
        None => {
            // Fallback: diagonal from squeezed covariance
            let squeezed = squeeze_matrix(mat, n);
            return MinimumError::new(squeezed, err.dcovar());
        }
    };

    // Squeeze the Hessian
    let squeezed_h = squeeze_matrix(&hessian, n);

    // Invert back to covariance
    match squeezed_h.try_inverse() {
        Some(cov) => MinimumError::new(cov, err.dcovar()),
        None => {
            // Fallback: use diagonal
            let new_dim = dim - 1;
            let mut diag = DMatrix::zeros(new_dim, new_dim);
            let mut ri = 0;
            for i in 0..dim {
                if i == n {
                    continue;
                }
                diag[(ri, ri)] = mat[(i, i)];
                ri += 1;
            }
            MinimumError::new(diag, err.dcovar())
        }
    }
}

/// Fallback: create diagonal MnUserCovariance from remaining elements.
fn make_diagonal_user_cov(cov: &MnUserCovariance, skip: usize) -> MnUserCovariance {
    let new_dim = cov.nrow() - 1;
    let mut result = MnUserCovariance::new(new_dim);
    let mut ri = 0;
    for i in 0..cov.nrow() {
        if i == skip {
            continue;
        }
        result.set(ri, ri, cov.get(i, i));
        ri += 1;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn squeeze_3x3_to_2x2() {
        // Remove middle row/col from:
        // [[1, 2, 3],
        //  [4, 5, 6],
        //  [7, 8, 9]]
        let mat = DMatrix::from_row_slice(3, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let squeezed = squeeze_matrix(&mat, 1);
        assert_eq!(squeezed.nrows(), 2);
        assert!((squeezed[(0, 0)] - 1.0).abs() < 1e-15);
        assert!((squeezed[(0, 1)] - 3.0).abs() < 1e-15);
        assert!((squeezed[(1, 0)] - 7.0).abs() < 1e-15);
        assert!((squeezed[(1, 1)] - 9.0).abs() < 1e-15);
    }

    #[test]
    fn squeeze_first_element() {
        let mat = DMatrix::from_row_slice(3, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let squeezed = squeeze_matrix(&mat, 0);
        assert_eq!(squeezed.nrows(), 2);
        assert!((squeezed[(0, 0)] - 5.0).abs() < 1e-15);
        assert!((squeezed[(0, 1)] - 6.0).abs() < 1e-15);
        assert!((squeezed[(1, 0)] - 8.0).abs() < 1e-15);
        assert!((squeezed[(1, 1)] - 9.0).abs() < 1e-15);
    }

    #[test]
    fn squeeze_error_preserves_structure() {
        // Create a valid 3x3 covariance (positive definite)
        let mut mat = DMatrix::zeros(3, 3);
        mat[(0, 0)] = 2.0;
        mat[(1, 1)] = 3.0;
        mat[(2, 2)] = 4.0;
        mat[(0, 1)] = 0.5;
        mat[(1, 0)] = 0.5;
        let err = MinimumError::new(mat, 0.0);

        let squeezed = squeeze_error(&err, 1);
        assert_eq!(squeezed.matrix().nrows(), 2);
    }
}
