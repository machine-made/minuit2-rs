//! Force a symmetric matrix to be positive-definite.
//!
//! Uses eigendecomposition to detect and correct non-positive-definite error
//! matrices by shifting eigenvalues.

use crate::precision::MnMachinePrecision;
use nalgebra::DMatrix;

/// Force `mat` to be positive-definite by shifting eigenvalues if needed.
///
/// Returns `(corrected_matrix, was_modified)`. The correction adds a small
/// amount to the diagonal until all eigenvalues are safely positive.
pub fn make_pos_def(mat: &DMatrix<f64>, prec: &MnMachinePrecision) -> (DMatrix<f64>, bool) {
    let n = mat.nrows();
    assert_eq!(n, mat.ncols(), "matrix must be square");

    if n == 0 {
        return (mat.clone(), false);
    }

    let epspdf = prec.eps2().max(1.0e-6);

    // Check diagonal elements first
    let mut dgmin = mat[(0, 0)];
    for i in 1..n {
        if mat[(i, i)] < dgmin {
            dgmin = mat[(i, i)];
        }
    }

    let mut err = mat.clone();
    let mut modified = false;

    // If minimum diagonal ≤ 0, shift all diagonals
    if dgmin <= 0.0 {
        let dg = 0.5 + epspdf - dgmin;
        for i in 0..n {
            err[(i, i)] += dg;
        }
        modified = true;
    }

    // Normalize: build correlation-like matrix p
    // s(i) = 1/sqrt(err(i,i))
    let mut s = vec![0.0; n];
    for i in 0..n {
        if err[(i, i)] > 0.0 {
            s[i] = 1.0 / err[(i, i)].sqrt();
        } else {
            s[i] = 1.0;
        }
    }

    let mut p = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            p[(i, j)] = err[(i, j)] * s[i] * s[j];
        }
    }

    // Eigendecomposition of the normalized matrix
    let eigen = p.symmetric_eigen();
    let eigenvalues = &eigen.eigenvalues;

    let mut pmin = eigenvalues[0];
    let mut pmax = eigenvalues[0].abs();
    for i in 1..n {
        if eigenvalues[i] < pmin {
            pmin = eigenvalues[i];
        }
        if eigenvalues[i].abs() > pmax {
            pmax = eigenvalues[i].abs();
        }
    }
    pmax = pmax.max(1.0);

    // If already positive-definite enough after any diagonal shift, return
    if pmin > epspdf * pmax {
        if modified {
            // Diagonal shift was applied — return the shifted matrix
            return (err, true);
        }
        return (mat.clone(), false);
    }

    // Shift: add padd to diagonal of eigenvalue matrix
    let padd = 0.001 * pmax - pmin;

    // Reconstruct: p_corrected = Q * diag(eigenvalues + padd) * Q^T
    // Then un-normalize back to err scale
    let q = &eigen.eigenvectors;
    let mut d = DMatrix::zeros(n, n);
    for i in 0..n {
        d[(i, i)] = eigenvalues[i] + padd;
    }

    let p_corrected = q * d * q.transpose();

    // Un-normalize: err(i,j) = p(i,j) / (s(i) * s(j))
    let mut result = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            result[(i, j)] = p_corrected[(i, j)] / (s[i] * s[j]);
        }
    }

    (result, true)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn already_posdef_unchanged() {
        // Identity matrix is already positive-definite
        let m = DMatrix::identity(3, 3);
        let prec = MnMachinePrecision::new();
        let (result, was_modified) = make_pos_def(&m, &prec);
        assert!(!was_modified);
        for i in 0..3 {
            for j in 0..3 {
                assert!((result[(i, j)] - m[(i, j)]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn non_posdef_gets_fixed() {
        // Matrix with a negative diagonal
        let mut m = DMatrix::identity(3, 3);
        m[(0, 0)] = -1.0;
        m[(0, 1)] = 0.5;
        m[(1, 0)] = 0.5;
        let prec = MnMachinePrecision::new();
        let (result, was_modified) = make_pos_def(&m, &prec);
        assert!(was_modified);

        // Check all eigenvalues of result are positive
        let eigen = result.symmetric_eigen();
        for ev in eigen.eigenvalues.iter() {
            assert!(*ev > 0.0, "eigenvalue {} should be positive", ev);
        }
    }

    #[test]
    fn diagonal_matrix_with_zero() {
        let mut m = DMatrix::zeros(2, 2);
        m[(0, 0)] = 1.0;
        m[(1, 1)] = 0.0; // zero eigenvalue = not positive-definite
        let prec = MnMachinePrecision::new();
        let (result, was_modified) = make_pos_def(&m, &prec);
        assert!(was_modified);

        let eigen = result.symmetric_eigen();
        for ev in eigen.eigenvalues.iter() {
            assert!(*ev > 0.0, "eigenvalue {} should be positive", ev);
        }
    }
}
