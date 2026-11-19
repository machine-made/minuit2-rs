//! Covariance from the Hessian supplied by the function.

use nalgebra::{DMatrix, DVector};

use crate::fcn::FCN;
use crate::minimum::error::MinimumError;
use crate::minimum::gradient::FunctionGradient;
use crate::minimum::state::MinimumState;
use crate::precision::MnMachinePrecision;
use crate::user_transformation::MnUserTransformation;

use super::calculator::HesseResult;

/// Follow ROOT MnHesse.cxx:100-168 without evaluating the function.
pub(super) fn calculate(
    fcn: &dyn FCN,
    state: &MinimumState,
    trafo: &MnUserTransformation,
) -> HesseResult {
    let n = state.parameters().vec().len();
    let external = trafo.transform(state.parameters().vec().as_slice());
    let supplied = fcn.hessian(&external);
    let n_ext = external.len();

    // An empty or malformed Hessian follows ROOT's failed-Hessian return.
    if supplied.is_empty()
        || supplied.len() != n_ext * (n_ext + 1) / 2
        || !supplied.iter().all(|entry| entry.is_finite())
    {
        return hesse_failed(state);
    }

    // ROOT AnalyticalGradientCalculator.cxx:75-89 uses external parameters and J_i H_ij J_j.
    // ExternalInternalGradientCalculator inherits this Hessian method, including the transformation.
    // Rust stores the upper triangle by columns instead of ROOT's full matrix.
    let mut hessian = DMatrix::zeros(n, n);
    for i in 0..n {
        let ext_i = trafo.ext_of_int(i);
        let jac_i = trafo.dint2ext(ext_i, state.parameters().vec()[i]);
        for j in i..n {
            let ext_j = trafo.ext_of_int(j);
            let jac_j = trafo.dint2ext(ext_j, state.parameters().vec()[j]);
            let (row, col) = (ext_i.min(ext_j), ext_i.max(ext_j));
            let value = jac_i * supplied[col * (col + 1) / 2 + row] * jac_j;
            hessian[(i, j)] = value;
            hessian[(j, i)] = value;
        }
    }

    // Limits far from zero make a Jacobian large enough for a finite Hessian to overflow in the product.
    // ROOT does not test the products and would carry the infinity into g2 and the covariance.
    // This check is the only difference: it returns the failed-Hessian state instead.
    if !hessian.iter().all(|entry| entry.is_finite()) {
        return hesse_failed(state);
    }

    // ROOT MnHesse.cxx:123-134 saves g2 before making the Hessian positive definite.
    let g2 = hessian.diagonal();
    // ROOT MnHesse.cxx:128 builds this gradient from two vectors only (FunctionGradient.h:35-38).
    // The saved step is dropped there, and the step vector is empty. Zeros of length n are the nearest defined value.
    let mut gradient = FunctionGradient::new(
        state.gradient().grad().clone(),
        g2.clone(),
        DVector::zeros(n),
    );
    // ROOT FunctionGradient.h:35-36 keeps this two-vector gradient analytical.
    gradient.set_analytical(true);
    let (positive, made_pos_def) = make_pos_def(hessian, trafo.precision());

    // ROOT MnHesse.cxx:138-148 keeps the old EDM if inversion fails.
    let Some(covariance) = invert(positive) else {
        let mut error =
            MinimumError::new(DMatrix::from_diagonal(&g2.map(|value| 1.0 / value)), 1.0);
        error.set_invert_failed(true);
        return HesseResult {
            state: MinimumState::new(
                state.parameters().clone(),
                error,
                gradient,
                state.edm(),
                state.nfcn(),
            ),
            hesse_failed: false,
            invert_failed: true,
            made_pos_def: false,
        };
    };

    // ROOT MnHesse.cxx:151-168 estimates EDM for either successful matrix status.
    let mut error = MinimumError::new(covariance, if made_pos_def { 1.0 } else { 0.0 });
    error.set_made_pos_def(made_pos_def);
    let grad = gradient.grad();
    let edm = 0.5 * grad.dot(&(error.matrix() * grad));
    HesseResult {
        state: MinimumState::new(
            state.parameters().clone(),
            error,
            gradient,
            edm,
            state.nfcn(),
        ),
        hesse_failed: false,
        invert_failed: false,
        made_pos_def,
    }
}

/// Return the failed-Hessian state of ROOT MnHesse.cxx:118-122.
///
/// The matrix is zero. Parameters, gradient, EDM and call count are those of the input state.
fn hesse_failed(state: &MinimumState) -> HesseResult {
    let n = state.parameters().vec().len();
    let mut error = MinimumError::new(DMatrix::zeros(n, n), 1.0);
    error.set_hesse_failed(true);
    HesseResult {
        state: MinimumState::new(
            state.parameters().clone(),
            error,
            state.gradient().clone(),
            state.edm(),
            state.nfcn(),
        ),
        hesse_failed: true,
        invert_failed: false,
        made_pos_def: false,
    }
}

/// Invert a positive-definite matrix, scaling it by a power of two when its entries are far from one.
///
/// For two to four rows nalgebra divides by a determinant of about `largest^n` or `smallest^n`.
/// This overflows or underflows for large or small entries. With two parameters and a diagonal of 1e200 the inverse is zero.
/// ROOT normalises the diagonal first (MnMatrix.cxx:33-41) and has no such limit.
/// Scaling by 2^k is exact, so the scaled inverse times 2^-k equals the true inverse up to normal rounding.
///
/// The scaled path is used when `n * log10(largest) > 290` or `n * log10(smallest nonzero diagonal) < -290`.
/// An `f64` overflows near 1e308 and has no full precision below 2.2e-308, so the bound keeps the determinant about 17 decades inside.
/// The sums in the cofactors add at most a factor of 24 for four rows.
/// Larger matrices use a pivoted elimination without a determinant, but the same bound is kept for all sizes.
/// Every other matrix takes the plain `try_inverse()` path, so its result does not change by a single bit.
#[expect(
    clippy::cast_precision_loss,
    reason = "the count is far below 2^53, so the conversion is exact"
)]
#[expect(
    clippy::cast_possible_truncation,
    reason = "the rounded binary exponent fits in i32"
)]
fn invert(matrix: DMatrix<f64>) -> Option<DMatrix<f64>> {
    let n = matrix.nrows() as f64;
    let largest = matrix
        .iter()
        .fold(0.0_f64, |acc, entry| acc.max(entry.abs()));
    let smallest_diagonal = matrix
        .diagonal()
        .iter()
        .map(|entry| entry.abs())
        .filter(|entry| *entry > 0.0)
        .fold(f64::INFINITY, f64::min);
    // A zero matrix or one without rows gives NaN or infinity here, and both comparisons are then false.
    let needs_scaling = n * largest.log10() > 290.0 || n * smallest_diagonal.log10() < -290.0;
    if !needs_scaling {
        return matrix.try_inverse();
    }

    // The limit of 1000 keeps both 2^k and 2^-k finite for any entry of `f64`.
    let k = (largest.log2().round() as i32).clamp(-1000, 1000);
    let scale = f64::powi(2.0, k);
    // When the scaling makes the small entries underflow, the plain inversion is still tried, as before.
    // An inverse that is not finite from either attempt is a failure, not a covariance of infinity.
    (&matrix / scale)
        .try_inverse()
        .map(|inverse| inverse / scale)
        .filter(|inverse| inverse.iter().all(|entry| entry.is_finite()))
        .or_else(|| matrix.try_inverse())
        .filter(|inverse| inverse.iter().all(|entry| entry.is_finite()))
}

/// Apply ROOT MnPosDef.cxx:35-102 to the supplied Hessian.
fn make_pos_def(mut matrix: DMatrix<f64>, precision: &MnMachinePrecision) -> (DMatrix<f64>, bool) {
    let n = matrix.nrows();
    if n == 0 {
        return (matrix, false);
    }
    if n == 1 {
        if matrix[(0, 0)] < precision.eps() {
            matrix[(0, 0)] = 1.0;
            return (matrix, true);
        }
        if matrix[(0, 0)] > precision.eps() {
            return (matrix, false);
        }
    }

    let epspdf = precision.eps2().max(1e-12);
    let minimum_diagonal = matrix.diagonal().min();
    let shift = if minimum_diagonal <= 0.0 {
        0.5 + epspdf - minimum_diagonal
    } else {
        0.0
    };
    for i in 0..n {
        matrix[(i, i)] += shift;
        if matrix[(i, i)] < 0.0 {
            matrix[(i, i)] = 1.0;
        }
    }
    let scale = matrix.diagonal().map(|value| 1.0 / value.sqrt());
    let normalized = DMatrix::from_fn(n, n, |i, j| matrix[(i, j)] * scale[i] * scale[j]);
    let eigenvalues = normalized.symmetric_eigen().eigenvalues;
    let pmin = eigenvalues.min();
    let pmax = eigenvalues.max().abs().max(1.0);
    if pmin > epspdf * pmax {
        // ROOT preserves the normal status when the initial diagonal shift suffices.
        return (matrix, false);
    }

    let addition = 0.001 * pmax - pmin;
    for i in 0..n {
        matrix[(i, i)] *= 1.0 + addition;
    }
    (matrix, true)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic numbers in [0, 1) from a linear congruential generator.
    #[expect(
        clippy::cast_precision_loss,
        reason = "shifting leaves at most 53 bits, so the conversion is exact"
    )]
    fn uniform(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (*state >> 11) as f64 / 9_007_199_254_740_992.0
    }

    /// Symmetric, diagonally dominant matrix with diagonal entries from 10^low to 10^high.
    #[expect(
        clippy::cast_precision_loss,
        reason = "the count is far below 2^53, so the conversion is exact"
    )]
    fn ordinary_matrix(n: usize, low: f64, high: f64, seed: u64) -> DMatrix<f64> {
        let mut state = seed;
        let diagonal: Vec<f64> = (0..n)
            .map(|_| 10f64.powf(low + (high - low) * uniform(&mut state)))
            .collect();
        let mut matrix = DMatrix::zeros(n, n);
        for i in 0..n {
            matrix[(i, i)] = diagonal[i];
            for j in (i + 1)..n {
                let value = 0.5 / n as f64
                    * (diagonal[i] * diagonal[j]).sqrt()
                    * (2.0 * uniform(&mut state) - 1.0);
                matrix[(i, j)] = value;
                matrix[(j, i)] = value;
            }
        }
        matrix
    }

    #[test]
    #[expect(
        clippy::cast_precision_loss,
        reason = "the count is far below 2^53, so the conversion is exact"
    )]
    fn ordinary_matrices_keep_the_plain_inverse_bit_for_bit() {
        for n in [1, 2, 3, 4, 5, 8, 13, 20] {
            for (low, high) in [(0.0, 10.0), (-5.0, 5.0)] {
                for seed in 1..=3 {
                    let matrix = ordinary_matrix(n, low, high, seed);
                    // These matrices must lie inside the plain path, or the comparison proves nothing.
                    let largest = matrix
                        .iter()
                        .fold(0.0_f64, |acc, entry| acc.max(entry.abs()));
                    let smallest = matrix
                        .diagonal()
                        .iter()
                        .fold(f64::INFINITY, |acc, entry| acc.min(entry.abs()));
                    assert!(n as f64 * largest.log10() <= 290.0);
                    assert!(n as f64 * smallest.log10() >= -290.0);

                    let direct = matrix.clone().try_inverse().expect("invertible");
                    let got = invert(matrix).expect("invertible");
                    // Exact: the plain path is the same call, so the bits must agree.
                    for (a, b) in got.iter().zip(direct.iter()) {
                        assert_eq!(a.to_bits(), b.to_bits(), "n = {n}, seed = {seed}");
                    }
                }
            }
        }
    }
}
