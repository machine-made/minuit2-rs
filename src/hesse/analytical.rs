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
    let mut hessian = DMatrix::zeros(n, n);

    // An empty or malformed Hessian follows ROOT's failed-Hessian return.
    if supplied.is_empty()
        || supplied.len() != n_ext * (n_ext + 1) / 2
        || !supplied.iter().all(|entry| entry.is_finite())
    {
        let mut error = MinimumError::new(hessian, 1.0);
        error.set_hesse_failed(true);
        return HesseResult {
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
        };
    }

    // ROOT AnalyticalGradientCalculator.cxx:75-89 uses external parameters and J_i H_ij J_j.
    // ExternalInternalGradientCalculator inherits this Hessian method, including the transformation.
    // Rust stores the upper triangle by columns instead of ROOT's full matrix.
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

    // ROOT MnHesse.cxx:123-134 saves g2 before making the Hessian positive definite.
    let g2 = hessian.diagonal();
    let mut gradient = FunctionGradient::new(
        state.gradient().grad().clone(),
        g2.clone(),
        DVector::zeros(n),
    );
    // ROOT FunctionGradient.h:35-36 keeps this two-vector gradient analytical.
    gradient.set_analytical(true);
    let (positive, made_pos_def) = make_pos_def(hessian, trafo.precision());

    // ROOT MnHesse.cxx:138-148 keeps the old EDM if inversion fails.
    let Some(covariance) = positive.try_inverse() else {
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
