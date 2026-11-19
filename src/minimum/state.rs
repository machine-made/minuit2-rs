//! State of the minimization at a particular iteration.
//!
//! Combines parameters, error matrix, and gradient at one point in the
//! iteration history.

use super::error::MinimumError;
use super::gradient::FunctionGradient;
use super::parameters::MinimumParameters;

#[derive(Debug, Clone)]
pub struct MinimumState {
    parameters: MinimumParameters,
    error: MinimumError,
    gradient: FunctionGradient,
    edm: f64,
    nfcn: usize,
}

impl MinimumState {
    #[must_use]
    pub fn new(
        parameters: MinimumParameters,
        error: MinimumError,
        gradient: FunctionGradient,
        edm: f64,
        nfcn: usize,
    ) -> Self {
        Self {
            parameters,
            error,
            gradient,
            edm,
            nfcn,
        }
    }

    /// Create a minimal state with just parameters (no gradient/error info).
    /// Used by Simplex which doesn't compute a Hessian.
    #[must_use]
    pub fn from_params_edm(parameters: MinimumParameters, edm: f64, nfcn: usize) -> Self {
        let n = parameters.vec().len();
        Self {
            edm,
            nfcn,
            parameters,
            error: MinimumError::from_diagonal(n),
            gradient: FunctionGradient::new(
                nalgebra::DVector::zeros(n),
                nalgebra::DVector::zeros(n),
                nalgebra::DVector::zeros(n),
            ),
        }
    }

    /// Get the parameters at this state.
    #[must_use]
    pub fn parameters(&self) -> &MinimumParameters {
        &self.parameters
    }

    /// Get the error matrix at this state.
    #[must_use]
    pub fn error(&self) -> &MinimumError {
        &self.error
    }

    /// Get the gradient at this state.
    #[must_use]
    pub fn gradient(&self) -> &FunctionGradient {
        &self.gradient
    }

    #[must_use]
    pub fn fval(&self) -> f64 {
        self.parameters.fval()
    }

    #[must_use]
    pub fn edm(&self) -> f64 {
        self.edm
    }

    #[must_use]
    pub fn nfcn(&self) -> usize {
        self.nfcn
    }

    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.parameters.is_valid()
    }

    #[must_use]
    pub fn has_parameters(&self) -> bool {
        true
    }

    #[must_use]
    pub fn has_covariance(&self) -> bool {
        self.error.is_available()
    }

    #[must_use]
    pub fn size(&self) -> usize {
        self.parameters.vec().len()
    }
}
