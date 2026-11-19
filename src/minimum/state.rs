//! State of the minimization at a particular iteration.
//!
//! Combines parameters, error matrix, and gradient at one point in the
//! iteration history.

use super::error::MinimumError;
use super::gradient::FunctionGradient;
use super::parameters::MinimumParameters;

/// Parameters, derivatives, inverse Hessian, objective distance, and call count at one iteration.
/// # Examples
///
/// ```
/// use minuit2::minimum::{parameters::MinimumParameters, state::MinimumState};
/// use nalgebra::DVector;
/// let p = MinimumParameters::new(DVector::from_vec(vec![0.0]), 0.0);
/// let state = MinimumState::from_params_edm(p, 0.01, 3);
/// assert_eq!(state.nfcn(), 3);
/// assert_eq!(state.size(), 1);
/// ```
#[derive(Debug, Clone)]
pub struct MinimumState {
    parameters: MinimumParameters,
    error: MinimumError,
    gradient: FunctionGradient,
    edm: f64,
    nfcn: usize,
}

impl MinimumState {
    /// Combine internal parameters, inverse Hessian, gradient, objective distance, and cumulative call count.
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

    /// Create a state with zero derivative vectors and an approximate identity inverse Hessian.
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

    /// Return the objective value.
    #[must_use]
    pub fn fval(&self) -> f64 {
        self.parameters.fval()
    }

    /// Return the estimated distance to the minimum, in objective units.
    #[must_use]
    pub fn edm(&self) -> f64 {
        self.edm
    }

    /// Return the number of objective evaluations.
    #[must_use]
    pub fn nfcn(&self) -> usize {
        self.nfcn
    }

    /// Return whether the stored parameters are valid.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.parameters.is_valid()
    }

    /// Return true; every state stores a parameter object.
    #[must_use]
    pub fn has_parameters(&self) -> bool {
        true
    }

    /// Return whether an inverse Hessian is available.
    #[must_use]
    pub fn has_covariance(&self) -> bool {
        self.error.is_available()
    }

    /// Return the number of internal parameters.
    #[must_use]
    pub fn size(&self) -> usize {
        self.parameters.vec().len()
    }
}
