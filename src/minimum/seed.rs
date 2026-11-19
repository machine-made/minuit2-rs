//! Initial state for a minimizer: starting point + gradient + error matrix.
//!
//! The seed is constructed by a seed generator and passed to a minimum builder.

use super::state::MinimumState;
use crate::user_transformation::MnUserTransformation;

/// Starting state and coordinate transformation supplied to a minimizer.
/// # Examples
///
/// ```
/// use minuit2::{MinuitParameter, MnStrategy, MnUserTransformation};
/// use minuit2::mn_fcn::MnFcn;
/// let trafo = MnUserTransformation::new(vec![MinuitParameter::new(0, "x", 0.0, 0.1)]);
/// let f = |p: &[f64]| (p[0] - 1.0).powi(2);
/// let fcn = MnFcn::new(&f, &trafo);
/// let strategy = MnStrategy::default();
/// let seed = minuit2::migrad::seed::MigradSeedGenerator::generate(&fcn, &trafo, &strategy);
/// assert!(seed.is_valid());
/// assert_eq!(seed.n_variable_params(), 1);
/// ```
#[derive(Debug, Clone)]
pub struct MinimumSeed {
    state: MinimumState,
    trafo: MnUserTransformation,
    is_valid: bool,
}

impl MinimumSeed {
    /// Store the starting state and its parameter transformation.
    #[must_use]
    pub fn new(state: MinimumState, trafo: MnUserTransformation) -> Self {
        let valid = state.is_valid();
        Self {
            state,
            trafo,
            is_valid: valid,
        }
    }

    /// Return the stored minimization state.
    #[must_use]
    pub fn state(&self) -> &MinimumState {
        &self.state
    }

    /// Return the parameter transformation.
    #[must_use]
    pub fn trafo(&self) -> &MnUserTransformation {
        &self.trafo
    }

    /// Return the starting internal parameters.
    #[must_use]
    pub fn parameters(&self) -> &super::parameters::MinimumParameters {
        self.state.parameters()
    }

    /// Return the starting internal inverse Hessian.
    #[must_use]
    pub fn error(&self) -> &super::error::MinimumError {
        self.state.error()
    }

    /// Return the starting internal derivatives.
    #[must_use]
    pub fn gradient(&self) -> &super::gradient::FunctionGradient {
        self.state.gradient()
    }

    /// Return the objective value.
    #[must_use]
    pub fn fval(&self) -> f64 {
        self.state.fval()
    }

    /// Return the estimated distance to the minimum, in objective units.
    #[must_use]
    pub fn edm(&self) -> f64 {
        self.state.edm()
    }

    /// Return the number of objective evaluations.
    #[must_use]
    pub fn nfcn(&self) -> usize {
        self.state.nfcn()
    }

    /// Return the precision used for coordinate transformations.
    #[must_use]
    pub fn precision(&self) -> &crate::precision::MnMachinePrecision {
        self.trafo.precision()
    }

    /// Return the validity captured from the starting state.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.is_valid
    }

    /// Return the number of variable parameters.
    #[must_use]
    pub fn n_variable_params(&self) -> usize {
        self.trafo.variable_parameters()
    }
}
