//! Initial state for a minimizer: starting point + gradient + error matrix.
//!
//! The seed is constructed by a seed generator and passed to a minimum builder.

use super::state::MinimumState;
use crate::user_transformation::MnUserTransformation;

#[derive(Debug, Clone)]
pub struct MinimumSeed {
    state: MinimumState,
    trafo: MnUserTransformation,
    is_valid: bool,
}

impl MinimumSeed {
    #[must_use]
    pub fn new(state: MinimumState, trafo: MnUserTransformation) -> Self {
        let valid = state.is_valid();
        Self {
            state,
            trafo,
            is_valid: valid,
        }
    }

    #[must_use]
    pub fn state(&self) -> &MinimumState {
        &self.state
    }

    #[must_use]
    pub fn trafo(&self) -> &MnUserTransformation {
        &self.trafo
    }

    #[must_use]
    pub fn parameters(&self) -> &super::parameters::MinimumParameters {
        self.state.parameters()
    }

    #[must_use]
    pub fn error(&self) -> &super::error::MinimumError {
        self.state.error()
    }

    #[must_use]
    pub fn gradient(&self) -> &super::gradient::FunctionGradient {
        self.state.gradient()
    }

    #[must_use]
    pub fn fval(&self) -> f64 {
        self.state.fval()
    }

    #[must_use]
    pub fn edm(&self) -> f64 {
        self.state.edm()
    }

    #[must_use]
    pub fn nfcn(&self) -> usize {
        self.state.nfcn()
    }

    #[must_use]
    pub fn precision(&self) -> &crate::precision::MnMachinePrecision {
        self.trafo.precision()
    }

    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.is_valid
    }

    #[must_use]
    pub fn n_variable_params(&self) -> usize {
        self.trafo.variable_parameters()
    }
}
