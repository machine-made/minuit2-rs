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
    pub fn new(state: MinimumState, trafo: MnUserTransformation) -> Self {
        let valid = state.is_valid();
        Self {
            state,
            trafo,
            is_valid: valid,
        }
    }

    pub fn state(&self) -> &MinimumState {
        &self.state
    }

    pub fn trafo(&self) -> &MnUserTransformation {
        &self.trafo
    }

    pub fn parameters(&self) -> &super::parameters::MinimumParameters {
        self.state.parameters()
    }

    pub fn error(&self) -> &super::error::MinimumError {
        self.state.error()
    }

    pub fn gradient(&self) -> &super::gradient::FunctionGradient {
        self.state.gradient()
    }

    pub fn fval(&self) -> f64 {
        self.state.fval()
    }

    pub fn edm(&self) -> f64 {
        self.state.edm()
    }

    pub fn nfcn(&self) -> usize {
        self.state.nfcn()
    }

    pub fn precision(&self) -> &crate::precision::MnMachinePrecision {
        self.trafo.precision()
    }

    pub fn is_valid(&self) -> bool {
        self.is_valid
    }

    pub fn n_variable_params(&self) -> usize {
        self.trafo.variable_parameters()
    }
}
