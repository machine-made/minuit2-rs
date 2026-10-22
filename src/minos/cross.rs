//! MnCross: result of a single crossing-point search.
//!
//! Contains the crossing multiplier, user parameter state, number of function
//! calls, and status flags.

use crate::user_parameter_state::MnUserParameterState;

/// Result of searching for one crossing point (upper or lower).
#[derive(Debug, Clone)]
pub struct MnCross {
    /// Crossing parameter multiplier.
    value: f64,
    /// Parameter state at the crossing.
    state: MnUserParameterState,
    /// Number of function calls used.
    nfcn: usize,
    /// Whether the crossing was found successfully.
    valid: bool,
    /// Whether the crossing is at a parameter limit.
    is_at_limit: bool,
    /// Whether the function call limit was reached.
    is_at_max_fcn: bool,
    /// Whether a new minimum was found during the search.
    new_minimum: bool,
}

impl MnCross {
    /// Successful crossing result.
    pub fn valid(value: f64, state: MnUserParameterState, nfcn: usize) -> Self {
        Self {
            value,
            state,
            nfcn,
            valid: true,
            is_at_limit: false,
            is_at_max_fcn: false,
            new_minimum: false,
        }
    }

    /// Crossing at a parameter limit.
    pub fn limit_reached(nfcn: usize) -> Self {
        Self {
            value: 0.0,
            state: MnUserParameterState::new(crate::user_parameters::MnUserParameters::new()),
            nfcn,
            valid: false,
            is_at_limit: true,
            is_at_max_fcn: false,
            new_minimum: false,
        }
    }

    /// Call limit reached.
    pub fn call_limit_reached(nfcn: usize) -> Self {
        Self {
            value: 0.0,
            state: MnUserParameterState::new(crate::user_parameters::MnUserParameters::new()),
            nfcn,
            valid: false,
            is_at_limit: false,
            is_at_max_fcn: true,
            new_minimum: false,
        }
    }

    /// A new minimum was found (original minimum is no longer valid).
    pub fn new_minimum_found(state: MnUserParameterState, nfcn: usize) -> Self {
        Self {
            value: 0.0,
            state,
            nfcn,
            valid: false,
            is_at_limit: false,
            is_at_max_fcn: false,
            new_minimum: true,
        }
    }

    /// Invalid result (generic failure).
    pub fn invalid(nfcn: usize) -> Self {
        Self {
            value: 0.0,
            state: MnUserParameterState::new(crate::user_parameters::MnUserParameters::new()),
            nfcn,
            valid: false,
            is_at_limit: false,
            is_at_max_fcn: false,
            new_minimum: false,
        }
    }

    /// Get the crossing value (multiplier).
    pub fn value(&self) -> f64 {
        self.value
    }

    pub fn state(&self) -> &MnUserParameterState {
        &self.state
    }

    pub fn nfcn(&self) -> usize {
        self.nfcn
    }

    /// Check if the crossing was found successfully.
    pub fn is_valid(&self) -> bool {
        self.valid
    }

    /// Check if the crossing is at a parameter limit.
    pub fn at_limit(&self) -> bool {
        self.is_at_limit
    }

    pub fn at_max_fcn(&self) -> bool {
        self.is_at_max_fcn
    }

    /// Check if a new minimum was found during the search.
    pub fn new_minimum(&self) -> bool {
        self.new_minimum
    }
}
