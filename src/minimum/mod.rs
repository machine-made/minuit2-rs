//! Top-level minimization result.
//!
//! Contains the full iteration history, final state, user transformation, and
//! validity flags.

pub mod error;
pub mod gradient;
pub mod parameters;
pub mod seed;
pub mod state;

use seed::MinimumSeed;
use state::MinimumState;

use crate::global_cc::global_correlation_coefficients;
use crate::user_parameter_state::MnUserParameterState;
use crate::user_parameters::MnUserParameters;

/// Result of a minimization.
#[derive(Debug, Clone)]
pub struct FunctionMinimum {
    seed: MinimumSeed,
    states: Vec<MinimumState>,
    up: f64,
    is_above_max_edm: bool,
    reached_call_limit: bool,
    user_state: MnUserParameterState,
}

impl FunctionMinimum {
    pub fn new(seed: MinimumSeed, states: Vec<MinimumState>, up: f64) -> Self {
        // Build user state from the final internal state
        let user_state = Self::build_user_state(&seed, states.last().unwrap_or(seed.state()), up);

        Self {
            seed,
            states,
            up,
            is_above_max_edm: false,
            reached_call_limit: false,
            user_state,
        }
    }

    /// Create a result that hit the call limit.
    pub fn with_call_limit(seed: MinimumSeed, states: Vec<MinimumState>, up: f64) -> Self {
        let user_state = Self::build_user_state(&seed, states.last().unwrap_or(seed.state()), up);
        Self {
            seed,
            states,
            up,
            is_above_max_edm: false,
            reached_call_limit: true,
            user_state,
        }
    }

    /// Create a result above max EDM.
    pub fn above_max_edm(seed: MinimumSeed, states: Vec<MinimumState>, up: f64) -> Self {
        let user_state = Self::build_user_state(&seed, states.last().unwrap_or(seed.state()), up);
        Self {
            seed,
            states,
            up,
            is_above_max_edm: true,
            reached_call_limit: false,
            user_state,
        }
    }

    fn build_user_state(seed: &MinimumSeed, last: &MinimumState, up: f64) -> MnUserParameterState {
        let trafo = seed.trafo();
        let internal = last.parameters().vec().as_slice();
        let external = trafo.transform(internal);
        let cov_is_valid = last.error().is_valid();
        let mut ext_cov_opt = if cov_is_valid {
            let mut cov = trafo.int2ext_covariance(internal, last.error().matrix());
            for v in cov.data_mut().iter_mut() {
                *v *= 2.0 * up;
            }
            Some(cov)
        } else {
            None
        };

        // Build MnUserParameters with updated values
        let mut uparams = MnUserParameters::new();
        for (i, p) in trafo.parameters().iter().enumerate() {
            if p.is_const() {
                uparams.add_const(p.name(), p.value());
                continue;
            }

            if p.is_fixed() {
                Self::add_parameter_from_state(&mut uparams, p, p.value(), p.error());
                uparams.fix(i);
                continue;
            }

            let err = if cov_is_valid {
                Self::transformed_error(trafo, i, internal, last, up)
            } else {
                p.error()
            };
            Self::add_parameter_from_state(&mut uparams, p, external[i], err);
        }

        let mut state = MnUserParameterState::new(uparams);
        state.set_fval(last.fval());
        state.set_edm(last.edm());
        state.set_nfcn(last.nfcn());
        state.set_valid(true);

        if let Some(ext_cov) = ext_cov_opt.take() {
            state.set_covariance(ext_cov.clone());
            let n = ext_cov.nrow();
            let mut cov_mat = nalgebra::DMatrix::zeros(n, n);
            for i in 0..n {
                for j in 0..n {
                    cov_mat[(i, j)] = ext_cov.get(i, j);
                }
            }
            let (gcc, _) = global_correlation_coefficients(&cov_mat);
            state.set_global_cc(gcc);
        }
        state
    }

    fn add_parameter_from_state(
        params: &mut MnUserParameters,
        p: &crate::parameter::MinuitParameter,
        value: f64,
        error: f64,
    ) {
        if p.has_limits() {
            params.add_limited(p.name(), value, error, p.lower_limit(), p.upper_limit());
        } else if p.has_lower_limit() {
            params.add_lower_limited(p.name(), value, error, p.lower_limit());
        } else if p.has_upper_limit() {
            params.add_upper_limited(p.name(), value, error, p.upper_limit());
        } else {
            params.add(p.name(), value, error);
        }
    }

    fn transformed_error(
        trafo: &crate::user_transformation::MnUserTransformation,
        i: usize,
        internal: &[f64],
        last: &MinimumState,
        up: f64,
    ) -> f64 {
        let int_i = trafo
            .int_of_ext(i)
            .expect("variable parameter must map to internal index");
        let sigma_int = (2.0 * up * last.error().matrix()[(int_i, int_i)]).sqrt();
        trafo.int2ext_error(i, internal[int_i], sigma_int)
    }

    // --- Accessors ---

    /// Get the initial seed.
    pub fn seed(&self) -> &MinimumSeed {
        &self.seed
    }

    pub fn states(&self) -> &[MinimumState] {
        &self.states
    }

    /// Get the final state.
    pub fn state(&self) -> &MinimumState {
        self.states.last().unwrap_or_else(|| self.seed.state())
    }

    /// Get the user-facing state (external space).
    pub fn user_state(&self) -> &MnUserParameterState {
        &self.user_state
    }

    /// Function value at the minimum.
    pub fn fval(&self) -> f64 {
        self.state().fval()
    }

    /// Estimated distance to minimum.
    pub fn edm(&self) -> f64 {
        self.state().edm()
    }

    /// Total function calls.
    pub fn nfcn(&self) -> usize {
        self.state().nfcn()
    }

    /// Error definition (Up value): 1.0 for chi-square, 0.5 for likelihood.
    pub fn up(&self) -> f64 {
        self.up
    }

    /// Whether the minimization converged properly.
    pub fn is_valid(&self) -> bool {
        self.state().is_valid() && !self.is_above_max_edm && !self.reached_call_limit
    }

    pub fn has_valid_parameters(&self) -> bool {
        self.state().has_parameters()
    }

    pub fn has_made_pos_def_covar(&self) -> bool {
        self.state().error().is_made_pos_def()
    }

    /// Check if the result is above the maximum EDM threshold.
    pub fn is_above_max_edm(&self) -> bool {
        self.is_above_max_edm
    }

    /// Check if the function call limit was reached.
    pub fn reached_call_limit(&self) -> bool {
        self.reached_call_limit
    }

    /// Parameter values in external (user) space.
    pub fn params(&self) -> Vec<f64> {
        self.seed
            .trafo()
            .transform(self.state().parameters().vec().as_slice())
    }

    /// Number of variable parameters.
    pub fn n_variable_params(&self) -> usize {
        self.seed.n_variable_params()
    }

    pub fn set_error_def(&mut self, up: f64) {
        self.up = up;
        let rebuilt = Self::build_user_state(&self.seed, self.state(), up);
        self.user_state = rebuilt;
    }

    /// Replace the user state (used by Hesse to inject covariance info).
    pub fn set_user_state(&mut self, state: MnUserParameterState) {
        self.user_state = state;
    }
}
