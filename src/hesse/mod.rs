//! MnHesse: accurate error analysis via full Hessian computation.
//!
//! After Migrad converges, MnHesse computes the full second-derivative matrix
//! (Hessian) to get accurate parameter errors and correlations. This is
//! especially important when Migrad's DFP-updated inverse Hessian is only
//! approximate.

pub mod calculator;
pub mod gradient;

use crate::application::default_max_fcn;
use crate::fcn::FCN;
use crate::global_cc::global_correlation_coefficients;
use crate::minimum::FunctionMinimum;
use crate::mn_fcn::MnFcn;
use crate::strategy::MnStrategy;
use crate::user_covariance::MnUserCovariance;
use crate::user_parameter_state::MnUserParameterState;

/// Builder for running Hesse error analysis.
pub struct MnHesse {
    strategy: MnStrategy,
    max_calls: Option<usize>,
}

impl MnHesse {
    /// Create a new Hesse error calculator with default settings.
    pub fn new() -> Self {
        Self {
            strategy: MnStrategy::default(),
            max_calls: None,
        }
    }

    /// Set strategy level (0=low, 1=medium, 2=high).
    pub fn with_strategy(mut self, level: u32) -> Self {
        self.strategy = MnStrategy::new(level);
        self
    }

    /// Set maximum number of function calls.
    pub fn with_max_calls(mut self, max: usize) -> Self {
        self.max_calls = Some(max);
        self
    }

    pub fn ncycles(&self) -> u32 {
        self.strategy.hessian_ncycles()
    }

    pub fn tolerstp(&self) -> f64 {
        self.strategy.hessian_step_tolerance()
    }

    pub fn toler_g2(&self) -> f64 {
        self.strategy.hessian_g2_tolerance()
    }

    /// Run Hesse on a minimization result.
    ///
    /// Returns a new FunctionMinimum with accurate covariance matrix.
    pub fn calculate(&self, fcn: &dyn FCN, minimum: &FunctionMinimum) -> FunctionMinimum {
        let trafo = minimum.seed().trafo();
        let n = trafo.variable_parameters();
        let maxcalls = self.max_calls.unwrap_or_else(|| default_max_fcn(n));

        let mn_fcn = MnFcn::new(fcn, trafo);
        let state = minimum.state();

        let result = calculator::calculate(&mn_fcn, state, trafo, &self.strategy, maxcalls);

        // Build new FunctionMinimum with the Hesse state
        let mut states = minimum.states().to_vec();
        states.push(result.state);

        if !states.last().is_some_and(|state| state.error().is_valid()) {
            return FunctionMinimum::above_max_edm(minimum.seed().clone(), states, minimum.up());
        }

        let mut min = FunctionMinimum::new(minimum.seed().clone(), states, minimum.up());
        // Update user state with covariance info
        let hesse_state = min.state();
        let user_state = build_user_state_with_covariance(
            minimum,
            hesse_state.error().matrix(),
            minimum.up(),
            trafo,
        );
        min.set_user_state(user_state);
        min
    }

    /// Compute errors and covariance without modifying the FunctionMinimum.
    ///
    /// Returns an MnUserParameterState with updated errors and covariance.
    pub fn calculate_errors(
        &self,
        fcn: &dyn FCN,
        minimum: &FunctionMinimum,
    ) -> MnUserParameterState {
        let trafo = minimum.seed().trafo();
        let n = trafo.variable_parameters();
        let maxcalls = self.max_calls.unwrap_or_else(|| default_max_fcn(n));

        let mn_fcn = MnFcn::new(fcn, trafo);
        let state = minimum.state();

        let result = calculator::calculate(&mn_fcn, state, trafo, &self.strategy, maxcalls);

        if !result.state.error().is_valid() {
            return minimum.user_state().clone();
        }

        build_user_state_with_covariance(
            minimum,
            result.state.error().matrix(),
            minimum.up(),
            trafo,
        )
    }
}

impl Default for MnHesse {
    fn default() -> Self {
        Self::new()
    }
}

/// Build an MnUserParameterState with covariance from the internal error matrix.
fn build_user_state_with_covariance(
    minimum: &FunctionMinimum,
    internal_cov: &nalgebra::DMatrix<f64>,
    up: f64,
    trafo: &crate::user_transformation::MnUserTransformation,
) -> MnUserParameterState {
    let mut user_state = minimum.user_state().clone();
    let n = trafo.variable_parameters();

    // Transform internal covariance to external covariance.
    // ROOT Minuit2 convention for user covariance is:
    //   V_user = 2 * up * V_int transformed to external coordinates.
    // where V_int is the internal inverse Hessian.
    let internal_params = minimum.state().parameters().vec();
    let mut ext_cov = MnUserCovariance::new(n);

    for i in 0..n {
        let ext_i = trafo.ext_of_int(i);
        let dxdi_i = trafo.dint2ext(ext_i, internal_params[i]);

        for j in i..n {
            let ext_j = trafo.ext_of_int(j);
            let dxdi_j = trafo.dint2ext(ext_j, internal_params[j]);

            let val = (2.0 * up) * dxdi_i * internal_cov[(i, j)] * dxdi_j;
            ext_cov.set(i, j, val);
        }
    }

    // Parameter errors are sqrt(diagonal(user covariance)).
    for i in 0..n {
        let ext_i = trafo.ext_of_int(i);
        let err = ext_cov.get(i, i).sqrt();
        user_state.set_error(ext_i, err);
    }

    // Set covariance
    user_state.set_covariance(ext_cov.clone());

    // Compute and set global correlation coefficients
    let mut cov_mat = nalgebra::DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            cov_mat[(i, j)] = ext_cov.get(i, j);
        }
    }
    let (gcc, _) = global_correlation_coefficients(&cov_mat);
    user_state.set_global_cc(gcc);

    user_state
}
