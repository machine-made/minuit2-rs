//! Public Migrad (variable-metric) minimizer API.
//!
//! `MnMigrad` is the user-facing entry point for quasi-Newton minimization
//! using the Davidon-Fletcher-Powell (DFP) update of the inverse Hessian.
//! Uses a builder pattern to configure parameters, then call `minimize()`.

pub mod builder;
pub mod minimizer;
pub mod seed;

use crate::application::{DEFAULT_TOLERANCE, default_max_fcn};
use crate::fcn::{FCN, FCNGradient};
use crate::minimum::FunctionMinimum;
use crate::mn_fcn::MnFcn;
use crate::strategy::MnStrategy;
use crate::user_parameters::MnUserParameters;

/// Builder for configuring and running Migrad minimization.
pub struct MnMigrad {
    params: MnUserParameters,
    strategy: MnStrategy,
    max_fcn: Option<usize>,
    tolerance: f64,
}

impl MnMigrad {
    /// Create a new Migrad minimizer with default settings.
    pub fn new() -> Self {
        Self {
            params: MnUserParameters::new(),
            strategy: MnStrategy::default(),
            max_fcn: None,
            tolerance: DEFAULT_TOLERANCE,
        }
    }

    /// Set the optimization strategy level.
    pub fn with_strategy(mut self, level: u32) -> Self {
        self.strategy = MnStrategy::new(level);
        self
    }

    /// Add a free parameter.
    pub fn add(mut self, name: impl Into<String>, value: f64, error: f64) -> Self {
        self.params.add(name, value, error);
        self
    }

    /// Add a parameter with both bounds.
    pub fn add_limited(
        mut self,
        name: impl Into<String>,
        value: f64,
        error: f64,
        lower: f64,
        upper: f64,
    ) -> Self {
        self.params.add_limited(name, value, error, lower, upper);
        self
    }

    /// Add a parameter with lower bound only.
    pub fn add_lower_limited(
        mut self,
        name: impl Into<String>,
        value: f64,
        error: f64,
        lower: f64,
    ) -> Self {
        self.params.add_lower_limited(name, value, error, lower);
        self
    }

    /// Add a parameter with upper bound only.
    pub fn add_upper_limited(
        mut self,
        name: impl Into<String>,
        value: f64,
        error: f64,
        upper: f64,
    ) -> Self {
        self.params.add_upper_limited(name, value, error, upper);
        self
    }

    /// Add a constant parameter.
    pub fn add_const(mut self, name: impl Into<String>, value: f64) -> Self {
        self.params.add_const(name, value);
        self
    }

    /// Fix parameter by index.
    pub fn fix(mut self, ext: usize) -> Self {
        self.params.fix(ext);
        self
    }

    /// Set maximum number of function calls. Default = 200 + 100*n + 5*n^2.
    pub fn max_fcn(mut self, max: usize) -> Self {
        self.max_fcn = Some(max);
        self
    }

    /// Set tolerance (relative to error_def). Default = 0.1.
    pub fn tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self
    }

    /// Run the minimization with numerical gradients (default).
    pub fn minimize(&self, fcn: &dyn FCN) -> FunctionMinimum {
        let n = self.params.variable_parameters();
        let max_fcn = self.max_fcn.unwrap_or_else(|| default_max_fcn(n));
        let trafo = self.params.trafo().clone();

        let mn_fcn = MnFcn::new(fcn, &trafo);
        minimizer::VariableMetricMinimizer::minimize(
            &mn_fcn,
            &trafo,
            &self.strategy,
            max_fcn,
            self.tolerance,
        )
    }

    /// Run the minimization with user-provided analytical gradients.
    ///
    /// Uses the analytical gradients provided by `FCNGradient::gradient()`.
    /// This typically requires fewer function evaluations than numerical differentiation.
    pub fn minimize_grad(&self, fcn: &dyn FCNGradient) -> FunctionMinimum {
        let n = self.params.variable_parameters();
        let max_fcn = self.max_fcn.unwrap_or_else(|| default_max_fcn(n));
        let trafo = self.params.trafo().clone();

        minimizer::VariableMetricMinimizer::minimize_with_gradient(
            fcn,
            &trafo,
            &self.strategy,
            max_fcn,
            self.tolerance,
        )
    }
}

impl Default for MnMigrad {
    fn default() -> Self {
        Self::new()
    }
}
