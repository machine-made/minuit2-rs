//! Public Minimize minimizer API (hybrid Simplex + Migrad).
//!
//! `MnMinimize` is a combined minimizer that uses a two-phase approach:
//! 1. Runs Simplex (derivative-free) to locate the approximate minimum
//! 2. Then runs Migrad (variable-metric) from that point for precise convergence
//!
//! This hybrid approach is robust for difficult functions and has fast convergence near the minimum.
//! Uses a builder pattern to configure parameters, then call `minimize()`.

use crate::application::{DEFAULT_TOLERANCE, default_max_fcn};
use crate::fcn::FCN;
use crate::migrad::MnMigrad;
use crate::minimum::FunctionMinimum;
use crate::parameter::MinuitParameter;
use crate::simplex::MnSimplex;
use crate::strategy::MnStrategy;
use crate::user_parameters::MnUserParameters;

/// Builder for configuring and running hybrid Simplex+Migrad minimization.
pub struct MnMinimize {
    params: MnUserParameters,
    strategy: MnStrategy,
    max_fcn: Option<usize>,
    tolerance: f64,
}

impl MnMinimize {
    /// Create a new Minimize minimizer with default strategy.
    pub fn new() -> Self {
        Self {
            params: MnUserParameters::new(),
            strategy: MnStrategy::default(),
            max_fcn: None,
            tolerance: DEFAULT_TOLERANCE,
        }
    }

    /// Set strategy level (0=low, 1=medium, 2=high).
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

    fn configure_simplex_from_params(simplex: MnSimplex, params: &MnUserParameters) -> MnSimplex {
        configure_builder_from_params(simplex, params)
    }

    fn configure_migrad_from_params(migrad: MnMigrad, params: &MnUserParameters) -> MnMigrad {
        configure_builder_from_params(migrad, params)
    }

    /// Run the combined minimization.
    ///
    /// Combined minimizer sequence:
    /// 1) Try Migrad first.
    /// 2) If Migrad fails, run Simplex with strategy 2.
    /// 3) If Simplex succeeds, run Migrad again from that point (strategy 2).
    /// 4) If second Migrad fails, return the Simplex minimum.
    pub fn minimize(&self, fcn: &dyn FCN) -> FunctionMinimum {
        let n = self.params.variable_parameters();
        let max_fcn = self.max_fcn.unwrap_or_else(|| default_max_fcn(n));

        // Attempt 1: Migrad with user-selected strategy.
        let migrad = Self::configure_migrad_from_params(
            MnMigrad::new().with_strategy(self.strategy.strategy()),
            &self.params,
        )
        .max_fcn(max_fcn)
        .tolerance(self.tolerance);
        let min = migrad.minimize(fcn);

        if min.is_valid() {
            return min;
        }

        // Fallback path: use the more robust strategy level.
        let fallback_strategy = 2_u32;
        let simplex = Self::configure_simplex_from_params(
            MnSimplex::new().with_strategy(fallback_strategy),
            &self.params,
        )
        .max_fcn(max_fcn)
        .tolerance(self.tolerance);
        let simplex_min = simplex.minimize(fcn);

        if !simplex_min.is_valid() {
            return simplex_min;
        }

        let migrad2 = Self::configure_migrad_from_params(
            MnMigrad::new().with_strategy(fallback_strategy),
            simplex_min.user_state().params(),
        )
        .max_fcn(max_fcn)
        .tolerance(self.tolerance);
        let min2 = migrad2.minimize(fcn);

        if min2.is_valid() { min2 } else { simplex_min }
    }
}

impl Default for MnMinimize {
    fn default() -> Self {
        Self::new()
    }
}

trait ParameterBuilder: Sized {
    fn add_const(self, name: impl Into<String>, value: f64) -> Self;
    fn add_limited(
        self,
        name: impl Into<String>,
        value: f64,
        error: f64,
        lower: f64,
        upper: f64,
    ) -> Self;
    fn add_lower_limited(self, name: impl Into<String>, value: f64, error: f64, lower: f64)
    -> Self;
    fn add_upper_limited(self, name: impl Into<String>, value: f64, error: f64, upper: f64)
    -> Self;
    fn add(self, name: impl Into<String>, value: f64, error: f64) -> Self;
    fn fix(self, ext: usize) -> Self;

    fn add_from_param(self, param: &MinuitParameter) -> Self {
        let configured = if param.is_const() {
            self.add_const(param.name(), param.value())
        } else if param.has_limits() {
            self.add_limited(
                param.name(),
                param.value(),
                param.error(),
                param.lower_limit(),
                param.upper_limit(),
            )
        } else if param.has_lower_limit() {
            self.add_lower_limited(
                param.name(),
                param.value(),
                param.error(),
                param.lower_limit(),
            )
        } else if param.has_upper_limit() {
            self.add_upper_limited(
                param.name(),
                param.value(),
                param.error(),
                param.upper_limit(),
            )
        } else {
            self.add(param.name(), param.value(), param.error())
        };

        if param.is_fixed() && !param.is_const() {
            configured.fix(param.number())
        } else {
            configured
        }
    }
}

impl ParameterBuilder for MnSimplex {
    fn add_const(self, name: impl Into<String>, value: f64) -> Self {
        MnSimplex::add_const(self, name, value)
    }

    fn add_limited(
        self,
        name: impl Into<String>,
        value: f64,
        error: f64,
        lower: f64,
        upper: f64,
    ) -> Self {
        MnSimplex::add_limited(self, name, value, error, lower, upper)
    }

    fn add_lower_limited(
        self,
        name: impl Into<String>,
        value: f64,
        error: f64,
        lower: f64,
    ) -> Self {
        MnSimplex::add_lower_limited(self, name, value, error, lower)
    }

    fn add_upper_limited(
        self,
        name: impl Into<String>,
        value: f64,
        error: f64,
        upper: f64,
    ) -> Self {
        MnSimplex::add_upper_limited(self, name, value, error, upper)
    }

    fn add(self, name: impl Into<String>, value: f64, error: f64) -> Self {
        MnSimplex::add(self, name, value, error)
    }

    fn fix(self, ext: usize) -> Self {
        MnSimplex::fix(self, ext)
    }
}

impl ParameterBuilder for MnMigrad {
    fn add_const(self, name: impl Into<String>, value: f64) -> Self {
        MnMigrad::add_const(self, name, value)
    }

    fn add_limited(
        self,
        name: impl Into<String>,
        value: f64,
        error: f64,
        lower: f64,
        upper: f64,
    ) -> Self {
        MnMigrad::add_limited(self, name, value, error, lower, upper)
    }

    fn add_lower_limited(
        self,
        name: impl Into<String>,
        value: f64,
        error: f64,
        lower: f64,
    ) -> Self {
        MnMigrad::add_lower_limited(self, name, value, error, lower)
    }

    fn add_upper_limited(
        self,
        name: impl Into<String>,
        value: f64,
        error: f64,
        upper: f64,
    ) -> Self {
        MnMigrad::add_upper_limited(self, name, value, error, upper)
    }

    fn add(self, name: impl Into<String>, value: f64, error: f64) -> Self {
        MnMigrad::add(self, name, value, error)
    }

    fn fix(self, ext: usize) -> Self {
        MnMigrad::fix(self, ext)
    }
}

fn configure_builder_from_params<B: ParameterBuilder>(
    mut builder: B,
    params: &MnUserParameters,
) -> B {
    for param in params.params() {
        builder = builder.add_from_param(param);
    }
    builder
}
