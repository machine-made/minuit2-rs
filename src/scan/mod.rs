//! MnScan / MnParameterScan: 1D parameter scan.
//!
//! Evaluates the function along one parameter direction, keeping others at
//! their minimum values.

pub mod builder;
pub mod minimizer;
pub mod seed;

use crate::application::default_max_fcn;
use crate::fcn::FCN;
use crate::minimum::FunctionMinimum;
use crate::mn_fcn::MnFcn;
use crate::strategy::MnStrategy;
use crate::user_parameters::MnUserParameters;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Builder for ROOT Minuit2's brute-force SCAn minimizer.
///
/// SCAn evaluates a coarse one-dimensional grid around the starting point for
/// each variable parameter and keeps the best point found.  Use it as a robust
/// pre-pass for bad starting values, often before running `MnMigrad`.  SCAn does
/// not calculate parameter errors or a covariance matrix; run `MnHesse` after it
/// when uncertainties are needed.
pub struct MnScanMinimizer {
    params: MnUserParameters,
    strategy: MnStrategy,
    max_fcn: Option<usize>,
}

impl MnScanMinimizer {
    /// Create a new SCAn minimizer with default settings.
    pub fn new() -> Self {
        Self {
            params: MnUserParameters::new(),
            strategy: MnStrategy::default(),
            max_fcn: None,
        }
    }

    /// Set the optimization strategy level used for seed generation.
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

    /// Set maximum number of function calls. Default = 200 + 100*n + 5*n^2.
    pub fn max_fcn(mut self, max: usize) -> Self {
        self.max_fcn = Some(max);
        self
    }

    /// Run the SCAn minimization.
    pub fn minimize(&self, fcn: &dyn FCN) -> FunctionMinimum {
        let n = self.params.variable_parameters();
        let max_fcn = self.max_fcn.unwrap_or_else(|| default_max_fcn(n));
        let trafo = self.params.trafo().clone();
        let mn_fcn = MnFcn::new(fcn, &trafo);
        minimizer::ScanMinimizer::minimize(&mn_fcn, &trafo, &self.strategy, max_fcn)
    }
}

impl Default for MnScanMinimizer {
    fn default() -> Self {
        Self::new()
    }
}

/// Low-level 1D parameter scan.
pub struct MnParameterScan<'a, F: FCN + ?Sized> {
    fcn: &'a F,
    params: MnUserParameters,
    fval: f64,
}

impl<'a, F: FCN + ?Sized> MnParameterScan<'a, F> {
    pub fn new(fcn: &'a F, params: MnUserParameters, fval: f64) -> Self {
        Self { fcn, params, fval }
    }

    /// Scan parameter `par` over `nsteps` points between `low` and `high`.
    ///
    /// If `low == high == 0.0`, auto-range to +/- 2*error.
    /// Returns `Vec<(parameter_value, function_value)>`.
    /// Updates internal fval and param value if a better point is found.
    pub fn scan(&mut self, par: usize, nsteps: usize, low: f64, high: f64) -> Vec<(f64, f64)> {
        self.scan_serial(par, nsteps, low, high)
    }

    /// Serial implementation of 1D scan.
    pub fn scan_serial(
        &mut self,
        par: usize,
        nsteps: usize,
        low: f64,
        high: f64,
    ) -> Vec<(f64, f64)> {
        let (nsteps, low, high, values) = self.setup_scan(par, nsteps, low, high);
        let result = self.scan_points(par, nsteps, low, high, values.as_slice());
        self.update_best(par, &result);

        result
    }

    /// Parallel implementation of 1D scan (requires `parallel` feature).
    #[cfg(feature = "parallel")]
    pub fn scan_parallel(
        &mut self,
        par: usize,
        nsteps: usize,
        low: f64,
        high: f64,
    ) -> Vec<(f64, f64)>
    where
        F: Sync,
    {
        let (nsteps, low, high, values) = self.setup_scan(par, nsteps, low, high);
        let result = self.scan_points_parallel(par, nsteps, low, high, values.as_slice());
        self.update_best(par, &result);

        result
    }

    fn update_best(&mut self, par: usize, result: &[(f64, f64)]) {
        if let Some((x_best, f_best)) = result.iter().copied().min_by(|a, b| a.1.total_cmp(&b.1))
            && f_best < self.fval
        {
            self.fval = f_best;
            self.params.set_value(par, x_best);
        }
    }

    fn scan_points(
        &self,
        par: usize,
        nsteps: usize,
        low: f64,
        high: f64,
        values: &[f64],
    ) -> Vec<(f64, f64)> {
        let step = (high - low) / nsteps as f64;
        let mut pars = values.to_vec();
        (0..=nsteps)
            .map(|i| {
                let x = low + i as f64 * step;
                pars[par] = x;
                let f = self.fcn.value(&pars);
                (x, f)
            })
            .collect()
    }

    #[cfg(feature = "parallel")]
    fn scan_point(&self, par: usize, x: f64, values: &[f64]) -> (f64, f64) {
        let mut pars = values.to_vec();
        pars[par] = x;
        let f = self.fcn.value(&pars);
        (x, f)
    }

    #[cfg(feature = "parallel")]
    fn scan_points_parallel(
        &self,
        par: usize,
        nsteps: usize,
        low: f64,
        high: f64,
        values: &[f64],
    ) -> Vec<(f64, f64)>
    where
        F: Sync,
    {
        let step = (high - low) / nsteps as f64;
        (0..=nsteps)
            .into_par_iter()
            .map(|i| {
                let x = low + i as f64 * step;
                self.scan_point(par, x, values)
            })
            .collect()
    }

    fn setup_scan(
        &self,
        par: usize,
        nsteps: usize,
        low: f64,
        high: f64,
    ) -> (usize, f64, f64, Vec<f64>) {
        let nsteps = nsteps.clamp(2, 101);
        let p = self.params.trafo().parameter(par);
        let val = p.value();
        let err = p.error();

        let (low, high) = if (low - high).abs() < 1e-15 {
            let lo = val - 2.0 * err;
            let hi = val + 2.0 * err;
            self.clamp_scan_bounds(
                lo,
                hi,
                p.has_lower_limit(),
                p.lower_limit(),
                p.has_upper_limit(),
                p.upper_limit(),
            )
        } else {
            self.clamp_scan_bounds(
                low,
                high,
                p.has_lower_limit(),
                p.lower_limit(),
                p.has_upper_limit(),
                p.upper_limit(),
            )
        };

        // Build parameter vector at minimum
        let nparams = self.params.len();
        let values: Vec<f64> = (0..nparams)
            .map(|i| self.params.trafo().parameter(i).value())
            .collect();

        (nsteps, low, high, values)
    }

    /// Current best function value (may have been updated by scan).
    pub fn fval(&self) -> f64 {
        self.fval
    }

    /// Current parameters (may have been updated by scan).
    pub fn params(&self) -> &MnUserParameters {
        &self.params
    }
}

/// High-level scan builder working with a FunctionMinimum.
pub struct MnScan<'a, F: FCN + ?Sized> {
    fcn: &'a F,
    minimum: &'a FunctionMinimum,
}

impl<'a, F: FCN + ?Sized> MnScan<'a, F> {
    /// Create a new high-level scan from a minimization result.
    pub fn new(fcn: &'a F, minimum: &'a FunctionMinimum) -> Self {
        Self { fcn, minimum }
    }

    /// Scan parameter `par` over `nsteps` points.
    ///
    /// If `low == high == 0.0`, auto-range to +/- 2*error.
    pub fn scan(&self, par: usize, nsteps: usize, low: f64, high: f64) -> Vec<(f64, f64)> {
        self.scan_serial(par, nsteps, low, high)
    }

    /// Serial scan implementation.
    pub fn scan_serial(&self, par: usize, nsteps: usize, low: f64, high: f64) -> Vec<(f64, f64)> {
        let mut scanner =
            MnParameterScan::new(self.fcn, self.build_user_parameters(), self.minimum.fval());
        scanner.scan_serial(par, nsteps, low, high)
    }

    /// Parallel scan implementation (requires `parallel` feature).
    #[cfg(feature = "parallel")]
    pub fn scan_parallel(&self, par: usize, nsteps: usize, low: f64, high: f64) -> Vec<(f64, f64)>
    where
        F: Sync,
    {
        let mut scanner =
            MnParameterScan::new(self.fcn, self.build_user_parameters(), self.minimum.fval());
        scanner.scan_parallel(par, nsteps, low, high)
    }

    fn build_user_parameters(&self) -> MnUserParameters {
        // Build MnUserParameters from the minimum
        let user_state = self.minimum.user_state();
        let nparams = user_state.len();

        let mut params = MnUserParameters::new();
        for i in 0..nparams {
            add_param_from_state(&mut params, user_state.parameter(i));
        }

        params
    }
}

fn add_param_from_state(params: &mut MnUserParameters, p: &crate::parameter::MinuitParameter) {
    if p.has_limits() {
        params.add_limited(
            p.name(),
            p.value(),
            p.error(),
            p.lower_limit(),
            p.upper_limit(),
        );
    } else if p.has_lower_limit() {
        params.add_lower_limited(p.name(), p.value(), p.error(), p.lower_limit());
    } else if p.has_upper_limit() {
        params.add_upper_limited(p.name(), p.value(), p.error(), p.upper_limit());
    } else if p.is_const() {
        params.add_const(p.name(), p.value());
    } else {
        params.add(p.name(), p.value(), p.error());
    }
}

impl<F: FCN + ?Sized> MnParameterScan<'_, F> {
    fn clamp_scan_bounds(
        &self,
        low: f64,
        high: f64,
        has_lower: bool,
        lower: f64,
        has_upper: bool,
        upper: f64,
    ) -> (f64, f64) {
        let low = if has_lower { low.max(lower) } else { low };
        let high = if has_upper { high.min(upper) } else { high };
        (low, high)
    }
}
