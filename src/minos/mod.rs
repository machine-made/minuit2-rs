//! MnMinos: asymmetric profile-likelihood errors.
//!
//! After Migrad, and optionally Hesse, MnMinos computes asymmetric errors by
//! finding where `F(x) = Fmin + Up` for each parameter.

pub mod cross;
pub mod function_cross;
pub mod minos_error;

pub use cross::MnCross;
pub use minos_error::MinosError;

use crate::application::default_max_fcn;
use crate::fcn::FCN;
use crate::minimum::FunctionMinimum;
use crate::strategy::MnStrategy;

/// Compute MINOS asymmetric errors.
pub struct MnMinos<'a> {
    fcn: &'a dyn FCN,
    minimum: &'a FunctionMinimum,
    strategy: MnStrategy,
    max_calls: Option<usize>,
    tolerance: f64,
}

impl<'a> MnMinos<'a> {
    /// Create a new MINOS error calculator.
    pub fn new(fcn: &'a dyn FCN, minimum: &'a FunctionMinimum) -> Self {
        Self {
            fcn,
            minimum,
            strategy: MnStrategy::default(),
            max_calls: None,
            tolerance: 0.1,
        }
    }

    /// Set strategy level.
    pub fn with_strategy(mut self, level: u32) -> Self {
        self.strategy = MnStrategy::new(level);
        self
    }

    /// Set maximum function calls.
    pub fn with_max_calls(mut self, max: usize) -> Self {
        self.max_calls = Some(max);
        self
    }

    /// Set tolerance for crossing convergence (default 0.1).
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self
    }

    /// Compute both upper and lower MINOS errors for parameter `par`.
    pub fn errors(&self, par: usize) -> (f64, f64) {
        let me = self.minos_error(par);
        (me.lower_error(), me.upper_error())
    }

    /// Alias for `errors`.
    pub fn minos(&self, par: usize) -> MinosError {
        self.minos_error(par)
    }

    /// Full MinosError (both crossings) for parameter `par`.
    pub fn minos_error(&self, par: usize) -> MinosError {
        let p = self.minimum.user_state().parameter(par);
        let min_val = p.value();
        let hesse_err = p.error();
        let lo = self.lower(par);
        let up = self.upper(par);
        MinosError::new(par, min_val, hesse_err, lo, up)
    }

    /// Lower crossing only.
    pub fn lower(&self, par: usize) -> MnCross {
        self.find_crossing(par, -1.0)
    }

    /// Upper crossing only.
    pub fn upper(&self, par: usize) -> MnCross {
        self.find_crossing(par, 1.0)
    }

    /// Alias for the lower crossing object.
    pub fn loval(&self, par: usize) -> MnCross {
        self.lower(par)
    }

    /// Alias for the upper crossing object.
    pub fn upval(&self, par: usize) -> MnCross {
        self.upper(par)
    }

    /// Helper with explicit crossing direction.
    pub fn find_cross_value(&self, dir: i32, par: usize, maxcalls: usize, toler: f64) -> MnCross {
        let direction = if dir < 0 { -1.0 } else { 1.0 };
        let nvar = self.minimum.n_variable_params();
        let default_calls = default_cross_calls(nvar);
        let maxcalls = if maxcalls == 0 {
            default_calls
        } else {
            maxcalls
        };

        let user_state = self.minimum.user_state();
        let p = user_state.parameter(par);
        let err = p.error();
        let val = p.value();
        let pdir = direction * err;
        let pmid = val + pdir;
        function_cross::find_crossing(
            self.fcn,
            self.minimum,
            par,
            pmid,
            pdir,
            toler,
            maxcalls,
            &self.strategy,
        )
    }

    fn find_crossing(&self, par: usize, direction: f64) -> MnCross {
        let nvar = self.minimum.n_variable_params();
        let maxcalls = self.max_calls.unwrap_or_else(|| default_cross_calls(nvar));

        let user_state = self.minimum.user_state();
        let p = user_state.parameter(par);
        let err = p.error();
        let val = p.value();

        // Check if parameter is fixed
        if p.is_fixed() || p.is_const() {
            return MnCross::invalid(0);
        }

        // The scan direction: parameter error scaled by direction
        let pdir = direction * err;

        // Starting point: current value + step in direction
        let pmid = val + pdir;

        // Check limits
        if direction > 0.0 && p.has_upper_limit() && pmid > p.upper_limit() {
            let pmid = p.upper_limit() - 1e-6 * (p.upper_limit() - val).abs().max(1e-10);
            return function_cross::find_crossing(
                self.fcn,
                self.minimum,
                par,
                pmid,
                pdir,
                self.tolerance,
                maxcalls,
                &self.strategy,
            );
        }

        if direction < 0.0 && p.has_lower_limit() && pmid < p.lower_limit() {
            let pmid = p.lower_limit() + 1e-6 * (val - p.lower_limit()).abs().max(1e-10);
            return function_cross::find_crossing(
                self.fcn,
                self.minimum,
                par,
                pmid,
                pdir,
                self.tolerance,
                maxcalls,
                &self.strategy,
            );
        }

        function_cross::find_crossing(
            self.fcn,
            self.minimum,
            par,
            pmid,
            pdir,
            self.tolerance,
            maxcalls,
            &self.strategy,
        )
    }
}

fn default_cross_calls(nvar: usize) -> usize {
    2 * (nvar + 1) * default_max_fcn(nvar)
}
