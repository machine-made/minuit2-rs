//! Display implementations for minimization results.
//!
//! Uses Rust's `Display` trait.

use std::fmt;

use crate::minimum::FunctionMinimum;

impl fmt::Display for FunctionMinimum {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "FunctionMinimum:")?;
        writeln!(f, "  valid:     {}", self.is_valid())?;
        writeln!(f, "  fval:      {:.6e}", self.fval())?;
        writeln!(f, "  edm:       {:.6e}", self.edm())?;
        writeln!(f, "  nfcn:      {}", self.nfcn())?;

        if self.is_above_max_edm() {
            writeln!(f, "  WARNING: EDM above maximum")?;
        }
        if self.reached_call_limit() {
            writeln!(f, "  WARNING: call limit reached")?;
        }

        writeln!(f, "  parameters:")?;
        let state = self.user_state();
        for i in 0..state.len() {
            let p = state.parameter(i);
            write!(
                f,
                "    {:>4}  {:>12}  {:>14.6e}  +/- {:>10.6e}",
                i,
                p.name(),
                p.value(),
                p.error(),
            )?;
            if p.is_fixed() {
                write!(f, "  (fixed)")?;
            }
            if p.has_limits() {
                write!(f, "  [{:.4e}, {:.4e}]", p.lower_limit(), p.upper_limit())?;
            } else if p.has_lower_limit() {
                write!(f, "  [{:.4e}, +inf)", p.lower_limit())?;
            } else if p.has_upper_limit() {
                write!(f, "  (-inf, {:.4e}]", p.upper_limit())?;
            }
            writeln!(f)?;
        }

        Ok(())
    }
}
