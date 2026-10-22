//! Seed generator for the SCAn minimizer.
//!
//! ROOT Minuit2's `ScanMinimizer` (v6-36-08
//! `math/minuit2/inc/Minuit2/ScanMinimizer.h`) composes the regular
//! `SimplexSeedGenerator` with `ScanBuilder`.  Keep the same composition here.

pub use crate::simplex::seed::SimplexSeedGenerator as ScanSeedGenerator;
