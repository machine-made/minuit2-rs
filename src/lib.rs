//! minuit2 — Pure Rust implementation of Minuit-style parameter optimization.
//!
//! # Contract
//!
//! This crate minimizes user functions and estimates parameter uncertainties with algorithms ported from ROOT's `Minuit2`.
//! Parameter vectors and covariance dimensions must match the configured parameters.
//! Fixed parameters keep their values, and parameter limits use internal coordinate transformations.
//! Rust minimizers report failure through validity flags on the minimum and covariance, rather than `Result`.
//! Callers must check these flags before using fitted values or uncertainties.
//!
//! With a deterministic function, the same inputs and build give bit-identical results on one CPU type.
//! The `parallel` feature evaluates scan points concurrently and returns them in input order.
//! Only the `python` feature permits unsafe code from pyo3 macro expansions.
//! No features are enabled by default.
//! The `python` feature provides Python bindings, and `trace` enables iteration tracing.
//!
//! Not proven: parity outside the 12 workloads of `reports/verification/diff_summary.md`
//!
//! # Quick Start — Migrad (recommended)
//!
//! ```
//! use minuit2::MnMigrad;
//!
//! let result = MnMigrad::new()
//!     .add("x", 0.0, 0.1)
//!     .add("y", 0.0, 0.1)
//!     .minimize(&|p: &[f64]| {
//!         (1.0 - p[0]).powi(2) + 100.0 * (p[1] - p[0] * p[0]).powi(2)
//!     });
//!
//! println!("{result}");
//! ```
//!
//! # Simplex (derivative-free)
//!
//! ```
//! use minuit2::MnSimplex;
//!
//! let result = MnSimplex::new()
//!     .add("x", 0.0, 0.1)
//!     .add("y", 0.0, 0.1)
//!     .minimize(&|p: &[f64]| {
//!         (1.0 - p[0]).powi(2) + 100.0 * (p[1] - p[0] * p[0]).powi(2)
//!     });
//!
//! println!("{result}");
//! ```

/// Repository used as the numerical reference baseline for parity work.
pub const REFERENCE_MINUIT2_REPOSITORY: &str = "https://github.com/root-project/root";
/// Minuit2 subtree used for reference-baseline metadata.
pub const REFERENCE_MINUIT2_SUBTREE: &str = "math/minuit2";
/// Target release baseline for functionality parity work.
pub const REFERENCE_MINUIT2_TAG: &str = "v6-36-08";
/// Pinned commit for the baseline release tag.
pub const REFERENCE_MINUIT2_COMMIT: &str = "a8ca1b23e38d7dbe0ff24027894ca0f2ad65f1bd";

pub mod application;
pub mod contours;
pub mod covariance_squeeze;
pub mod fcn;
pub mod global_cc;
pub mod gradient;
pub mod hesse;
pub mod linesearch;
pub mod migrad;
pub mod minimize;
pub mod minimum;
pub mod minos;
pub mod mn_fcn;
pub mod parabola;
pub mod parameter;
pub mod posdef;
pub mod precision;
pub mod print;
pub mod scan;
pub mod simplex;
pub mod strategy;
pub mod transform;
pub mod user_covariance;
pub mod user_parameter_state;
pub mod user_parameters;
pub mod user_transformation;

#[cfg(feature = "python")]
pub mod python;

// Re-exports for convenience
pub use contours::MnContours;
pub use fcn::{FCN, FCNGradient};
pub use hesse::MnHesse;
pub use migrad::MnMigrad;
pub use minimize::MnMinimize;
pub use minimum::FunctionMinimum;
pub use minos::MnMinos;
pub use parameter::MinuitParameter;
pub use precision::MnMachinePrecision;
pub use scan::{MnScan, MnScanMinimizer};
pub use simplex::MnSimplex;
pub use strategy::MnStrategy;
pub use user_covariance::MnUserCovariance;
pub use user_parameter_state::MnUserParameterState;
pub use user_parameters::MnUserParameters;
pub use user_transformation::MnUserTransformation;
