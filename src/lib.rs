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
//! The crate contains no unsafe code, and the lint denies it in every module, the `python` feature included.
//! No features are enabled by default.
//! The `python` feature provides Python bindings, and `trace` enables iteration tracing.
//!
//! Not proven: bit-identity to ROOT is not claimed.
//! Measured distances in `parity/parity_ulp.json` cover 15 workloads and 218 fields, with 4 failed fields.
//! The maximum ULP distance is 301375601036874150.
//! Of the failed fields, 2 are fixed-parameter call counts with manifest waivers; 2 are signed zeros in covariance entries.
//! The large ULP distance comes from EDM in the two fixed-parameter workloads.
//! ROOT reports about 1e-25 and Rust about 1.5e-5, within the manifest absolute tolerance of 5e-05.
//! Assertions and unchecked indexing still panic on invalid caller input, as do some `unwrap` and `expect` calls.
//! A typed error interface remains open, so the crate does not claim production kernel status.
//!
//! # Quick Start — Migrad (recommended)
//!
//! ```
//! use minuit2::MnMigrad;
//!
//! let result = MnMigrad::new()
//!     .add("x", 0.0, 0.1)
//!     .add("y", 0.0, 0.1)
//!     // Tighten the objective stopping tolerance so both parameter errors stay below 1e-3.
//!     .tolerance(1e-8)
//!     .minimize(&|p: &[f64]| {
//!         (1.0 - p[0]).powi(2) + 100.0 * (p[1] - p[0] * p[0]).powi(2)
//!     });
//!
//! assert!(result.is_valid());
//! assert!(result.params().iter().all(|x| (x - 1.0).abs() < 1e-3));
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
//!     // Tighten the simplex spread tolerance so both parameter errors stay below 1e-3.
//!     .tolerance(1e-8)
//!     .minimize(&|p: &[f64]| {
//!         (1.0 - p[0]).powi(2) + 100.0 * (p[1] - p[0] * p[0]).powi(2)
//!     });
//!
//! assert!(result.is_valid());
//! assert!(result.params().iter().all(|x| (x - 1.0).abs() < 1e-3));
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
/// Objective function and derivative contracts.
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
/// Individual external parameter definitions.
pub mod parameter;
pub mod posdef;
/// Machine precision used for numerical steps.
pub mod precision;
pub mod print;
pub mod scan;
pub mod simplex;
/// Derivative refinement strategy settings.
pub mod strategy;
pub mod transform;
/// Packed symmetric covariance storage in external coordinates.
pub mod user_covariance;
pub mod user_parameter_state;
pub mod user_parameters;
pub mod user_transformation;

/// Python bindings for fitting, error analysis, profiles, and contours.
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
