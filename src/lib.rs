//! minuit2 — Pure Rust implementation of Minuit-style parameter optimization.
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
