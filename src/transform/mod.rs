//! Parameter transformations between external (user) and internal (optimizer) spaces.
//!
//! Bounded parameters are transformed to an unbounded internal space so the
//! optimizer can search freely. Three transforms cover all bound combinations:
//! - `SinTransform`: both upper and lower bounds
//! - `SqrtLowTransform`: lower bound only
//! - `SqrtUpTransform`: upper bound only

pub mod sin;
pub mod sqrt_low;
pub mod sqrt_up;

pub use sin::SinTransform;
pub use sqrt_low::SqrtLowTransform;
pub use sqrt_up::SqrtUpTransform;

/// Common interface for parameter transformations.
/// Interface for parameter transformations (bounded <-> unbounded).
pub trait ParameterTransform {
    /// Transform internal value to external (user) space.
    fn int2ext(&self, value: f64, upper: f64, lower: f64) -> f64;

    /// Transform external value to internal (optimizer) space.
    fn ext2int(
        &self,
        value: f64,
        upper: f64,
        lower: f64,
        precision: &crate::precision::MnMachinePrecision,
    ) -> f64;

    /// Derivative d(external)/d(internal).
    fn dint2ext(&self, value: f64, upper: f64, lower: f64) -> f64;
}
