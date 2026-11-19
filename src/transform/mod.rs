//! Parameter transformations between external (user) and internal (optimizer) spaces.
//!
//! Bounded parameters are transformed to an unbounded internal space so the
//! optimizer can search freely. Three transforms cover all bound combinations:
//! - `SinTransform`: both upper and lower bounds
//! - `SqrtLowTransform`: lower bound only
//! - `SqrtUpTransform`: upper bound only

/// Coordinate transformation for parameters with two bounds.
pub mod sin;
/// Coordinate transformation for a lower bound.
pub mod sqrt_low;
/// Coordinate transformation for an upper bound.
pub mod sqrt_up;

pub use sin::SinTransform;
pub use sqrt_low::SqrtLowTransform;
pub use sqrt_up::SqrtUpTransform;

/// Common interface for parameter transformations.
/// Interface for parameter transformations (bounded <-> unbounded).
/// # Examples
///
/// ```
/// use minuit2::transform::{SinTransform, ParameterTransform};
/// use minuit2::MnMachinePrecision;
/// let t = SinTransform;
/// let internal = t.ext2int(2.0, 4.0, 0.0, &MnMachinePrecision::new());
/// assert!((t.int2ext(internal, 4.0, 0.0) - 2.0).abs() < 1e-12);
/// ```
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
