//! Gradient calculators.
//!
//! The `GradientCalculator` trait defines the interface. Concrete impls:
//! - `InitialGradientCalculator`: computes a first gradient estimate from step sizes
//! - `Numerical2PGradientCalculator`: two-point central differences
//! - `AnalyticalGradientCalculator`: user-provided analytical gradients

pub mod analytical;
pub mod initial;
pub mod numerical;

pub use analytical::AnalyticalGradientCalculator;
pub use initial::InitialGradientCalculator;
pub use numerical::Numerical2PGradientCalculator;

use crate::minimum::gradient::FunctionGradient;
use crate::minimum::parameters::MinimumParameters;
use crate::mn_fcn::MnFcn;
use crate::user_transformation::MnUserTransformation;

/// Trait for computing function gradients.
/// # Examples
///
/// ```
/// use minuit2::gradient::GradientCalculator;
/// use minuit2::minimum::{gradient::FunctionGradient, parameters::MinimumParameters};
/// use minuit2::{MnUserTransformation, mn_fcn::MnFcn};
/// struct ZeroGradient;
/// impl GradientCalculator for ZeroGradient {
///     fn compute(&self, _: &MnFcn, p: &MinimumParameters, _: &MnUserTransformation) -> FunctionGradient {
///         FunctionGradient::analytical(nalgebra::DVector::zeros(p.vec().len()))
///     }
/// }
/// let t = MnUserTransformation::new(vec![]);
/// let f = |_: &[f64]| 0.0;
/// let fcn = MnFcn::new(&f, &t);
/// let p = MinimumParameters::new(nalgebra::DVector::zeros(0), 0.0);
/// assert!(ZeroGradient.compute(&fcn, &p, &t).grad().is_empty());
/// ```
pub trait GradientCalculator {
    /// Compute internal derivatives at the supplied internal parameter point.
    fn compute(
        &self,
        fcn: &MnFcn,
        params: &MinimumParameters,
        trafo: &MnUserTransformation,
    ) -> FunctionGradient;
}
