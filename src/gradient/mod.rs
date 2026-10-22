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
pub trait GradientCalculator {
    fn compute(
        &self,
        fcn: &MnFcn,
        params: &MinimumParameters,
        trafo: &MnUserTransformation,
    ) -> FunctionGradient;
}
