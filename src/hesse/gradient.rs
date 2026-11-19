//! Hessian-based gradient refinement.
//!
//! Start from the diagonal sweep's gradient, curvature and steps.
//! Follow ROOT's HessianGradientCalculator.cxx:68-132 for probe sizes and convergence.

use crate::minimum::gradient::FunctionGradient;
use crate::minimum::parameters::MinimumParameters;
use crate::mn_fcn::MnFcn;
use crate::strategy::MnStrategy;
use crate::user_transformation::MnUserTransformation;

/// Refine internal first derivatives using the Hessian sweep curvature and steps.
/// # Examples
///
/// ```
/// use minuit2::{MinuitParameter, MnStrategy, MnUserTransformation};
/// use minuit2::mn_fcn::MnFcn;
/// let trafo = MnUserTransformation::new(vec![MinuitParameter::new(0, "x", 0.0, 0.1)]);
/// let f = |p: &[f64]| (p[0] - 1.0).powi(2);
/// let fcn = MnFcn::new(&f, &trafo);
/// let strategy = MnStrategy::default();
/// let seed = minuit2::migrad::seed::MigradSeedGenerator::generate(&fcn, &trafo, &strategy);
/// use minuit2::hesse::gradient::HessianGradientCalculator;
/// let g = HessianGradientCalculator::compute(&fcn, seed.parameters(), &trafo, &strategy, seed.gradient());
/// assert!((g.grad()[0] + 2.0).abs() < 1e-3);
/// ```
pub struct HessianGradientCalculator;

impl HessianGradientCalculator {
    /// Return the maximum number of gradient refinement cycles from the strategy.
    #[must_use]
    pub fn ncycle(strategy: &MnStrategy) -> u32 {
        strategy.hessian_gradient_ncycles()
    }

    /// Refine first derivatives while keeping the diagonal sweep's curvature estimates.
    ///
    /// # Panics
    /// Panics if derivative vectors or the coordinate mapping do not cover the internal parameters.
    pub fn compute(
        fcn: &MnFcn,
        params: &MinimumParameters,
        trafo: &MnUserTransformation,
        strategy: &MnStrategy,
        gradient: &FunctionGradient,
    ) -> FunctionGradient {
        let precision = trafo.precision();
        let eps2 = precision.eps2();
        let dfmin = 4.0 * eps2 * (params.fval().abs() + fcn.up());
        let mut x = params.vec().clone();
        let mut grad = gradient.grad().clone();
        let mut gstep = gradient.gstep().clone();
        let g2 = gradient.g2();

        for i in 0..x.len() {
            let xi = x[i];
            let minimum_step = 4.0 * eps2 * (xi + eps2);
            let epspri = eps2 + (grad[i] * eps2).abs();
            let optstp = (dfmin / (g2[i].abs() + epspri)).sqrt();
            let mut d = (0.2 * gstep[i].abs()).min(optstp).max(minimum_step);
            let mut previous_change = 10000.0;

            for cycle in 0..Self::ncycle(strategy) {
                x[i] = xi + d;
                let fp = fcn.call(x.as_slice());
                x[i] = xi - d;
                let fm = fcn.call(x.as_slice());
                x[i] = xi;

                let old = grad[i];
                let new = (fp - fm) / (2.0 * d);
                let roundoff = precision.eps() * (fp.abs() + fm.abs()) / d;
                // ROOT stops on an exactly zero gradient; classify() avoids a float comparison
                // and keeps both signed zeros, with no epsilon.
                if matches!(new.classify(), std::num::FpCategory::Zero) {
                    break;
                }
                let change = ((old - new) / new).abs();
                if change > previous_change && cycle > 1 {
                    break;
                }
                previous_change = change;
                grad[i] = new;
                gstep[i] = d;

                if change < 0.05 || (old - new).abs() < roundoff || d < minimum_step {
                    break;
                }
                d *= 0.2;
            }
        }

        FunctionGradient::new(grad, g2.clone(), gstep)
    }
}
