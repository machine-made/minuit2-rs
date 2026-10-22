//! Hessian-based gradient refinement.
//!
//! Refines the gradient using second-derivative information from the Hessian
//! diagonal. The step sizes use `g2` from the Hessian computation.

use nalgebra::DVector;

use crate::minimum::gradient::FunctionGradient;
use crate::minimum::parameters::MinimumParameters;
use crate::mn_fcn::MnFcn;
use crate::strategy::MnStrategy;
use crate::user_transformation::MnUserTransformation;

pub struct HessianGradientCalculator;

impl HessianGradientCalculator {
    pub fn ncycle(strategy: &MnStrategy) -> u32 {
        strategy.hessian_gradient_ncycles()
    }

    /// Refine gradient using Hessian diagonal information.
    ///
    /// `hessian_g2` contains the second derivatives from the Hessian diagonal step.
    /// `hessian_gstep` contains the step sizes used during Hessian computation.
    pub fn compute(
        fcn: &MnFcn,
        params: &MinimumParameters,
        trafo: &MnUserTransformation,
        strategy: &MnStrategy,
        hessian_g2: &DVector<f64>,
        hessian_gstep: &DVector<f64>,
    ) -> FunctionGradient {
        let n = trafo.variable_parameters();
        let eps2 = trafo.precision().eps2();
        let fcnmin = params.fval();
        let dfmin = 8.0 * eps2 * (fcnmin.abs() + fcn.up());
        let vrysml = 8.0 * eps2 * eps2;

        let x = params.vec();
        let ncycles = Self::ncycle(strategy);
        let step_tol = strategy.grad_step_tol();
        let grad_tol = strategy.grad_tol();

        let mut grad = DVector::zeros(n);
        let mut g2 = DVector::zeros(n);
        let mut gstep = DVector::zeros(n);

        for i in 0..n {
            let ext_idx = trafo.ext_of_int(i);
            let xi = x[i];
            let p = &trafo.parameters()[ext_idx];
            let has_limits = p.has_limits() || p.has_lower_limit() || p.has_upper_limit();

            // Use Hessian g2 and gstep as starting point
            let mut g2i = hessian_g2[i];
            let mut gstepi = hessian_gstep[i].max(vrysml);

            for cycle in 0..ncycles {
                let optstp = (dfmin / (g2i.abs() + eps2)).sqrt();
                let mut step = optstp.max(0.1 * gstepi.abs());

                if has_limits {
                    step = step.min(0.5);
                }

                let stpmax = 10.0 * gstepi.abs();
                let stpmin = vrysml.max(8.0 * eps2 * xi.abs());
                step = if stpmax >= stpmin {
                    step.clamp(stpmin, stpmax)
                } else {
                    stpmin
                };

                let stepb4 = gstepi;
                let grdb4 = grad[i];

                gstepi = step;

                let mut xp = x.clone();
                let mut xm = x.clone();
                xp[i] = xi + step;
                xm[i] = xi - step;

                let fp = fcn.call(xp.as_slice());
                let fm = fcn.call(xm.as_slice());

                let grdi = 0.5 * (fp - fm) / step;
                let g2i_new = (fp + fm - 2.0 * fcnmin) / (step * step);

                grad[i] = grdi;
                g2[i] = g2i_new;
                gstep[i] = gstepi;
                g2i = g2i_new;

                if cycle > 0 {
                    let step_change = (gstepi - stepb4).abs() / gstepi.abs();
                    if step_change < step_tol {
                        break;
                    }
                    let grad_change = (grdi - grdb4).abs() / (grdi.abs() + dfmin / step);
                    if grad_change < grad_tol {
                        break;
                    }
                }
            }
        }

        FunctionGradient::new(grad, g2, gstep)
    }
}
