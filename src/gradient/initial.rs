//! Initial gradient calculator from parameter step sizes.
//!
//! Produces a rough heuristic gradient estimate from the user-provided
//! parameter errors and the error definition (`Up` value). It does not evaluate
//! the FCN.
//!
//! The key formula: g2 = 2 * error_def / dirin^2, grad = g2 * dirin,
//! where dirin is the internal-space step derived from the external error.

use nalgebra::DVector;

use crate::minimum::gradient::FunctionGradient;
use crate::minimum::parameters::MinimumParameters;
use crate::mn_fcn::MnFcn;
use crate::strategy::MnStrategy;
use crate::user_transformation::MnUserTransformation;

pub struct InitialGradientCalculator {
    _strategy: MnStrategy,
}

impl InitialGradientCalculator {
    pub fn new(strategy: MnStrategy) -> Self {
        Self {
            _strategy: strategy,
        }
    }

    pub fn compute(
        &self,
        fcn: &MnFcn,
        params: &MinimumParameters,
        trafo: &MnUserTransformation,
    ) -> FunctionGradient {
        let n = trafo.variable_parameters();
        let mut grad = DVector::zeros(n);
        let mut g2 = DVector::zeros(n);
        let mut gstep = DVector::zeros(n);

        let eps2 = trafo.precision().eps2();
        let error_def = fcn.error_def();

        for i in 0..n {
            let ext_idx = trafo.ext_of_int(i);
            let var = params.vec()[i]; // internal value
            let werr = trafo.parameters()[ext_idx].error();
            let p = &trafo.parameters()[ext_idx];

            // Convert internal → external, add werr, convert back
            let sav = trafo.int2ext(ext_idx, var);

            // Forward step: external + werr, clamped to upper limit
            let mut sav_plus = sav + werr;
            if p.has_upper_limit() && sav_plus > p.upper_limit() {
                sav_plus = p.upper_limit();
            }
            let var_plus = trafo.ext2int(ext_idx, sav_plus);
            let vplu = var_plus - var;

            // Backward step: external - werr, clamped to lower limit
            let mut sav_minus = sav - werr;
            if p.has_lower_limit() && sav_minus < p.lower_limit() {
                sav_minus = p.lower_limit();
            }
            let var_minus = trafo.ext2int(ext_idx, sav_minus);
            let vmin = var_minus - var;

            // Minimum step size from machine precision
            let gsmin = 8.0 * eps2 * (var.abs() + eps2);

            // Direction magnitude: average of forward/backward internal steps
            let dirin = (0.5 * (vplu.abs() + vmin.abs())).max(gsmin);

            // Heuristic gradient: assumes parabolic shape with curvature from error_def
            let g2i = 2.0 * error_def / (dirin * dirin);
            let grdi = g2i * dirin;
            let mut gstepi = (gsmin).max(0.1 * dirin);

            // For limited parameters, cap step at 0.5
            if p.has_limits() && gstepi > 0.5 {
                gstepi = 0.5;
            }

            grad[i] = grdi;
            g2[i] = g2i;
            gstep[i] = gstepi;
        }

        FunctionGradient::new(grad, g2, gstep)
    }
}
