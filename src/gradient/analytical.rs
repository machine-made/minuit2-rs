//! Analytical gradient calculator from user-provided gradients.
//!
//! Takes gradients computed by the user (in external parameter space) and
//! transforms them to internal parameter space using the chain rule:
//!   `g_int[i] = g_ext[i] * dext/dint[i]`
//!
//! Also computes g2 (second derivative heuristic) and gstep (step sizes) using
//! the same logic as InitialGradientCalculator.

use nalgebra::DVector;

use crate::fcn::{FCN, FCNGradient};
use crate::minimum::gradient::FunctionGradient;
use crate::minimum::parameters::MinimumParameters;
use crate::user_transformation::MnUserTransformation;

pub struct AnalyticalGradientCalculator;

impl AnalyticalGradientCalculator {
    /// Compute gradient from user-provided analytical gradient.
    ///
    /// Takes the user's gradient (in external parameter space), transforms it
    /// to internal space using the chain rule, and provides g2 and gstep heuristics.
    pub fn compute(
        fcn: &dyn FCNGradient,
        trafo: &MnUserTransformation,
        params: &MinimumParameters,
    ) -> FunctionGradient {
        let n = trafo.variable_parameters();
        let eps2 = trafo.precision().eps2();

        // Get the external parameter values
        let internal_vec = params.vec();
        let external_vals = trafo.transform(internal_vec.as_slice());

        // Call user's analytical gradient (returns gradient in external space)
        let ext_gradient = fcn.gradient(&external_vals);

        // Transform gradient from external to internal space using chain rule:
        // g_int[i] = g_ext[ext_idx] * dext/dint[i]
        let mut grad = DVector::zeros(n);
        let mut g2 = DVector::zeros(n);
        let mut gstep = DVector::zeros(n);

        let error_def = fcn.error_def();

        for i in 0..n {
            let ext_idx = trafo.ext_of_int(i);
            let int_val = internal_vec[i];

            // Get the derivative of the transform
            let dext_dint = trafo.dint2ext(ext_idx, int_val);

            // Transform gradient: chain rule
            let g_ext = ext_gradient[ext_idx];
            let g_int = g_ext * dext_dint;

            // Compute g2 heuristic (second derivative estimate)
            // Use the same heuristic as InitialGradientCalculator:
            // g2 = 2 * error_def / dirin^2
            let werr = trafo.parameters()[ext_idx].error();
            let sav = trafo.int2ext(ext_idx, int_val);

            // Forward step: external + werr, clamped to upper limit
            let mut sav_plus = sav + werr;
            let p = &trafo.parameters()[ext_idx];
            if p.has_upper_limit() && sav_plus > p.upper_limit() {
                sav_plus = p.upper_limit();
            }
            let var_plus = trafo.ext2int(ext_idx, sav_plus);
            let vplu = var_plus - int_val;

            // Backward step: external - werr, clamped to lower limit
            let mut sav_minus = sav - werr;
            if p.has_lower_limit() && sav_minus < p.lower_limit() {
                sav_minus = p.lower_limit();
            }
            let var_minus = trafo.ext2int(ext_idx, sav_minus);
            let vmin = var_minus - int_val;

            // Minimum step size from machine precision
            let gsmin = 8.0 * eps2 * (int_val.abs() + eps2);

            // Direction magnitude: average of forward/backward internal steps
            let dirin = (0.5 * (vplu.abs() + vmin.abs())).max(gsmin);

            // Heuristic g2: assumes parabolic shape
            let g2i = 2.0 * error_def / (dirin * dirin);
            let mut gstepi = gsmin.max(0.1 * dirin);

            // For limited parameters, cap step at 0.5
            if p.has_limits() && gstepi > 0.5 {
                gstepi = 0.5;
            }

            grad[i] = g_int;
            g2[i] = g2i;
            gstep[i] = gstepi;
        }

        // Return analytical gradient with computed g2 and gstep
        let mut result = FunctionGradient::new(grad, g2, gstep);
        result.set_analytical(true);
        result
    }

    pub fn can_compute_g2(fcn: &dyn FCN) -> bool {
        fcn.has_g2() || fcn.has_hessian()
    }

    pub fn can_compute_hessian(fcn: &dyn FCN) -> bool {
        fcn.has_hessian()
    }

    pub fn g2(fcn: &dyn FCN, parameters: &[f64]) -> Option<Vec<f64>> {
        if fcn.has_g2() {
            Some(fcn.g2(parameters))
        } else if fcn.has_hessian() {
            let h = fcn.hessian(parameters);
            if h.is_empty() {
                None
            } else {
                let mut n = 0usize;
                while n * (n + 1) / 2 < h.len() {
                    n += 1;
                }
                if n * (n + 1) / 2 != h.len() {
                    return None;
                }
                let mut out = Vec::new();
                for i in 0..n {
                    let diag = i * (i + 3) / 2;
                    out.push(h[diag]);
                }
                Some(out)
            }
        } else {
            None
        }
    }

    pub fn hessian(fcn: &dyn FCN, parameters: &[f64]) -> Option<Vec<f64>> {
        if fcn.has_hessian() {
            Some(fcn.hessian(parameters))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fcn::{FCN, FCNGradient};
    use crate::parameter::MinuitParameter;

    /// Quadratic with analytical gradient
    struct QuadraticWithGrad;

    impl FCN for QuadraticWithGrad {
        fn value(&self, p: &[f64]) -> f64 {
            // f(x,y) = x² + 4y²
            p[0] * p[0] + 4.0 * p[1] * p[1]
        }
    }

    impl FCNGradient for QuadraticWithGrad {
        fn gradient(&self, p: &[f64]) -> Vec<f64> {
            // df/dx = 2x, df/dy = 8y
            vec![2.0 * p[0], 8.0 * p[1]]
        }
    }

    #[test]
    fn analytical_gradient_quadratic() {
        let params = vec![
            MinuitParameter::new(0, "x", 3.0, 0.1),
            MinuitParameter::new(1, "y", 2.0, 0.1),
        ];
        let trafo = MnUserTransformation::new(params);

        // Evaluate at (3, 2)
        let x = DVector::from_vec(vec![3.0, 2.0]);
        let f_val = 9.0 + 16.0; // 25
        let min_params = MinimumParameters::new(x, f_val);

        let grad = AnalyticalGradientCalculator::compute(&QuadraticWithGrad, &trafo, &min_params);

        // df/dx = 6, df/dy = 16
        assert!(
            (grad.grad()[0] - 6.0).abs() < 1e-10,
            "dfdx should be ~6.0, got {}",
            grad.grad()[0]
        );
        assert!(
            (grad.grad()[1] - 16.0).abs() < 1e-10,
            "dfdy should be ~16.0, got {}",
            grad.grad()[1]
        );

        // Check that g2 and gstep are computed (non-zero)
        assert!(grad.g2()[0] > 0.0, "g2[0] should be positive");
        assert!(grad.g2()[1] > 0.0, "g2[1] should be positive");
        assert!(grad.gstep()[0] > 0.0, "gstep[0] should be positive");
        assert!(grad.gstep()[1] > 0.0, "gstep[1] should be positive");

        // Check that it's marked as analytical
        assert!(
            grad.is_analytical(),
            "gradient should be marked as analytical"
        );
    }

    #[test]
    fn analytical_gradient_bounded_param() {
        // Test with bounded parameter: [0, 10]
        let params = vec![
            MinuitParameter::with_limits(0, "x", 5.0, 0.1, 0.0, 10.0),
            MinuitParameter::new(1, "y", 2.0, 0.1),
        ];
        let trafo = MnUserTransformation::new(params);

        // Evaluate at (5, 2)
        let x = DVector::from_vec(vec![0.0, 2.0]); // internal values
        let min_params = MinimumParameters::new(x, 100.0);

        let grad = AnalyticalGradientCalculator::compute(&QuadraticWithGrad, &trafo, &min_params);

        // Should not panic and should produce valid gradients
        assert!(grad.is_valid());
        assert!(grad.is_analytical());
        assert!(grad.grad().len() == 2);
        assert!(grad.g2().len() == 2);
        assert!(grad.gstep().len() == 2);
    }
}
