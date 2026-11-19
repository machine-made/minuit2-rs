//! Analytical gradient from the gradient of the function.
//!
//! Each internal derivative is the external derivative times dx/dint.
//! The seed stores a Hessian diagonal when a Hessian exists.
//! The seed stores a user g2 when no Hessian exists.
//! After a step, the escape prefers a finite user g2.
//! A function with no g2 and no Hessian keeps the step-size estimate.

use nalgebra::DVector;

use crate::fcn::{FCN, FCNGradient};
use crate::minimum::gradient::FunctionGradient;
use crate::minimum::parameters::MinimumParameters;
use crate::user_transformation::MnUserTransformation;

/// Convert user derivatives to internal coordinates using the parameter transformation.
/// # Examples
///
/// ```
/// use minuit2::{FCN, FCNGradient};
/// struct Square;
/// impl FCN for Square {
///     fn value(&self, p: &[f64]) -> f64 { p[0] * p[0] }
/// }
/// impl FCNGradient for Square {
///     fn gradient(&self, p: &[f64]) -> Vec<f64> { vec![2.0 * p[0]] }
/// }
/// use minuit2::{MinuitParameter, MnUserTransformation};
/// use minuit2::minimum::parameters::MinimumParameters;
/// use minuit2::gradient::AnalyticalGradientCalculator;
/// let t = MnUserTransformation::new(vec![MinuitParameter::new(0, "x", 3.0, 0.1)]);
/// let p = MinimumParameters::new(nalgebra::DVector::from_vec(vec![3.0]), 9.0);
/// let g = AnalyticalGradientCalculator::compute(&Square, &t, &p);
/// assert!((g.grad()[0] - 6.0).abs() < 1e-12);
/// ```
pub struct AnalyticalGradientCalculator;

impl AnalyticalGradientCalculator {
    /// Transform the gradient of the function into internal coordinates.
    ///
    /// The calculator multiplies each external derivative by dx/dint.
    /// A user g2 is stored only when the function has no Hessian.
    /// A short or non-finite g2 is not stored.
    /// The step-size estimate stays in that case.
    ///
    /// # Panics
    /// Panics if the parameter vector or returned external gradient is too short for the configured parameters.
    #[expect(
        clippy::manual_midpoint,
        reason = "keeps the operation order of the ROOT reference"
    )]
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

        // operator() stores a user g2 only when no Hessian exists.
        // A short vector or a non-finite entry is not a supplied g2.
        // The step-size estimate stays in that case.
        let supplied_g2 = if fcn.has_g2() && !fcn.has_hessian() {
            let ext_g2 = fcn.g2(&external_vals);
            if g2_is_supplied(&ext_g2, trafo, n) {
                Some(ext_g2)
            } else {
                None
            }
        } else {
            None
        };

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

            // dirin is the internal step from the parameter error.
            // gstep uses dirin. g2 uses it only when the function has no g2.
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

            // A stored g2 uses the factor dxdi squared.
            // The step-size estimate stays when no g2 is stored.
            let g2i = if let Some(ext_g2) = &supplied_g2 {
                dext_dint * dext_dint * ext_g2[ext_idx]
            } else {
                2.0 * error_def / (dirin * dirin)
            };
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

    /// Return whether the function advertises diagonal second derivatives.
    pub fn can_compute_g2(fcn: &dyn FCN) -> bool {
        fcn.has_g2() || fcn.has_hessian()
    }

    /// Return whether the function advertises a packed Hessian.
    pub fn can_compute_hessian(fcn: &dyn FCN) -> bool {
        fcn.has_hessian()
    }

    /// Return the external second derivative.
    ///
    /// A finite user g2 wins.
    /// The calculator skips a NaN or an infinity.
    /// The calculator reads the Hessian diagonal when no user g2 remains.
    pub fn g2(fcn: &dyn FCN, parameters: &[f64]) -> Option<Vec<f64>> {
        if fcn.has_g2() {
            let ext = fcn.g2(parameters);
            if !ext.is_empty() && ext.iter().all(|v| v.is_finite()) {
                return Some(ext);
            }
        }
        Self::hessian_diagonal(fcn, parameters)
    }

    /// Return the external Hessian diagonal.
    pub(crate) fn hessian_diagonal(fcn: &dyn FCN, parameters: &[f64]) -> Option<Vec<f64>> {
        if !fcn.has_hessian() {
            return None;
        }
        packed_diagonal(&fcn.hessian(parameters))
    }

    /// Return the external g2 that the escape stores after a step.
    ///
    /// A finite user g2 that covers every variable parameter wins.
    /// The calculator reads the Hessian diagonal when that g2 is absent.
    pub(crate) fn external_second_derivative(
        fcn: &dyn FCN,
        trafo: &MnUserTransformation,
        parameters: &[f64],
    ) -> Option<Vec<f64>> {
        let n = trafo.variable_parameters();
        if fcn.has_g2() {
            let ext = fcn.g2(parameters);
            if g2_is_supplied(&ext, trafo, n) {
                return Some(ext);
            }
        }
        Self::hessian_diagonal(fcn, parameters)
    }

    /// Return the packed Hessian of the function.
    pub fn hessian(fcn: &dyn FCN, parameters: &[f64]) -> Option<Vec<f64>> {
        if fcn.has_hessian() {
            Some(fcn.hessian(parameters))
        } else {
            None
        }
    }
}

/// Return true when every variable parameter has a finite g2 entry.
fn g2_is_supplied(ext_g2: &[f64], trafo: &MnUserTransformation, n_int: usize) -> bool {
    ext_g2.iter().all(|v| v.is_finite()) && (0..n_int).all(|i| trafo.ext_of_int(i) < ext_g2.len())
}

/// Return the diagonal of a packed upper triangle.
fn packed_diagonal(h: &[f64]) -> Option<Vec<f64>> {
    if h.is_empty() {
        return None;
    }
    let mut n = 0usize;
    while n * (n + 1) / 2 < h.len() {
        n += 1;
    }
    if n * (n + 1) / 2 != h.len() {
        return None;
    }
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let diag = i * (i + 3) / 2;
        out.push(h[diag]);
    }
    Some(out)
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
        assert_eq!(grad.grad().len(), 2);
        assert_eq!(grad.g2().len(), 2);
        assert_eq!(grad.gstep().len(), 2);
    }
}
