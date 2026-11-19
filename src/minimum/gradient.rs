//! Function gradient in internal parameter space.
//!
//! Stores the gradient vector plus auxiliary vectors `g2` (second derivative
//! estimates) and `gstep` (step sizes used for numerical differentiation).

use nalgebra::DVector;

/// Internal first derivatives, diagonal second derivatives, and numerical differentiation steps.
/// # Examples
///
/// ```
/// use minuit2::minimum::gradient::FunctionGradient;
/// use nalgebra::DVector;
/// let g = FunctionGradient::analytical(DVector::from_vec(vec![4.0]));
/// assert!(g.is_analytical());
/// assert!((g.grad()[0] - 4.0).abs() < 1e-12);
/// ```
#[derive(Debug, Clone)]
pub struct FunctionGradient {
    /// First derivatives ∂`f/∂p_i`.
    grad: DVector<f64>,
    /// Second derivative estimates ∂`²f/∂p_i²`.
    g2: DVector<f64>,
    /// Step sizes used to compute the gradient.
    gstep: DVector<f64>,
    /// Whether this is a valid, converged gradient.
    valid: bool,
    /// Whether this is an analytical (user-provided) gradient.
    analytical: bool,
}

impl FunctionGradient {
    /// Numerical gradient with grad, g2, and gstep.
    #[must_use]
    pub fn new(grad: DVector<f64>, g2: DVector<f64>, gstep: DVector<f64>) -> Self {
        Self {
            grad,
            g2,
            gstep,
            valid: true,
            analytical: false,
        }
    }

    /// Analytical gradient (only grad vector, no g2/gstep).
    #[must_use]
    pub fn analytical(grad: DVector<f64>) -> Self {
        let n = grad.len();
        Self {
            grad,
            g2: DVector::zeros(n),
            gstep: DVector::zeros(n),
            valid: true,
            analytical: true,
        }
    }

    /// Return first derivatives with respect to internal parameters.
    #[must_use]
    pub fn grad(&self) -> &DVector<f64> {
        &self.grad
    }

    /// Return diagonal second derivatives with respect to internal parameters.
    #[must_use]
    pub fn g2(&self) -> &DVector<f64> {
        &self.g2
    }

    /// Return numerical differentiation steps in internal coordinates.
    #[must_use]
    pub fn gstep(&self) -> &DVector<f64> {
        &self.gstep
    }

    /// Return the stored gradient validity flag.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.valid
    }

    /// Return whether the gradient came from analytical derivatives.
    #[must_use]
    pub fn is_analytical(&self) -> bool {
        self.analytical
    }

    /// Set the stored validity flag.
    pub fn set_valid(&mut self, valid: bool) {
        self.valid = valid;
    }

    /// Set whether the gradient came from analytical derivatives.
    pub fn set_analytical(&mut self, analytical: bool) {
        self.analytical = analytical;
    }
}
