//! Function gradient in internal parameter space.
//!
//! Stores the gradient vector plus auxiliary vectors `g2` (second derivative
//! estimates) and `gstep` (step sizes used for numerical differentiation).

use nalgebra::DVector;

#[derive(Debug, Clone)]
pub struct FunctionGradient {
    /// First derivatives ∂f/∂p_i.
    grad: DVector<f64>,
    /// Second derivative estimates ∂²f/∂p_i².
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

    pub fn grad(&self) -> &DVector<f64> {
        &self.grad
    }

    pub fn g2(&self) -> &DVector<f64> {
        &self.g2
    }

    pub fn gstep(&self) -> &DVector<f64> {
        &self.gstep
    }

    pub fn is_valid(&self) -> bool {
        self.valid
    }

    pub fn is_analytical(&self) -> bool {
        self.analytical
    }

    pub fn set_valid(&mut self, valid: bool) {
        self.valid = valid;
    }

    pub fn set_analytical(&mut self, analytical: bool) {
        self.analytical = analytical;
    }
}
