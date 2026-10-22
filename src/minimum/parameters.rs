//! Internal parameter values at a point in the minimization.
//!
//! Stores the parameter vector in internal space, step sizes, and function
//! value at this point.

use nalgebra::DVector;

#[derive(Debug, Clone)]
pub struct MinimumParameters {
    /// Parameter values in internal space.
    vec: DVector<f64>,
    /// Step sizes used to reach this point.
    step: DVector<f64>,
    /// Function value at this point.
    fval: f64,
    /// Whether this is a valid point.
    valid: bool,
    /// Whether step has a direction (vs. just magnitudes).
    has_step: bool,
}

impl MinimumParameters {
    /// Create with just parameter values and function value.
    pub fn new(vec: DVector<f64>, fval: f64) -> Self {
        let n = vec.len();
        Self {
            vec,
            step: DVector::zeros(n),
            fval,
            valid: true,
            has_step: false,
        }
    }

    /// Create with parameter values, step direction, and function value.
    pub fn with_step(vec: DVector<f64>, step: DVector<f64>, fval: f64) -> Self {
        Self {
            vec,
            step,
            fval,
            valid: true,
            has_step: true,
        }
    }

    pub fn vec(&self) -> &DVector<f64> {
        &self.vec
    }

    pub fn step(&self) -> &DVector<f64> {
        &self.step
    }

    pub fn fval(&self) -> f64 {
        self.fval
    }

    pub fn is_valid(&self) -> bool {
        self.valid
    }

    pub fn has_step(&self) -> bool {
        self.has_step
    }

    pub fn has_step_size(&self) -> bool {
        self.has_step()
    }
}
