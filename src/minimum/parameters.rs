//! Internal parameter values at a point in the minimization.
//!
//! Stores the parameter vector in internal space, step sizes, and function
//! value at this point.

use nalgebra::DVector;

/// Internal parameter values, objective value, and an optional step vector.
/// # Examples
///
/// ```
/// use minuit2::minimum::parameters::MinimumParameters;
/// use nalgebra::DVector;
/// let p = MinimumParameters::new(DVector::from_vec(vec![2.0]), 4.0);
/// assert!(p.is_valid());
/// assert!((p.fval() - 4.0).abs() < 1e-12);
/// ```
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
    #[must_use]
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
    #[must_use]
    pub fn with_step(vec: DVector<f64>, step: DVector<f64>, fval: f64) -> Self {
        Self {
            vec,
            step,
            fval,
            valid: true,
            has_step: true,
        }
    }

    /// Return the internal parameter vector.
    #[must_use]
    pub fn vec(&self) -> &DVector<f64> {
        &self.vec
    }

    /// Return the stored step vector in internal coordinates.
    #[must_use]
    pub fn step(&self) -> &DVector<f64> {
        &self.step
    }

    /// Return the objective value.
    #[must_use]
    pub fn fval(&self) -> f64 {
        self.fval
    }

    /// Return the stored point validity flag.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.valid
    }

    /// Return whether a step vector was supplied.
    #[must_use]
    pub fn has_step(&self) -> bool {
        self.has_step
    }

    /// Return whether a step vector was supplied.
    #[must_use]
    pub fn has_step_size(&self) -> bool {
        self.has_step()
    }
}
