//! Call-counting FCN wrapper that operates in internal parameter space.
//!
//! Takes internal parameter vectors, transforms them to external space via
//! `MnUserTransformation`, calls the user's FCN, and counts every call.

use std::cell::Cell;

use crate::fcn::FCN;
use crate::user_transformation::MnUserTransformation;

/// FCN wrapper that counts calls and operates in internal parameter space.
pub struct MnFcn<'a> {
    fcn: &'a dyn FCN,
    trafo: &'a MnUserTransformation,
    num_calls: Cell<usize>,
}

impl<'a> MnFcn<'a> {
    /// Create a new MnFcn wrapper.
    pub fn new(fcn: &'a dyn FCN, trafo: &'a MnUserTransformation) -> Self {
        Self {
            fcn,
            trafo,
            num_calls: Cell::new(0),
        }
    }

    /// Evaluate the function given internal-space parameters.
    /// Transforms to external space, then calls the user's FCN.
    pub fn call(&self, internal: &[f64]) -> f64 {
        self.num_calls.set(self.num_calls.get() + 1);
        let external = self.trafo.transform(internal);
        self.fcn.value(&external)
    }

    /// Evaluate an FCN with already transformed (external) parameters.
    pub fn call_with_transformed_params(&self, external: &[f64]) -> f64 {
        self.call_external(external)
    }

    /// Evaluate without applying the internal->external transformation.
    pub fn call_without_doing_trafo(&self, external: &[f64]) -> f64 {
        self.call_external(external)
    }

    /// Get the total number of function calls made.
    pub fn num_of_calls(&self) -> usize {
        self.num_calls.get()
    }

    /// Get the error definition from the user's FCN.
    pub fn error_def(&self) -> f64 {
        self.fcn.error_def()
    }

    /// Get the Up value (error definition).
    pub fn up(&self) -> f64 {
        self.error_def()
    }

    fn call_external(&self, external: &[f64]) -> f64 {
        self.num_calls.set(self.num_calls.get() + 1);
        self.fcn.value(external)
    }
}
