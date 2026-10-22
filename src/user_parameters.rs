//! High-level parameter collection with name-based lookup.
//!
//! Delegates index mapping and transforms to `MnUserTransformation`. Provides
//! the user-facing API for adding, fixing, releasing, and bounding parameters.

use std::collections::HashMap;

use crate::parameter::MinuitParameter;
use crate::user_transformation::MnUserTransformation;

#[derive(Debug, Clone)]
pub struct MnUserParameters {
    trafo: MnUserTransformation,
    name_map: HashMap<String, usize>,
}

impl MnUserParameters {
    /// Create an empty parameter collection.
    pub fn new() -> Self {
        Self {
            trafo: MnUserTransformation::new(Vec::new()),
            name_map: HashMap::new(),
        }
    }

    /// Get the internal transformation object.
    pub fn trafo(&self) -> &MnUserTransformation {
        &self.trafo
    }

    pub fn trafo_mut(&mut self) -> &mut MnUserTransformation {
        &mut self.trafo
    }

    /// Add a free parameter. Returns external index.
    pub fn add(&mut self, name: impl Into<String>, value: f64, error: f64) -> usize {
        let name = name.into();
        let ext = self.trafo.parameters_len();
        let param = MinuitParameter::new(ext, &name, value, error);
        self.trafo.add(param);
        self.name_map.insert(name, ext);
        ext
    }

    /// Add a parameter with both bounds.
    pub fn add_limited(
        &mut self,
        name: impl Into<String>,
        value: f64,
        error: f64,
        lower: f64,
        upper: f64,
    ) -> usize {
        let name = name.into();
        let ext = self.trafo.parameters_len();
        let param = MinuitParameter::with_limits(ext, &name, value, error, lower, upper);
        self.trafo.add(param);
        self.name_map.insert(name, ext);
        ext
    }

    /// Add a parameter with lower bound only.
    pub fn add_lower_limited(
        &mut self,
        name: impl Into<String>,
        value: f64,
        error: f64,
        lower: f64,
    ) -> usize {
        let name = name.into();
        let ext = self.trafo.parameters_len();
        let param = MinuitParameter::with_lower_limit(ext, &name, value, error, lower);
        self.trafo.add(param);
        self.name_map.insert(name, ext);
        ext
    }

    /// Add a parameter with upper bound only.
    pub fn add_upper_limited(
        &mut self,
        name: impl Into<String>,
        value: f64,
        error: f64,
        upper: f64,
    ) -> usize {
        let name = name.into();
        let ext = self.trafo.parameters_len();
        let param = MinuitParameter::with_upper_limit(ext, &name, value, error, upper);
        self.trafo.add(param);
        self.name_map.insert(name, ext);
        ext
    }

    /// Add a constant parameter (fixed, never released).
    pub fn add_const(&mut self, name: impl Into<String>, value: f64) -> usize {
        let name = name.into();
        let ext = self.trafo.parameters_len();
        let param = MinuitParameter::constant(ext, &name, value);
        self.trafo.add(param);
        self.name_map.insert(name, ext);
        ext
    }

    /// Fix parameter by external index.
    pub fn fix(&mut self, ext: usize) {
        self.trafo.fix(ext);
    }

    /// Release parameter by external index.
    pub fn release(&mut self, ext: usize) {
        self.trafo.release(ext);
    }

    /// Set value by external index.
    pub fn set_value(&mut self, ext: usize, val: f64) {
        self.trafo.parameter_mut(ext).set_value(val);
    }

    /// Set error by external index.
    pub fn set_error(&mut self, ext: usize, err: f64) {
        self.trafo.parameter_mut(ext).set_error(err);
    }

    /// Set limits by external index.
    pub fn set_limits(&mut self, ext: usize, lower: f64, upper: f64) {
        self.trafo.parameter_mut(ext).set_limits(lower, upper);
    }

    pub fn set_lower_limit(&mut self, ext: usize, lower: f64) {
        self.trafo.parameter_mut(ext).set_lower_limit(lower);
    }

    pub fn set_upper_limit(&mut self, ext: usize, upper: f64) {
        self.trafo.parameter_mut(ext).set_upper_limit(upper);
    }

    /// Remove limits by external index.
    pub fn remove_limits(&mut self, ext: usize) {
        self.trafo.parameter_mut(ext).remove_limits();
    }

    pub fn set_name(&mut self, ext: usize, name: impl Into<String>) {
        let old = self.trafo.parameter(ext).name().to_string();
        let new = name.into();
        self.trafo.parameter_mut(ext).set_name(new.clone());
        self.name_map.remove(&old);
        self.name_map.insert(new, ext);
    }

    pub fn set_precision(&mut self, eps: f64) {
        self.trafo.precision_mut().set_precision(eps);
    }

    /// Lookup external index by name.
    pub fn index(&self, name: &str) -> Option<usize> {
        self.name_map.get(name).copied()
    }

    /// Get parameter by name.
    pub fn parameter(&self, name: &str) -> Option<&MinuitParameter> {
        self.name_map.get(name).map(|&i| self.trafo.parameter(i))
    }

    /// Get parameter value by name.
    pub fn value(&self, name: &str) -> Option<f64> {
        self.parameter(name).map(|p| p.value())
    }

    /// Get parameter error by name.
    pub fn error(&self, name: &str) -> Option<f64> {
        self.parameter(name).map(|p| p.error())
    }

    pub fn errors(&self) -> Vec<f64> {
        self.trafo.parameters().iter().map(|p| p.error()).collect()
    }

    /// Number of total parameters.
    pub fn len(&self) -> usize {
        self.trafo.parameters_len()
    }

    /// Whether there are no parameters.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Number of variable (non-fixed) parameters.
    pub fn variable_parameters(&self) -> usize {
        self.trafo.variable_parameters()
    }

    /// All parameter references.
    pub fn params(&self) -> &[MinuitParameter] {
        self.trafo.parameters()
    }
}

impl Default for MnUserParameters {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_and_lookup() {
        let mut p = MnUserParameters::new();
        p.add("x", 1.0, 0.1);
        p.add("y", 2.0, 0.2);
        assert_eq!(p.len(), 2);
        assert_eq!(p.index("x"), Some(0));
        assert_eq!(p.index("y"), Some(1));
        assert!((p.value("x").unwrap() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn fix_reduces_variable() {
        let mut p = MnUserParameters::new();
        p.add("x", 1.0, 0.1);
        p.add("y", 2.0, 0.2);
        assert_eq!(p.variable_parameters(), 2);
        p.fix(0);
        assert_eq!(p.variable_parameters(), 1);
        p.release(0);
        assert_eq!(p.variable_parameters(), 2);
    }

    #[test]
    fn set_value_and_error() {
        let mut p = MnUserParameters::new();
        p.add("x", 1.0, 0.1);
        p.set_value(0, 42.0);
        p.set_error(0, 0.5);
        assert!((p.value("x").unwrap() - 42.0).abs() < 1e-15);
        assert!((p.error("x").unwrap() - 0.5).abs() < 1e-15);
    }

    #[test]
    fn set_name_updates_lookup_map() {
        let mut p = MnUserParameters::new();
        p.add("x", 1.0, 0.1);
        p.set_name(0, "alpha");
        assert_eq!(p.index("x"), None);
        assert_eq!(p.index("alpha"), Some(0));
    }

    #[test]
    fn set_lower_and_upper_limits() {
        let mut p = MnUserParameters::new();
        p.add("x", 1.0, 0.1);
        p.set_lower_limit(0, -2.0);
        p.set_upper_limit(0, 3.0);
        let x = p.parameter("x").expect("x must exist");
        assert!(x.has_lower_limit());
        assert!(x.has_upper_limit());
        assert!((x.lower_limit() + 2.0).abs() < 1e-15);
        assert!((x.upper_limit() - 3.0).abs() < 1e-15);
    }

    #[test]
    fn set_precision_propagates_to_transformation() {
        let mut p = MnUserParameters::new();
        p.add("x", 1.0, 0.1);
        p.set_precision(1.0e-12);
        assert!((p.trafo().precision().eps() - 1.0e-12).abs() < 1.0e-24);
    }
}
