//! Complete parameter state: parameters + function value + EDM + nfcn + covariance.
//!
//! State object returned to the user after minimization, containing fitted
//! values, errors, and optionally the covariance matrix.

use crate::parameter::MinuitParameter;
use crate::user_covariance::MnUserCovariance;
use crate::user_parameters::MnUserParameters;

/// External parameters, fit statistics, and optional covariance and correlation coefficients.
/// # Examples
///
/// ```
/// use minuit2::{MnUserParameters, MnUserParameterState};
/// let mut params = MnUserParameters::new();
/// params.add("x", 2.0, 0.1);
/// let state = MnUserParameterState::new(params);
/// assert_eq!(state.variable_parameters(), 1);
/// assert!(!state.has_covariance());
/// ```
#[derive(Debug, Clone)]
pub struct MnUserParameterState {
    params: MnUserParameters,
    covariance: Option<MnUserCovariance>,
    global_cc: Option<Vec<f64>>,
    fval: f64,
    edm: f64,
    nfcn: usize,
    is_valid: bool,
    covariance_valid: bool,
}

impl MnUserParameterState {
    /// Create initial state from user parameters (pre-minimization).
    #[must_use]
    pub fn new(params: MnUserParameters) -> Self {
        Self {
            params,
            covariance: None,
            global_cc: None,
            fval: 0.0,
            edm: 0.0,
            nfcn: 0,
            is_valid: false,
            covariance_valid: false,
        }
    }

    /// Return the stored user parameters.
    #[must_use]
    pub fn params(&self) -> &MnUserParameters {
        &self.params
    }

    /// Return mutable access to the user parameters.
    pub fn params_mut(&mut self) -> &mut MnUserParameters {
        &mut self.params
    }

    /// Get the function value.
    #[must_use]
    pub fn fval(&self) -> f64 {
        self.fval
    }

    /// Set the objective value.
    pub fn set_fval(&mut self, fval: f64) {
        self.fval = fval;
    }

    /// Get the estimated distance to minimum.
    #[must_use]
    pub fn edm(&self) -> f64 {
        self.edm
    }

    /// Set the estimated distance to the minimum, in objective units.
    pub fn set_edm(&mut self, edm: f64) {
        self.edm = edm;
    }

    /// Get the total number of function calls.
    #[must_use]
    pub fn nfcn(&self) -> usize {
        self.nfcn
    }

    /// Set the number of objective evaluations.
    pub fn set_nfcn(&mut self, nfcn: usize) {
        self.nfcn = nfcn;
    }

    /// Check if the minimization result is valid.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.is_valid
    }

    /// Set the stored validity flag.
    pub fn set_valid(&mut self, valid: bool) {
        self.is_valid = valid;
    }

    /// Check if the covariance matrix is available.
    #[must_use]
    pub fn has_covariance(&self) -> bool {
        self.covariance_valid
    }

    /// Get the covariance matrix if available.
    #[must_use]
    pub fn covariance(&self) -> Option<&MnUserCovariance> {
        self.covariance.as_ref()
    }

    /// Store a covariance matrix and mark it available.
    pub fn set_covariance(&mut self, cov: MnUserCovariance) {
        self.covariance_valid = true;
        self.covariance = Some(cov);
    }

    /// Add to a covariance entry, allocating a zero matrix for variable parameters if needed.
    ///
    /// # Panics
    /// Panics if the computed packed index is outside the covariance storage.
    pub fn add_covariance(&mut self, i: usize, j: usize, value: f64) {
        if self.covariance.is_none() {
            self.covariance = Some(MnUserCovariance::new(self.params.variable_parameters()));
            self.covariance_valid = true;
        }
        if let Some(cov) = self.covariance.as_mut() {
            let current = cov.get(i, j);
            cov.set(i, j, current + value);
        }
    }

    /// Return one when covariance is marked available, otherwise zero.
    #[must_use]
    pub fn covariance_status(&self) -> i32 {
        i32::from(self.covariance_valid)
    }

    /// Return a clone of the stored covariance matrix; this method does not invert it.
    #[must_use]
    pub fn hessian(&self) -> Option<MnUserCovariance> {
        self.covariance.clone()
    }

    /// Get the global correlation coefficients if available.
    #[must_use]
    pub fn global_cc(&self) -> Option<&[f64]> {
        self.global_cc.as_deref()
    }

    /// Store the dimensionless global correlation coefficients.
    pub fn set_global_cc(&mut self, gcc: Vec<f64>) {
        self.global_cc = Some(gcc);
    }

    // --- Delegation to MnUserParameters ---

    /// Add a free parameter with an external value and initial error; return its external index.
    pub fn add(&mut self, name: impl Into<String>, value: f64, error: f64) -> usize {
        self.params.add(name, value, error)
    }

    /// Add a parameter with two external bounds; return its external index.
    pub fn add_limited(
        &mut self,
        name: impl Into<String>,
        value: f64,
        error: f64,
        lower: f64,
        upper: f64,
    ) -> usize {
        self.params.add_limited(name, value, error, lower, upper)
    }

    /// Exclude the parameter from subsequent minimization steps.
    ///
    /// # Panics
    /// Panics if the external parameter index is out of range.
    pub fn fix(&mut self, ext: usize) {
        self.params.fix(ext);
    }

    /// Release a fixed parameter unless it is constant.
    ///
    /// # Panics
    /// Panics if the external parameter index is out of range.
    pub fn release(&mut self, ext: usize) {
        self.params.release(ext);
    }

    /// Set the parameter value in external coordinates.
    ///
    /// # Panics
    /// Panics if the external parameter index is out of range.
    pub fn set_value(&mut self, ext: usize, val: f64) {
        self.params.set_value(ext, val);
    }

    /// Set the parameter error in external coordinates.
    ///
    /// # Panics
    /// Panics if the external parameter index is out of range.
    pub fn set_error(&mut self, ext: usize, err: f64) {
        self.params.set_error(ext, err);
    }

    /// Set both bounds in external coordinates.
    ///
    /// # Panics
    /// Panics if the external parameter index is out of range.
    /// Panics if the lower bound is not less than the upper bound.
    pub fn set_limits(&mut self, ext: usize, lower: f64, upper: f64) {
        self.params.set_limits(ext, lower, upper);
    }

    /// Set the lower bound in external coordinates.
    ///
    /// # Panics
    /// Panics if the external parameter index is out of range.
    pub fn set_lower_limit(&mut self, ext: usize, lower: f64) {
        self.params.set_lower_limit(ext, lower);
    }

    /// Set the upper bound in external coordinates.
    ///
    /// # Panics
    /// Panics if the external parameter index is out of range.
    pub fn set_upper_limit(&mut self, ext: usize, upper: f64) {
        self.params.set_upper_limit(ext, upper);
    }

    /// Remove both parameter bounds.
    ///
    /// # Panics
    /// Panics if the external parameter index is out of range.
    pub fn remove_limits(&mut self, ext: usize) {
        self.params.remove_limits(ext);
    }

    /// Set the parameter name.
    ///
    /// # Panics
    /// Panics if the external parameter index is out of range.
    pub fn set_name(&mut self, ext: usize, name: impl Into<String>) {
        self.params.set_name(ext, name);
    }

    /// Set the epsilon used to compute numerical transformation steps.
    pub fn set_precision(&mut self, eps: f64) {
        self.params.set_precision(eps);
    }

    /// Return the parameter value in external coordinates.
    #[must_use]
    pub fn value(&self, name: &str) -> Option<f64> {
        self.params.value(name)
    }

    /// Return the parameter error in external coordinates.
    #[must_use]
    pub fn error(&self, name: &str) -> Option<f64> {
        self.params.error(name)
    }

    /// Return parameter errors in external order and external units.
    #[must_use]
    pub fn errors(&self) -> Vec<f64> {
        self.params.errors()
    }

    /// Find a parameter by name, returning its external index if present.
    #[must_use]
    pub fn index(&self, name: &str) -> Option<usize> {
        self.params.index(name)
    }

    /// Convert one internal coordinate to its external value.
    ///
    /// # Panics
    /// Panics if the internal parameter index is out of range.
    #[must_use]
    pub fn int2ext(&self, int: usize, internal: f64) -> f64 {
        let ext = self.params.trafo().ext_of_int(int);
        self.params.trafo().int2ext(ext, internal)
    }

    /// Convert an external value to its internal coordinate.
    ///
    /// # Panics
    /// Panics if the external parameter index is out of range.
    #[must_use]
    pub fn ext2int(&self, ext: usize, value: f64) -> f64 {
        self.params.trafo().ext2int(ext, value)
    }

    /// Return the internal index of a variable external parameter, or `None` for a fixed or constant parameter.
    ///
    /// # Panics
    /// Panics if the external parameter index is out of range.
    #[must_use]
    pub fn int_of_ext(&self, ext: usize) -> Option<usize> {
        self.params.trafo().int_of_ext(ext)
    }

    /// Return the external index corresponding to an internal index.
    ///
    /// # Panics
    /// Panics if the internal parameter index is out of range.
    #[must_use]
    pub fn ext_of_int(&self, int: usize) -> usize {
        self.params.trafo().ext_of_int(int)
    }

    /// Return a parameter by external index.
    ///
    /// # Panics
    /// Panics if the external parameter index is out of range.
    #[must_use]
    pub fn parameter(&self, ext: usize) -> &MinuitParameter {
        self.params.trafo().parameter(ext)
    }

    /// Return the number of variable parameters.
    #[must_use]
    pub fn variable_parameters(&self) -> usize {
        self.params.variable_parameters()
    }

    /// Return the number of parameters, including fixed and constant parameters.
    #[must_use]
    pub fn len(&self) -> usize {
        self.params.len()
    }

    /// Return whether there are no parameters.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.params.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn state_delegates_limit_operations() {
        let mut params = MnUserParameters::new();
        params.add("x", 1.0, 0.1);
        let mut state = MnUserParameterState::new(params);

        state.set_lower_limit(0, -1.0);
        state.set_upper_limit(0, 2.0);
        let p = state.parameter(0);
        assert!(p.has_lower_limit());
        assert!(p.has_upper_limit());
    }

    #[test]
    fn state_name_and_index_roundtrip() {
        let mut params = MnUserParameters::new();
        params.add("x", 1.0, 0.1);
        let mut state = MnUserParameterState::new(params);
        state.set_name(0, "alpha");
        assert_eq!(state.index("alpha"), Some(0));
        assert_eq!(state.index("x"), None);
    }

    #[test]
    fn state_internal_external_mapping() {
        let mut params = MnUserParameters::new();
        params.add("x", 3.0, 0.1);
        let state = MnUserParameterState::new(params);

        assert_eq!(state.int_of_ext(0), Some(0));
        assert_eq!(state.ext_of_int(0), 0);
        assert!((state.ext2int(0, 3.0) - 3.0).abs() < 1e-15);
        assert!((state.int2ext(0, 3.0) - 3.0).abs() < 1e-15);
    }
}
