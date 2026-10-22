//! External / internal parameter mapping.
//!
//! Manages the index mapping between user (external) parameter indices and
//! optimizer (internal) indices. Fixed/const parameters are excluded from
//! the internal space. Bounded parameters go through transforms.
//!
//! Owns parameter transforms and internal/external index mappings.

use crate::parameter::MinuitParameter;
use crate::precision::MnMachinePrecision;
use crate::transform::{ParameterTransform, SinTransform, SqrtLowTransform, SqrtUpTransform};
use crate::user_covariance::MnUserCovariance;
use nalgebra::DMatrix;

#[derive(Debug, Clone)]
pub struct MnUserTransformation {
    precision: MnMachinePrecision,
    parameters: Vec<MinuitParameter>,
    /// Maps internal index → external index.
    int_of_ext: Vec<usize>,
    /// For each external param: Some(internal_index) if variable, None if fixed.
    ext_of_int: Vec<Option<usize>>,
    cache: Vec<f64>,
}

impl MnUserTransformation {
    pub fn new(params: Vec<MinuitParameter>) -> Self {
        let n = params.len();
        let mut ext_of_int = vec![None; n];
        let mut int_of_ext = Vec::new();

        for (ext, p) in params.iter().enumerate() {
            if !p.is_fixed() {
                ext_of_int[ext] = Some(int_of_ext.len());
                int_of_ext.push(ext);
            }
        }

        let cache = vec![0.0; n];

        Self {
            precision: MnMachinePrecision::new(),
            parameters: params,
            int_of_ext,
            ext_of_int,
            cache,
        }
    }

    pub fn precision(&self) -> &MnMachinePrecision {
        &self.precision
    }

    pub fn precision_mut(&mut self) -> &mut MnMachinePrecision {
        &mut self.precision
    }

    /// Number of variable (non-fixed) parameters = internal dimension.
    pub fn variable_parameters(&self) -> usize {
        self.int_of_ext.len()
    }

    /// Total number of parameters (including fixed).
    pub fn parameters_len(&self) -> usize {
        self.parameters.len()
    }

    pub fn parameters(&self) -> &[MinuitParameter] {
        &self.parameters
    }

    pub fn parameter(&self, ext: usize) -> &MinuitParameter {
        &self.parameters[ext]
    }

    pub fn parameter_mut(&mut self, ext: usize) -> &mut MinuitParameter {
        &mut self.parameters[ext]
    }

    pub fn params(&self) -> Vec<f64> {
        self.parameters.iter().map(|p| p.value()).collect()
    }

    pub fn errors(&self) -> Vec<f64> {
        self.parameters.iter().map(|p| p.error()).collect()
    }

    pub fn value(&self, ext: usize) -> f64 {
        self.parameters[ext].value()
    }

    pub fn error(&self, ext: usize) -> f64 {
        self.parameters[ext].error()
    }

    pub fn index(&self, name: &str) -> Option<usize> {
        self.parameters.iter().position(|p| p.name() == name)
    }

    pub fn find_index(&self, name: &str) -> Option<usize> {
        self.index(name)
    }

    pub fn set_value(&mut self, ext: usize, val: f64) {
        self.parameters[ext].set_value(val);
    }

    pub fn set_error(&mut self, ext: usize, err: f64) {
        self.parameters[ext].set_error(err);
    }

    pub fn set_limits(&mut self, ext: usize, lower: f64, upper: f64) {
        self.parameters[ext].set_limits(lower, upper);
    }

    pub fn set_lower_limit(&mut self, ext: usize, lower: f64) {
        self.parameters[ext].set_lower_limit(lower);
    }

    pub fn set_upper_limit(&mut self, ext: usize, upper: f64) {
        self.parameters[ext].set_upper_limit(upper);
    }

    pub fn remove_limits(&mut self, ext: usize) {
        self.parameters[ext].remove_limits();
    }

    pub fn set_name(&mut self, ext: usize, name: impl Into<String>) {
        self.parameters[ext].set_name(name);
    }

    pub fn set_precision(&mut self, eps: f64) {
        self.precision.set_precision(eps);
    }

    /// External index → internal index. Returns None if fixed.
    pub fn int_of_ext(&self, ext: usize) -> Option<usize> {
        self.ext_of_int[ext]
    }

    /// Internal index → external index.
    pub fn ext_of_int(&self, int: usize) -> usize {
        self.int_of_ext[int]
    }

    /// Transform a full internal vector to external values.
    pub fn transform(&self, internal: &[f64]) -> Vec<f64> {
        let mut result = Vec::with_capacity(self.parameters.len());
        for (ext, p) in self.parameters.iter().enumerate() {
            if p.is_fixed() {
                result.push(p.value());
            } else {
                let int_idx = self.ext_of_int[ext].unwrap();
                let int_val = internal[int_idx];
                result.push(self.int2ext(ext, int_val));
            }
        }
        result
    }

    /// Transform a single internal value to external.
    pub fn int2ext(&self, ext: usize, internal: f64) -> f64 {
        let p = &self.parameters[ext];
        if p.has_limits() {
            SinTransform.int2ext(internal, p.upper_limit(), p.lower_limit())
        } else if p.has_lower_limit() {
            SqrtLowTransform.int2ext(internal, p.upper_limit(), p.lower_limit())
        } else if p.has_upper_limit() {
            SqrtUpTransform.int2ext(internal, p.upper_limit(), p.lower_limit())
        } else {
            internal
        }
    }

    /// Transform a single external value to internal.
    pub fn ext2int(&self, ext: usize, value: f64) -> f64 {
        let p = &self.parameters[ext];
        if p.has_limits() {
            SinTransform.ext2int(value, p.upper_limit(), p.lower_limit(), &self.precision)
        } else if p.has_lower_limit() {
            SqrtLowTransform.ext2int(value, p.upper_limit(), p.lower_limit(), &self.precision)
        } else if p.has_upper_limit() {
            SqrtUpTransform.ext2int(value, p.upper_limit(), p.lower_limit(), &self.precision)
        } else {
            value
        }
    }

    /// Derivative d(external)/d(internal) for parameter `ext`.
    pub fn dint2ext(&self, ext: usize, internal: f64) -> f64 {
        let p = &self.parameters[ext];
        if p.has_limits() {
            SinTransform.dint2ext(internal, p.upper_limit(), p.lower_limit())
        } else if p.has_lower_limit() {
            SqrtLowTransform.dint2ext(internal, p.upper_limit(), p.lower_limit())
        } else if p.has_upper_limit() {
            SqrtUpTransform.dint2ext(internal, p.upper_limit(), p.lower_limit())
        } else {
            1.0
        }
    }

    /// Derivative d(internal)/d(external) for parameter `ext`.
    pub fn dext2int(&self, ext: usize, internal: f64) -> f64 {
        let d = self.dint2ext(ext, internal);
        if d.abs() > self.precision.eps2() {
            1.0 / d
        } else {
            0.0
        }
    }

    /// Internal error from external error, accounting for transform derivative.
    pub fn int2ext_error(&self, ext: usize, internal: f64, err: f64) -> f64 {
        let dx = self.dint2ext(ext, internal);
        if dx > 0.0 { err / dx } else { err }
    }

    pub fn int2ext_covariance(&self, internal: &[f64], cov: &DMatrix<f64>) -> MnUserCovariance {
        let n = self.variable_parameters();
        assert_eq!(internal.len(), n, "internal vector size mismatch");
        assert_eq!(cov.nrows(), n, "covariance row size mismatch");
        assert_eq!(cov.ncols(), n, "covariance col size mismatch");

        let jac: Vec<f64> = (0..n)
            .map(|int| {
                let ext = self.ext_of_int(int);
                self.dint2ext(ext, internal[int])
            })
            .collect();

        let mut data = Vec::with_capacity(n * (n + 1) / 2);
        for i in 0..n {
            for j in 0..=i {
                data.push(jac[i] * cov[(i, j)] * jac[j]);
            }
        }

        MnUserCovariance::from_vec(data, n)
    }

    pub fn ext2int_covariance(&self, internal: &[f64], cov: &DMatrix<f64>) -> MnUserCovariance {
        let n = self.variable_parameters();
        assert_eq!(internal.len(), n, "internal vector size mismatch");
        assert_eq!(cov.nrows(), n, "covariance row size mismatch");
        assert_eq!(cov.ncols(), n, "covariance col size mismatch");

        let jac_inv: Vec<f64> = (0..n)
            .map(|int| {
                let ext = self.ext_of_int(int);
                self.dext2int(ext, internal[int])
            })
            .collect();

        let mut data = Vec::with_capacity(n * (n + 1) / 2);
        for i in 0..n {
            for j in 0..=i {
                data.push(jac_inv[i] * cov[(i, j)] * jac_inv[j]);
            }
        }

        MnUserCovariance::from_vec(data, n)
    }

    /// Add a new variable parameter. Returns external index.
    pub fn add(&mut self, param: MinuitParameter) -> usize {
        let ext = self.parameters.len();
        let is_fixed = param.is_fixed();
        self.parameters.push(param);

        if is_fixed {
            self.ext_of_int.push(None);
        } else {
            self.ext_of_int.push(Some(self.int_of_ext.len()));
            self.int_of_ext.push(ext);
        }
        self.cache.push(0.0);
        ext
    }

    /// Fix parameter at external index. Removes from internal space.
    pub fn fix(&mut self, ext: usize) {
        self.parameters[ext].fix();
        self.rebuild_index();
    }

    /// Release parameter at external index. Adds back to internal space.
    pub fn release(&mut self, ext: usize) {
        self.parameters[ext].release();
        self.rebuild_index();
    }

    fn rebuild_index(&mut self) {
        self.int_of_ext.clear();
        for (ext, p) in self.parameters.iter().enumerate() {
            if !p.is_fixed() {
                self.ext_of_int[ext] = Some(self.int_of_ext.len());
                self.int_of_ext.push(ext);
            } else {
                self.ext_of_int[ext] = None;
            }
        }
    }

    /// Build internal parameter vector from current external values.
    pub fn initial_internal_values(&self) -> Vec<f64> {
        self.int_of_ext
            .iter()
            .map(|&ext| {
                let p = &self.parameters[ext];
                self.ext2int(ext, p.value())
            })
            .collect()
    }

    /// Build internal error vector from current external errors.
    pub fn initial_internal_errors(&self) -> Vec<f64> {
        self.int_of_ext
            .iter()
            .map(|&ext| {
                let p = &self.parameters[ext];
                let int_val = self.ext2int(ext, p.value());
                let dx = self.dint2ext(ext, int_val);
                if dx > 0.0 { p.error() / dx } else { p.error() }
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn variable_count() {
        let params = vec![
            MinuitParameter::new(0, "x", 1.0, 0.1),
            MinuitParameter::new(1, "y", 2.0, 0.2),
        ];
        let t = MnUserTransformation::new(params);
        assert_eq!(t.variable_parameters(), 2);
        assert_eq!(t.parameters_len(), 2);
    }

    #[test]
    fn fix_reduces_variable_count() {
        let params = vec![
            MinuitParameter::new(0, "x", 1.0, 0.1),
            MinuitParameter::new(1, "y", 2.0, 0.2),
        ];
        let mut t = MnUserTransformation::new(params);
        t.fix(0);
        assert_eq!(t.variable_parameters(), 1);
        assert_eq!(t.ext_of_int(0), 1); // internal 0 → external 1 (y)
    }

    #[test]
    fn bounded_transform_roundtrip() {
        let params = vec![MinuitParameter::with_limits(0, "x", 5.0, 0.1, 0.0, 10.0)];
        let t = MnUserTransformation::new(params);
        let int_val = t.ext2int(0, 5.0);
        let back = t.int2ext(0, int_val);
        assert!((back - 5.0).abs() < 1e-12);
    }

    #[test]
    fn unbounded_passthrough() {
        let pi = std::f64::consts::PI;
        let params = vec![MinuitParameter::new(0, "x", pi, 0.1)];
        let t = MnUserTransformation::new(params);
        assert!((t.ext2int(0, pi) - pi).abs() < 1e-15);
        assert!((t.int2ext(0, pi) - pi).abs() < 1e-15);
    }

    #[test]
    fn find_index_by_name() {
        let params = vec![
            MinuitParameter::new(0, "x", 1.0, 0.1),
            MinuitParameter::new(1, "y", 2.0, 0.2),
        ];
        let t = MnUserTransformation::new(params);
        assert_eq!(t.find_index("x"), Some(0));
        assert_eq!(t.find_index("y"), Some(1));
        assert_eq!(t.find_index("z"), None);
    }

    #[test]
    fn set_precision_updates_machine_precision() {
        let params = vec![MinuitParameter::new(0, "x", 1.0, 0.1)];
        let mut t = MnUserTransformation::new(params);
        t.set_precision(1.0e-12);
        assert!((t.precision().eps() - 1.0e-12).abs() < 1.0e-24);
    }

    #[test]
    fn int2ext_covariance_identity_for_unbounded() {
        let params = vec![
            MinuitParameter::new(0, "x", 1.0, 0.1),
            MinuitParameter::new(1, "y", 2.0, 0.2),
        ];
        let t = MnUserTransformation::new(params);
        let internal = t.initial_internal_values();
        let cov = DMatrix::from_row_slice(2, 2, &[1.0, 0.2, 0.2, 4.0]);
        let ucov = t.int2ext_covariance(&internal, &cov);
        assert!((ucov.get(0, 0) - 1.0).abs() < 1e-12);
        assert!((ucov.get(0, 1) - 0.2).abs() < 1e-12);
        assert!((ucov.get(1, 1) - 4.0).abs() < 1e-12);
    }
}
