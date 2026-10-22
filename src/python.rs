use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::exceptions::{
    PyIndexError, PyKeyError, PyNotImplementedError, PyRuntimeError, PyTypeError, PyValueError,
};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use std::collections::{HashMap, HashSet};

use crate::{FCN, FunctionMinimum, MnContours, MnHesse, MnMigrad, MnMinos, MnSimplex};

// Aliases for the numpy array handles returned by the scan/profile/contour
// methods (keeps their signatures readable and clippy::type_complexity quiet).
type Arr1f<'py> = Bound<'py, PyArray1<f64>>;
type Arr1b<'py> = Bound<'py, PyArray1<bool>>;
type Arr2f<'py> = Bound<'py, PyArray2<f64>>;

// ============================================================================
// FCN Wrapper
// ============================================================================

struct PythonFCN {
    fcn: Py<PyAny>,
    errordef: f64,
}

impl FCN for PythonFCN {
    fn error_def(&self) -> f64 {
        self.errordef
    }

    fn value(&self, par: &[f64]) -> f64 {
        Python::attach(|py| {
            let Ok(args) = PyTuple::new(py, par) else {
                return f64::INFINITY;
            };
            match self.fcn.call(py, &args, None) {
                Ok(val) => val.extract::<f64>(py).unwrap_or(f64::INFINITY),
                Err(e) => {
                    e.print(py);
                    f64::INFINITY
                }
            }
        })
    }
}

// ============================================================================
// Minuit Class
// ============================================================================

struct StoredMError {
    name: String,
    lower: f64,
    upper: f64,
    is_valid: bool,
    lower_valid: bool,
    upper_valid: bool,
    at_lower_limit: bool,
    at_upper_limit: bool,
    at_lower_max_fcn: bool,
    at_upper_max_fcn: bool,
    lower_new_min: bool,
    upper_new_min: bool,
    nfcn: usize,
    min: f64,
}

#[pyclass(frozen, skip_from_py_object)]
#[derive(Clone)]
struct FMin {
    fval: f64,
    edm: f64,
    nfcn: usize,
    errordef: f64,
    is_valid: bool,
    has_valid_parameters: bool,
    has_covariance: bool,
    has_made_posdef_covar: bool,
    is_above_max_edm: bool,
    has_reached_call_limit: bool,
}

#[pymethods]
impl FMin {
    #[getter]
    fn fval(&self) -> f64 {
        self.fval
    }

    #[getter]
    fn edm(&self) -> f64 {
        self.edm
    }

    #[getter]
    fn nfcn(&self) -> usize {
        self.nfcn
    }

    #[getter]
    fn errordef(&self) -> f64 {
        self.errordef
    }

    #[getter]
    fn is_valid(&self) -> bool {
        self.is_valid
    }

    #[getter]
    fn has_valid_parameters(&self) -> bool {
        self.has_valid_parameters
    }

    #[getter]
    fn has_covariance(&self) -> bool {
        self.has_covariance
    }

    #[getter]
    fn has_made_posdef_covar(&self) -> bool {
        self.has_made_posdef_covar
    }

    #[getter]
    fn is_above_max_edm(&self) -> bool {
        self.is_above_max_edm
    }

    #[getter]
    fn has_reached_call_limit(&self) -> bool {
        self.has_reached_call_limit
    }
}

impl FMin {
    fn from_minimum(m: &FunctionMinimum) -> Self {
        Self {
            fval: m.fval(),
            edm: m.edm(),
            nfcn: m.nfcn(),
            errordef: m.up(),
            is_valid: m.is_valid(),
            has_valid_parameters: m.has_valid_parameters(),
            has_covariance: m.user_state().covariance().is_some(),
            has_made_posdef_covar: m.has_made_pos_def_covar(),
            is_above_max_edm: m.is_above_max_edm(),
            has_reached_call_limit: m.reached_call_limit(),
        }
    }
}

#[pyclass(frozen)]
struct MError {
    name: String,
    lower: f64,
    upper: f64,
    is_valid: bool,
    lower_valid: bool,
    upper_valid: bool,
    at_lower_limit: bool,
    at_upper_limit: bool,
    at_lower_max_fcn: bool,
    at_upper_max_fcn: bool,
    lower_new_min: bool,
    upper_new_min: bool,
    nfcn: usize,
    min: f64,
}

#[pymethods]
impl MError {
    #[getter]
    fn name(&self) -> &str {
        &self.name
    }

    #[getter]
    fn lower(&self) -> f64 {
        self.lower
    }

    #[getter]
    fn upper(&self) -> f64 {
        self.upper
    }

    #[getter]
    fn is_valid(&self) -> bool {
        self.is_valid
    }

    #[getter]
    fn lower_valid(&self) -> bool {
        self.lower_valid
    }

    #[getter]
    fn upper_valid(&self) -> bool {
        self.upper_valid
    }

    #[getter]
    fn at_lower_limit(&self) -> bool {
        self.at_lower_limit
    }

    #[getter]
    fn at_upper_limit(&self) -> bool {
        self.at_upper_limit
    }

    #[getter]
    fn at_lower_max_fcn(&self) -> bool {
        self.at_lower_max_fcn
    }

    #[getter]
    fn at_upper_max_fcn(&self) -> bool {
        self.at_upper_max_fcn
    }

    #[getter]
    fn lower_new_min(&self) -> bool {
        self.lower_new_min
    }

    #[getter]
    fn upper_new_min(&self) -> bool {
        self.upper_new_min
    }

    #[getter]
    fn nfcn(&self) -> usize {
        self.nfcn
    }

    #[getter]
    fn min(&self) -> f64 {
        self.min
    }
}

impl MError {
    fn from_stored(s: &StoredMError) -> Self {
        Self {
            name: s.name.clone(),
            lower: s.lower,
            upper: s.upper,
            is_valid: s.is_valid,
            lower_valid: s.lower_valid,
            upper_valid: s.upper_valid,
            at_lower_limit: s.at_lower_limit,
            at_upper_limit: s.at_upper_limit,
            at_lower_max_fcn: s.at_lower_max_fcn,
            at_upper_max_fcn: s.at_upper_max_fcn,
            lower_new_min: s.lower_new_min,
            upper_new_min: s.upper_new_min,
            nfcn: s.nfcn,
            min: s.min,
        }
    }
}

#[pyclass(frozen)]
struct Param {
    number: usize,
    name: String,
    value: f64,
    error: f64,
    is_fixed: bool,
    lower_limit: Option<f64>,
    upper_limit: Option<f64>,
    has_limits: bool,
}

#[pymethods]
impl Param {
    #[getter]
    fn number(&self) -> usize {
        self.number
    }

    #[getter]
    fn name(&self) -> &str {
        &self.name
    }

    #[getter]
    fn value(&self) -> f64 {
        self.value
    }

    #[getter]
    fn error(&self) -> f64 {
        self.error
    }

    #[getter]
    fn is_fixed(&self) -> bool {
        self.is_fixed
    }

    #[getter]
    fn lower_limit(&self) -> Option<f64> {
        self.lower_limit
    }

    #[getter]
    fn upper_limit(&self) -> Option<f64> {
        self.upper_limit
    }

    #[getter]
    fn has_limits(&self) -> bool {
        self.has_limits
    }
}

// ============================================================================
// Parameter views (iminuit-style indexable, mutable, persisting accessors)
// ============================================================================

/// Map a non-finite limit (±inf) or `None` to `None` (= unbounded that side).
fn finite_or_none(v: Option<f64>) -> Option<f64> {
    match v {
        Some(x) if x.is_finite() => Some(x),
        _ => None,
    }
}

/// Resolve an int (incl. negative) or str key to a parameter name.
fn resolve_param_name(names: &[String], key: &Bound<'_, PyAny>) -> PyResult<String> {
    if let Ok(i) = key.extract::<isize>() {
        let n = names.len() as isize;
        let idx = if i < 0 { i + n } else { i };
        if idx < 0 || idx >= n {
            return Err(PyIndexError::new_err("parameter index out of range"));
        }
        Ok(names[idx as usize].clone())
    } else if let Ok(s) = key.extract::<String>() {
        if names.iter().any(|x| x == &s) {
            Ok(s)
        } else {
            Err(PyKeyError::new_err(s))
        }
    } else {
        Err(PyTypeError::new_err("parameter key must be int or str"))
    }
}

fn list_iter(list: Bound<'_, PyList>) -> PyResult<Py<PyAny>> {
    Ok(list.into_any().try_iter()?.into_any().unbind())
}

#[pyclass]
struct ValueView {
    owner: Py<Minuit>,
}

#[pymethods]
impl ValueView {
    fn __getitem__(&self, py: Python<'_>, key: Bound<'_, PyAny>) -> PyResult<f64> {
        let m = self.owner.borrow(py);
        let name = resolve_param_name(&m.names, &key)?;
        Ok(*m.values.get(&name).unwrap_or(&0.0))
    }
    fn __setitem__(&self, py: Python<'_>, key: Bound<'_, PyAny>, value: f64) -> PyResult<()> {
        let mut m = self.owner.borrow_mut(py);
        let name = resolve_param_name(&m.names, &key)?;
        m.values.insert(name, value);
        Ok(())
    }
    fn __len__(&self, py: Python<'_>) -> usize {
        self.owner.borrow(py).names.len()
    }
    fn __iter__(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let m = self.owner.borrow(py);
        let vals: Vec<f64> = m
            .names
            .iter()
            .map(|n| *m.values.get(n).unwrap_or(&0.0))
            .collect();
        list_iter(PyList::new(py, vals)?)
    }
    fn to_dict(&self, py: Python<'_>) -> HashMap<String, f64> {
        let m = self.owner.borrow(py);
        m.names
            .iter()
            .map(|n| (n.clone(), *m.values.get(n).unwrap_or(&0.0)))
            .collect()
    }
    fn __repr__(&self, py: Python<'_>) -> String {
        format!("<ValueView {:?}>", self.to_dict(py))
    }
}

#[pyclass]
struct ErrorView {
    owner: Py<Minuit>,
}

#[pymethods]
impl ErrorView {
    fn __getitem__(&self, py: Python<'_>, key: Bound<'_, PyAny>) -> PyResult<f64> {
        let m = self.owner.borrow(py);
        let name = resolve_param_name(&m.names, &key)?;
        Ok(*m.errors.get(&name).unwrap_or(&0.1))
    }
    fn __setitem__(&self, py: Python<'_>, key: Bound<'_, PyAny>, value: f64) -> PyResult<()> {
        let mut m = self.owner.borrow_mut(py);
        let name = resolve_param_name(&m.names, &key)?;
        m.errors.insert(name, value);
        Ok(())
    }
    fn __len__(&self, py: Python<'_>) -> usize {
        self.owner.borrow(py).names.len()
    }
    fn __iter__(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let m = self.owner.borrow(py);
        let vals: Vec<f64> = m
            .names
            .iter()
            .map(|n| *m.errors.get(n).unwrap_or(&0.1))
            .collect();
        list_iter(PyList::new(py, vals)?)
    }
    fn to_dict(&self, py: Python<'_>) -> HashMap<String, f64> {
        let m = self.owner.borrow(py);
        m.names
            .iter()
            .map(|n| (n.clone(), *m.errors.get(n).unwrap_or(&0.1)))
            .collect()
    }
    fn __repr__(&self, py: Python<'_>) -> String {
        format!("<ErrorView {:?}>", self.to_dict(py))
    }
}

#[pyclass]
struct FixedView {
    owner: Py<Minuit>,
}

#[pymethods]
impl FixedView {
    fn __getitem__(&self, py: Python<'_>, key: Bound<'_, PyAny>) -> PyResult<bool> {
        let m = self.owner.borrow(py);
        let name = resolve_param_name(&m.names, &key)?;
        Ok(m.fixed.contains(&name))
    }
    fn __setitem__(&self, py: Python<'_>, key: Bound<'_, PyAny>, value: bool) -> PyResult<()> {
        let mut m = self.owner.borrow_mut(py);
        let name = resolve_param_name(&m.names, &key)?;
        if value {
            m.fixed.insert(name);
        } else {
            m.fixed.remove(&name);
        }
        Ok(())
    }
    fn __len__(&self, py: Python<'_>) -> usize {
        self.owner.borrow(py).names.len()
    }
    fn __iter__(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let m = self.owner.borrow(py);
        let vals: Vec<bool> = m.names.iter().map(|n| m.fixed.contains(n)).collect();
        list_iter(PyList::new(py, vals)?)
    }
    fn to_dict(&self, py: Python<'_>) -> HashMap<String, bool> {
        let m = self.owner.borrow(py);
        m.names
            .iter()
            .map(|n| (n.clone(), m.fixed.contains(n)))
            .collect()
    }
    fn __repr__(&self, py: Python<'_>) -> String {
        format!("<FixedView {:?}>", self.to_dict(py))
    }
}

#[pyclass]
struct LimitView {
    owner: Py<Minuit>,
}

impl LimitView {
    fn pair(m: &Minuit, name: &str) -> (f64, f64) {
        match m.limits.get(name) {
            Some((lo, hi)) => (lo.unwrap_or(f64::NEG_INFINITY), hi.unwrap_or(f64::INFINITY)),
            None => (f64::NEG_INFINITY, f64::INFINITY),
        }
    }
}

#[pymethods]
impl LimitView {
    fn __getitem__(&self, py: Python<'_>, key: Bound<'_, PyAny>) -> PyResult<(f64, f64)> {
        let m = self.owner.borrow(py);
        let name = resolve_param_name(&m.names, &key)?;
        Ok(Self::pair(&m, &name))
    }
    fn __setitem__(
        &self,
        py: Python<'_>,
        key: Bound<'_, PyAny>,
        value: Bound<'_, PyAny>,
    ) -> PyResult<()> {
        let (lo, hi) = if value.is_none() {
            (None, None)
        } else {
            let (l, u): (Option<f64>, Option<f64>) = value.extract()?;
            (finite_or_none(l), finite_or_none(u))
        };
        let mut m = self.owner.borrow_mut(py);
        let name = resolve_param_name(&m.names, &key)?;
        if lo.is_none() && hi.is_none() {
            m.limits.remove(&name);
        } else {
            m.limits.insert(name, (lo, hi));
        }
        Ok(())
    }
    fn __len__(&self, py: Python<'_>) -> usize {
        self.owner.borrow(py).names.len()
    }
    fn __iter__(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let m = self.owner.borrow(py);
        let vals: Vec<(f64, f64)> = m.names.iter().map(|n| Self::pair(&m, n)).collect();
        list_iter(PyList::new(py, vals)?)
    }
    fn to_dict(&self, py: Python<'_>) -> HashMap<String, (f64, f64)> {
        let m = self.owner.borrow(py);
        m.names
            .iter()
            .map(|n| (n.clone(), Self::pair(&m, n)))
            .collect()
    }
    fn __repr__(&self, py: Python<'_>) -> String {
        format!("<LimitView {:?}>", self.to_dict(py))
    }
}

#[pyclass(name = "Minuit")]
struct Minuit {
    fcn: Py<PyAny>,
    names: Vec<String>,
    values: HashMap<String, f64>,
    errors: HashMap<String, f64>,
    limits: HashMap<String, (Option<f64>, Option<f64>)>,
    fixed: HashSet<String>,
    last_minimum: Option<FunctionMinimum>,
    strategy: u32,
    tolerance: f64,
    max_calls: Option<usize>,
    errordef: f64,
    merrors: HashMap<String, StoredMError>,
    init_values: HashMap<String, f64>,
    init_errors: HashMap<String, f64>,
    init_fixed: HashSet<String>,
    init_limits: HashMap<String, (Option<f64>, Option<f64>)>,
    scan_fmin: Option<FMin>,
}

#[pymethods]
impl Minuit {
    #[classattr]
    const LEAST_SQUARES: f64 = 1.0;
    #[classattr]
    const LIKELIHOOD: f64 = 0.5;

    #[new]
    #[pyo3(signature = (fcn, *args, name=None, **kwds))]
    fn new(
        py: Python<'_>,
        fcn: Py<PyAny>,
        args: Vec<f64>,
        name: Option<Vec<String>>,
        kwds: Option<Bound<'_, PyDict>>,
    ) -> PyResult<Self> {
        let names = if let Some(name) = name {
            name
        } else {
            let mut introspected = Vec::new();
            if let Ok(inspect) = py.import("inspect")
                && let Ok(sig) = inspect.call_method1("signature", (&fcn,))
                && let Ok(params_map) = sig.getattr("parameters")
                && let Ok(keys) = params_map.call_method0("keys")
                && let Ok(iter) = keys.try_iter()
            {
                for param_name in iter {
                    let param_name = param_name?;
                    let param = params_map.call_method1("__getitem__", (&param_name,))?;
                    let kind = param.getattr("kind")?.extract::<u8>()?;
                    if kind != 2 && kind != 4 {
                        introspected.push(param_name.extract::<String>()?);
                    }
                }
            }
            if introspected.is_empty() {
                let mut fallback = Vec::new();
                if let Some(p) = &kwds {
                    for (key, _) in p.iter() {
                        fallback.push(key.extract::<String>()?);
                    }
                }
                fallback
            } else {
                introspected
            }
        };
        let mut values = HashMap::new();
        let mut errors = HashMap::new();
        let fixed = HashSet::new();
        let limits = HashMap::new();

        for name in &names {
            values.insert(name.clone(), 0.0);
            errors.insert(name.clone(), 0.1);
        }

        if args.len() > names.len() {
            return Err(PyValueError::new_err("too many positional arguments"));
        }
        for (i, val) in args.into_iter().enumerate() {
            values.insert(names[i].clone(), val);
        }

        if let Some(p) = kwds {
            for (name, value) in p.iter() {
                let name_str = name.extract::<String>()?;
                if !values.contains_key(&name_str) {
                    return Err(PyValueError::new_err(format!(
                        "unknown parameter: {}",
                        name_str
                    )));
                }
                let val = value.extract::<f64>()?;
                values.insert(name_str, val);
            }
        }

        let init_values = values.clone();
        let init_errors = errors.clone();
        let init_fixed = fixed.clone();
        let init_limits = limits.clone();

        Ok(Minuit {
            fcn,
            names,
            values,
            errors,
            limits,
            fixed,
            last_minimum: None,
            strategy: 1,
            tolerance: 0.1,
            max_calls: None,
            errordef: 1.0,
            merrors: HashMap::new(),
            init_values,
            init_errors,
            init_fixed,
            init_limits,
            scan_fmin: None,
        })
    }

    #[getter]
    fn get_errordef(&self) -> f64 {
        self.errordef
    }

    #[setter]
    fn set_errordef(&mut self, v: f64) {
        self.errordef = v;
    }

    #[getter]
    fn get_strategy(&self) -> u32 {
        self.strategy
    }

    #[setter]
    fn set_strategy(&mut self, v: u32) -> PyResult<()> {
        if v > 2 {
            return Err(PyValueError::new_err("strategy must be 0, 1, or 2"));
        }
        self.strategy = v;
        Ok(())
    }

    #[getter]
    fn get_tol(&self) -> f64 {
        self.tolerance
    }

    #[setter]
    fn set_tol(&mut self, v: f64) {
        self.tolerance = v;
    }

    #[getter]
    fn get_values(slf: Bound<'_, Self>) -> ValueView {
        ValueView {
            owner: slf.unbind(),
        }
    }

    #[setter]
    fn set_values(&mut self, values: HashMap<String, f64>) {
        for (k, v) in values {
            if self.values.contains_key(&k) {
                self.values.insert(k, v);
            }
        }
    }

    #[getter]
    fn get_errors(slf: Bound<'_, Self>) -> ErrorView {
        ErrorView {
            owner: slf.unbind(),
        }
    }

    #[setter]
    fn set_errors(&mut self, errors: HashMap<String, f64>) {
        for (k, v) in errors {
            if self.errors.contains_key(&k) {
                self.errors.insert(k, v);
            }
        }
    }

    #[getter]
    fn get_limits(slf: Bound<'_, Self>) -> LimitView {
        LimitView {
            owner: slf.unbind(),
        }
    }

    #[setter]
    fn set_limits(&mut self, limits: Bound<'_, PyDict>) -> PyResult<()> {
        for (key, value) in limits.iter() {
            let name = key.extract::<String>()?;
            if !self.names.contains(&name) {
                continue;
            }
            if value.is_none() {
                self.limits.remove(&name);
            } else if let Ok(tuple) = value.cast::<PyTuple>()
                && tuple.len() == 2
            {
                let l = finite_or_none(tuple.get_item(0)?.extract::<Option<f64>>()?);
                let u = finite_or_none(tuple.get_item(1)?.extract::<Option<f64>>()?);
                if l.is_none() && u.is_none() {
                    self.limits.remove(&name);
                } else {
                    self.limits.insert(name, (l, u));
                }
            }
        }
        Ok(())
    }

    #[getter]
    fn get_fixed(slf: Bound<'_, Self>) -> FixedView {
        FixedView {
            owner: slf.unbind(),
        }
    }

    #[setter]
    fn set_fixed(&mut self, fixed: Vec<String>) {
        self.fixed.clear();
        for name in fixed {
            if self.names.contains(&name) {
                self.fixed.insert(name);
            }
        }
    }

    #[getter]
    fn get_fval(&self) -> Option<f64> {
        self.scan_fmin
            .as_ref()
            .map(|m| m.fval)
            .or_else(|| self.last_minimum.as_ref().map(|m| m.fval()))
    }

    #[getter]
    fn get_valid(&self) -> Option<bool> {
        self.scan_fmin
            .as_ref()
            .map(|m| m.is_valid)
            .or_else(|| self.last_minimum.as_ref().map(|m| m.is_valid()))
    }

    #[getter]
    fn get_covariance(&self) -> Option<Vec<Vec<f64>>> {
        if self.scan_fmin.is_none()
            && let Some(min) = &self.last_minimum
            && let Some(cov) = min.user_state().covariance()
        {
            let n = cov.nrow();
            let mut matrix = Vec::with_capacity(n);
            for r in 0..n {
                let mut row = Vec::with_capacity(n);
                for c in 0..n {
                    row.push(cov.get(r, c));
                }
                matrix.push(row);
            }
            return Some(matrix);
        }
        None
    }

    #[getter]
    fn get_global_cc(&self) -> Option<Vec<f64>> {
        self.last_minimum
            .as_ref()
            .and_then(|m| m.user_state().global_cc())
            .map(|s| s.to_vec())
    }

    #[getter]
    fn get_fmin(&self) -> Option<FMin> {
        self.scan_fmin
            .clone()
            .or_else(|| self.last_minimum.as_ref().map(FMin::from_minimum))
    }

    #[getter]
    fn get_params(&self) -> Vec<Param> {
        self.names
            .iter()
            .enumerate()
            .map(|(i, name)| {
                let value = *self.values.get(name).unwrap_or(&0.0);
                let error = *self.errors.get(name).unwrap_or(&0.1);
                let is_fixed = self.fixed.contains(name);
                let (lower_limit, upper_limit) =
                    self.limits.get(name).copied().unwrap_or((None, None));
                let has_limits = lower_limit.is_some() || upper_limit.is_some();
                Param {
                    number: i,
                    name: name.clone(),
                    value,
                    error,
                    is_fixed,
                    lower_limit,
                    upper_limit,
                    has_limits,
                }
            })
            .collect()
    }

    #[getter]
    fn get_merrors(&self, py: Python<'_>) -> PyResult<Py<PyDict>> {
        let dict = PyDict::new(py);
        for (k, v) in &self.merrors {
            dict.set_item(k, Py::new(py, MError::from_stored(v))?)?;
        }
        Ok(dict.into())
    }

    #[getter]
    fn get_nfcn(&self) -> usize {
        self.scan_fmin
            .as_ref()
            .map(|m| m.nfcn)
            .or_else(|| self.last_minimum.as_ref().map(|m| m.nfcn()))
            .unwrap_or(0)
    }

    #[getter]
    fn get_npar(&self) -> usize {
        self.names.len()
    }

    #[getter]
    fn get_nfit(&self) -> usize {
        self.names
            .iter()
            .filter(|n| !self.fixed.contains(*n))
            .count()
    }

    #[getter]
    fn get_parameters(&self) -> Vec<String> {
        self.names.clone()
    }

    #[getter]
    fn get_accurate(&self) -> Option<bool> {
        self.last_minimum
            .as_ref()
            .map(|m| m.is_valid() && !m.has_made_pos_def_covar())
    }

    fn reset(mut slf: PyRefMut<'_, Self>) -> Py<Minuit> {
        slf.values = slf.init_values.clone();
        slf.errors = slf.init_errors.clone();
        slf.fixed = slf.init_fixed.clone();
        slf.limits = slf.init_limits.clone();
        slf.last_minimum = None;
        slf.scan_fmin = None;
        slf.merrors.clear();
        slf.into()
    }

    #[pyo3(signature = (key, value))]
    fn fixto(
        mut slf: PyRefMut<'_, Self>,
        key: Bound<'_, PyAny>,
        value: f64,
    ) -> PyResult<Py<Minuit>> {
        let name: String = if let Ok(idx) = key.extract::<usize>() {
            slf.names
                .get(idx)
                .ok_or_else(|| {
                    PyValueError::new_err(format!("parameter index {} out of range", idx))
                })?
                .clone()
        } else {
            let s = key.extract::<String>()?;
            if !slf.names.contains(&s) {
                return Err(PyValueError::new_err(format!("unknown parameter: {}", s)));
            }
            s
        };
        slf.values.insert(name.clone(), value);
        slf.fixed.insert(name);
        Ok(slf.into())
    }

    fn migrad(mut slf: PyRefMut<'_, Self>, py: Python<'_>) -> PyResult<Py<Minuit>> {
        {
            let fcn = PythonFCN {
                fcn: slf.fcn.clone_ref(py),
                errordef: slf.errordef,
            };
            slf.merrors.clear();
            slf.scan_fmin = None;
            let minimizer = slf.build_migrad();
            let result = minimizer.minimize(&fcn);
            slf.update_state_from_result(&result);
            slf.last_minimum = Some(result);
        }
        Ok(slf.into())
    }

    fn simplex(mut slf: PyRefMut<'_, Self>, py: Python<'_>) -> PyResult<Py<Minuit>> {
        {
            let fcn = PythonFCN {
                fcn: slf.fcn.clone_ref(py),
                errordef: slf.errordef,
            };
            slf.merrors.clear();
            slf.scan_fmin = None;
            let minimizer = slf.build_simplex();
            let result = minimizer.minimize(&fcn);
            slf.update_state_from_result(&result);
            slf.last_minimum = Some(result);
        }
        Ok(slf.into())
    }

    fn hesse(mut slf: PyRefMut<'_, Self>, py: Python<'_>) -> PyResult<Py<Minuit>> {
        if slf.last_minimum.is_none() {
            return Err(PyRuntimeError::new_err("Run migrad/simplex first"));
        }

        {
            let min = slf.last_minimum.as_ref().unwrap().clone();
            let fcn = PythonFCN {
                fcn: slf.fcn.clone_ref(py),
                errordef: slf.errordef,
            };
            let mut hesse = MnHesse::new().with_strategy(slf.strategy);
            if let Some(max) = slf.max_calls {
                hesse = hesse.with_max_calls(max);
            }
            // Note: unlike migrad/simplex, hesse refines the covariance without
            // invalidating existing Minos errors, so merrors is preserved here
            // (matches iminuit behavior).
            let result = hesse.calculate(&fcn, &min);
            slf.update_state_from_result(&result);
            slf.last_minimum = Some(result);
        }
        Ok(slf.into())
    }

    #[pyo3(signature = (*parameters))]
    fn minos(
        mut slf: PyRefMut<'_, Self>,
        py: Python<'_>,
        parameters: Vec<String>,
    ) -> PyResult<Py<Minuit>> {
        if slf.last_minimum.is_none() {
            return Err(PyRuntimeError::new_err("Run migrad/simplex first"));
        }

        {
            let min = slf.last_minimum.as_ref().unwrap().clone();
            let selected: Vec<(usize, String)> = if parameters.is_empty() {
                slf.names
                    .iter()
                    .enumerate()
                    .filter(|(_, name)| !slf.fixed.contains(*name))
                    .map(|(i, name)| (i, name.clone()))
                    .collect()
            } else {
                parameters
                    .iter()
                    .filter_map(|name| {
                        if slf.fixed.contains(name) {
                            None
                        } else {
                            slf.names
                                .iter()
                                .position(|n| n == name)
                                .map(|idx| (idx, name.clone()))
                        }
                    })
                    .collect()
            };
            let fcn = PythonFCN {
                fcn: slf.fcn.clone_ref(py),
                errordef: slf.errordef,
            };
            let minos = MnMinos::new(&fcn, &min);
            let mut results = Vec::with_capacity(selected.len());
            for (i, name) in selected {
                let err = minos.minos_error(i);
                results.push(StoredMError {
                    name,
                    // lower_error() is already negative (matches iminuit's
                    // signed-lower convention); store it as-is.
                    lower: err.lower_error(),
                    upper: err.upper_error(),
                    is_valid: err.is_valid(),
                    lower_valid: err.lower_valid(),
                    upper_valid: err.upper_valid(),
                    at_lower_limit: err.at_lower_limit(),
                    at_upper_limit: err.at_upper_limit(),
                    at_lower_max_fcn: err.at_lower_max_fcn(),
                    at_upper_max_fcn: err.at_upper_max_fcn(),
                    lower_new_min: err.lower_new_min(),
                    upper_new_min: err.upper_new_min(),
                    nfcn: err.nfcn(),
                    min: err.min(),
                });
            }
            for err in results {
                slf.merrors.insert(err.name.clone(), err);
            }
        }
        Ok(slf.into())
    }

    #[pyo3(signature = (vname, *, size=100, bound=2.0, subtract_min=false))]
    fn profile<'py>(
        &self,
        py: Python<'py>,
        vname: String,
        size: usize,
        bound: f64,
        subtract_min: bool,
    ) -> PyResult<(Arr1f<'py>, Arr1f<'py>)> {
        let _min = self
            .last_minimum
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("Run migrad/simplex first"))?;
        let idx = self
            .names
            .iter()
            .position(|n| *n == vname)
            .ok_or_else(|| PyValueError::new_err(format!("Parameter '{}' not found", vname)))?;
        let v = *self.values.get(&vname).unwrap_or(&0.0);
        let e = self
            .errors
            .get(&vname)
            .copied()
            .unwrap_or(0.1)
            .abs()
            .max(1e-10);
        let low = v - bound * e;
        let high = v + bound * e;
        let fcn = PythonFCN {
            fcn: self.fcn.clone_ref(py),
            errordef: self.errordef,
        };
        // Build the grid directly (others held at their fitted values) so the
        // result has exactly `size` points. MnParameterScan clamps nsteps to
        // [2, 101], so it cannot honor an arbitrary `size`.
        let base_pars: Vec<f64> = self
            .names
            .iter()
            .map(|n| *self.values.get(n).unwrap_or(&0.0))
            .collect();
        let mut xs = Vec::with_capacity(size);
        let mut fs = Vec::with_capacity(size);
        for i in 0..size {
            let x = if size > 1 {
                low + i as f64 * (high - low) / (size - 1) as f64
            } else {
                v
            };
            let mut pars = base_pars.clone();
            pars[idx] = x;
            xs.push(x);
            fs.push(fcn.value(&pars));
        }
        if subtract_min {
            let fmin_val = fs.iter().copied().fold(f64::INFINITY, f64::min);
            for f in &mut fs {
                *f -= fmin_val;
            }
        }
        Ok((PyArray1::from_vec(py, xs), PyArray1::from_vec(py, fs)))
    }

    #[pyo3(signature = (vname, *, size=30, bound=2.0, subtract_min=false))]
    fn mnprofile<'py>(
        &self,
        py: Python<'py>,
        vname: String,
        size: usize,
        bound: f64,
        subtract_min: bool,
    ) -> PyResult<(Arr1f<'py>, Arr1f<'py>, Arr1b<'py>)> {
        let _min = self
            .last_minimum
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("Run migrad/simplex first"))?;
        self.names
            .iter()
            .position(|n| *n == vname)
            .ok_or_else(|| PyValueError::new_err(format!("Parameter '{}' not found", vname)))?;
        let v = *self.values.get(&vname).unwrap_or(&0.0);
        let e = self
            .errors
            .get(&vname)
            .copied()
            .unwrap_or(0.1)
            .abs()
            .max(1e-10);
        let low = v - bound * e;
        let high = v + bound * e;
        let step = if size > 1 {
            (high - low) / (size - 1) as f64
        } else {
            0.0
        };
        let mut xs = Vec::with_capacity(size);
        let mut fs = Vec::with_capacity(size);
        let mut oks = Vec::with_capacity(size);
        for i in 0..size {
            let xi = if size > 1 { low + i as f64 * step } else { v };
            xs.push(xi);
            let fcn = PythonFCN {
                fcn: self.fcn.clone_ref(py),
                errordef: self.errordef,
            };
            let result = self.build_migrad_with_const(&vname, xi).minimize(&fcn);
            fs.push(result.fval());
            oks.push(result.is_valid());
        }
        if subtract_min {
            let fmin_val = fs.iter().copied().fold(f64::INFINITY, f64::min);
            for f in &mut fs {
                *f -= fmin_val;
            }
        }
        Ok((
            PyArray1::from_vec(py, xs),
            PyArray1::from_vec(py, fs),
            PyArray1::from_vec(py, oks),
        ))
    }

    #[pyo3(signature = (x, y, *, cl=None, size=100))]
    fn mncontour<'py>(
        &self,
        py: Python<'py>,
        x: String,
        y: String,
        cl: Option<f64>,
        size: usize,
    ) -> PyResult<Arr2f<'py>> {
        if cl.is_some() {
            return Err(PyNotImplementedError::new_err(
                "cl scaling not yet implemented; mncontour uses the errordef-defined level",
            ));
        }
        let min = self
            .last_minimum
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("Run migrad/simplex first"))?;
        let idx_x = self
            .names
            .iter()
            .position(|n| *n == x)
            .ok_or_else(|| PyValueError::new_err(format!("Parameter '{}' not found", x)))?;
        let idx_y = self
            .names
            .iter()
            .position(|n| *n == y)
            .ok_or_else(|| PyValueError::new_err(format!("Parameter '{}' not found", y)))?;
        let fcn = PythonFCN {
            fcn: self.fcn.clone_ref(py),
            errordef: self.errordef,
        };
        let contours = MnContours::new(&fcn, min).with_strategy(self.strategy);
        let mut pts = contours.points(idx_x, idx_y, size);
        if let Some(first) = pts.first().copied() {
            pts.push(first);
        }
        let nrows = pts.len();
        let mut arr = Array2::<f64>::zeros((nrows, 2));
        for (i, (px, py_val)) in pts.iter().enumerate() {
            arr[[i, 0]] = *px;
            arr[[i, 1]] = *py_val;
        }
        Ok(arr.into_pyarray(py))
    }

    #[pyo3(signature = (x, y, *, size=50, bound=2.0, subtract_min=false))]
    fn contour<'py>(
        &self,
        py: Python<'py>,
        x: String,
        y: String,
        size: usize,
        bound: f64,
        subtract_min: bool,
    ) -> PyResult<(Arr1f<'py>, Arr1f<'py>, Arr2f<'py>)> {
        let min = self
            .last_minimum
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("Run migrad/simplex first"))?;
        let idx_x = self
            .names
            .iter()
            .position(|n| *n == x)
            .ok_or_else(|| PyValueError::new_err(format!("Parameter '{}' not found", x)))?;
        let idx_y = self
            .names
            .iter()
            .position(|n| *n == y)
            .ok_or_else(|| PyValueError::new_err(format!("Parameter '{}' not found", y)))?;
        let vx = *self.values.get(&x).unwrap_or(&0.0);
        let ex = self.errors.get(&x).copied().unwrap_or(0.1).abs().max(1e-10);
        let vy = *self.values.get(&y).unwrap_or(&0.0);
        let ey = self.errors.get(&y).copied().unwrap_or(0.1).abs().max(1e-10);
        let x_low = vx - bound * ex;
        let x_high = vx + bound * ex;
        let y_low = vy - bound * ey;
        let y_high = vy + bound * ey;
        let xg: Vec<f64> = (0..size)
            .map(|i| {
                if size > 1 {
                    x_low + i as f64 * (x_high - x_low) / (size - 1) as f64
                } else {
                    vx
                }
            })
            .collect();
        let yg: Vec<f64> = (0..size)
            .map(|j| {
                if size > 1 {
                    y_low + j as f64 * (y_high - y_low) / (size - 1) as f64
                } else {
                    vy
                }
            })
            .collect();
        let fmin_val = if subtract_min { min.fval() } else { 0.0 };
        let fcn = PythonFCN {
            fcn: self.fcn.clone_ref(py),
            errordef: self.errordef,
        };
        let mut fval2d = Array2::<f64>::zeros((size, size));
        let base_pars: Vec<f64> = self
            .names
            .iter()
            .map(|n| *self.values.get(n).unwrap_or(&0.0))
            .collect();
        for (i, xi) in xg.iter().enumerate() {
            for (j, yj) in yg.iter().enumerate() {
                let mut pars = base_pars.clone();
                pars[idx_x] = *xi;
                pars[idx_y] = *yj;
                fval2d[[i, j]] = fcn.value(&pars) - fmin_val;
            }
        }
        Ok((
            PyArray1::from_vec(py, xg),
            PyArray1::from_vec(py, yg),
            fval2d.into_pyarray(py),
        ))
    }

    #[pyo3(signature = (ncall=None))]
    fn scan(
        mut slf: PyRefMut<'_, Self>,
        py: Python<'_>,
        ncall: Option<usize>,
    ) -> PyResult<Py<Minuit>> {
        let nfit = slf.names.iter().filter(|n| !slf.fixed.contains(*n)).count();
        let ncall = ncall
            .or(slf.max_calls)
            // iminuit's _migrad_maxcall heuristic: 200 + 100*nfit + 5*nfit^2
            .unwrap_or(200 + 100 * nfit + 5 * nfit * nfit);
        let nstep = if nfit == 0 {
            1
        } else {
            ((ncall as f64).powf(1.0 / nfit as f64) as usize).max(1)
        };
        let fcn = PythonFCN {
            fcn: slf.fcn.clone_ref(py),
            errordef: slf.errordef,
        };
        let mut grids = Vec::with_capacity(slf.names.len());
        for name in &slf.names {
            let v = *slf.values.get(name).unwrap_or(&0.0);
            if slf.fixed.contains(name) {
                grids.push(vec![v]);
                continue;
            }
            let e = slf.errors.get(name).copied().unwrap_or(0.1).abs();
            let (lo, hi) = slf.limits.get(name).copied().unwrap_or((None, None));
            let low = lo.unwrap_or(v - e);
            let high = hi.unwrap_or(v + e);
            let grid = if low == high || nstep == 1 {
                vec![low]
            } else {
                // iminuit uses nstep = int(ncall ** (1 / nfit)) and evaluates
                // np.linspace(low, high, nstep), including both boundaries.
                (0..nstep)
                    .map(|i| low + i as f64 * (high - low) / (nstep - 1) as f64)
                    .collect()
            };
            grids.push(grid);
        }
        let mut point = vec![0.0; slf.names.len()];
        let mut best = point.clone();
        let mut best_f = f64::INFINITY;
        let mut nfcn = 0usize;
        fn run(
            ipar: usize,
            grids: &[Vec<f64>],
            point: &mut [f64],
            best: &mut Vec<f64>,
            best_f: &mut f64,
            nfcn: &mut usize,
            fcn: &PythonFCN,
        ) {
            if ipar == grids.len() {
                let f = fcn.value(point);
                *nfcn += 1;
                if f < *best_f {
                    *best_f = f;
                    best.clone_from_slice(point);
                }
            } else {
                for &x in &grids[ipar] {
                    point[ipar] = x;
                    run(ipar + 1, grids, point, best, best_f, nfcn, fcn);
                }
            }
        }
        run(
            0,
            &grids,
            &mut point,
            &mut best,
            &mut best_f,
            &mut nfcn,
            &fcn,
        );
        for (name, value) in slf.names.clone().into_iter().zip(best.iter().copied()) {
            slf.values.insert(name, value);
        }
        let mut edm = 0.0;
        for (i, name) in slf.names.iter().enumerate() {
            if slf.fixed.contains(name) {
                continue;
            }
            let step = slf.errors.get(name).copied().unwrap_or(0.1).abs().max(1.0) * 1e-4;
            let mut lo = best.clone();
            let mut hi = best.clone();
            lo[i] -= step;
            hi[i] += step;
            if let Some((lower, upper)) = slf.limits.get(name).copied() {
                if let Some(lower) = lower {
                    lo[i] = lo[i].max(lower);
                }
                if let Some(upper) = upper {
                    hi[i] = hi[i].min(upper);
                }
            }
            if hi[i] == lo[i] {
                continue;
            }
            let flo = fcn.value(&lo);
            let fhi = fcn.value(&hi);
            nfcn += 2;
            let grad = (fhi - flo) / (hi[i] - lo[i]);
            let hess = (fhi - 2.0 * best_f + flo) / (0.5 * (hi[i] - lo[i])).powi(2);
            if hess.is_finite() && hess > 0.0 {
                edm += 0.5 * grad * grad / hess;
            } else {
                edm = f64::INFINITY;
                break;
            }
        }
        let edm_goal = slf.tolerance * slf.errordef * 0.002;
        let is_valid = edm <= edm_goal;
        slf.last_minimum = None;
        slf.scan_fmin = Some(FMin {
            fval: best_f,
            edm,
            nfcn,
            errordef: slf.errordef,
            is_valid,
            has_valid_parameters: true,
            has_covariance: false,
            has_made_posdef_covar: false,
            is_above_max_edm: !is_valid,
            has_reached_call_limit: false,
        });
        slf.merrors.clear();
        Ok(slf.into())
    }
}

impl Minuit {
    fn build_migrad(&self) -> MnMigrad {
        let mut m = MnMigrad::new()
            .with_strategy(self.strategy)
            .tolerance(self.tolerance);
        if let Some(max) = self.max_calls {
            m = m.max_fcn(max);
        }
        for name in &self.names {
            let val = *self.values.get(name).unwrap_or(&0.0);
            let err = *self.errors.get(name).unwrap_or(&0.1);
            if self.fixed.contains(name) {
                m = m.add_const(name, val);
            } else {
                m = match self.limits.get(name) {
                    Some((Some(l), Some(u))) => m.add_limited(name, val, err, *l, *u),
                    Some((Some(l), None)) => m.add_lower_limited(name, val, err, *l),
                    Some((None, Some(u))) => m.add_upper_limited(name, val, err, *u),
                    _ => m.add(name, val, err),
                };
            }
        }
        m
    }

    fn build_simplex(&self) -> MnSimplex {
        let mut m = MnSimplex::new()
            .with_strategy(self.strategy)
            .tolerance(self.tolerance);
        if let Some(max) = self.max_calls {
            m = m.max_fcn(max);
        }
        for name in &self.names {
            let val = *self.values.get(name).unwrap_or(&0.0);
            let err = *self.errors.get(name).unwrap_or(&0.1);
            if self.fixed.contains(name) {
                m = m.add_const(name, val);
            } else {
                m = match self.limits.get(name) {
                    Some((Some(l), Some(u))) => m.add_limited(name, val, err, *l, *u),
                    Some((Some(l), None)) => m.add_lower_limited(name, val, err, *l),
                    Some((None, Some(u))) => m.add_upper_limited(name, val, err, *u),
                    _ => m.add(name, val, err),
                };
            }
        }
        m
    }

    fn build_migrad_with_const(&self, fixed_name: &str, fixed_val: f64) -> MnMigrad {
        let mut m = MnMigrad::new()
            .with_strategy(self.strategy)
            .tolerance(self.tolerance);
        if let Some(max) = self.max_calls {
            m = m.max_fcn(max);
        }
        for name in &self.names {
            if name == fixed_name {
                m = m.add_const(name, fixed_val);
            } else {
                let val = *self.values.get(name).unwrap_or(&0.0);
                let err = *self.errors.get(name).unwrap_or(&0.1);
                if self.fixed.contains(name) {
                    m = m.add_const(name, val);
                } else {
                    m = match self.limits.get(name) {
                        Some((Some(l), Some(u))) => m.add_limited(name, val, err, *l, *u),
                        Some((Some(l), None)) => m.add_lower_limited(name, val, err, *l),
                        Some((None, Some(u))) => m.add_upper_limited(name, val, err, *u),
                        _ => m.add(name, val, err),
                    };
                }
            }
        }
        m
    }

    fn update_state_from_result(&mut self, result: &FunctionMinimum) {
        let user_state = result.user_state();
        let params = user_state.params();
        for i in 0..params.len() {
            let p = params.trafo().parameter(i);
            if let Some(name) = self.names.get(i) {
                self.values.insert(name.clone(), p.value());
                self.errors.insert(name.clone(), p.error());
            }
        }
    }
}

// ============================================================================
// PyModule: _minuit2
// ============================================================================

#[pymodule]
fn _minuit2(m: &Bound<'_, pyo3::types::PyModule>) -> PyResult<()> {
    m.add_class::<Minuit>()?;
    m.add_class::<FMin>()?;
    m.add_class::<Param>()?;
    m.add_class::<MError>()?;
    m.add_class::<ValueView>()?;
    m.add_class::<ErrorView>()?;
    m.add_class::<FixedView>()?;
    m.add_class::<LimitView>()?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
