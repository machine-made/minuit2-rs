/// A single minimization parameter with name, value, error, and optional bounds.
///
/// Parameters can be free, fixed, or constant.
/// "Constant" means permanently fixed (never released during minimization).
#[derive(Debug, Clone)]
pub struct MinuitParameter {
    num: usize,
    name: String,
    value: f64,
    error: f64,
    is_const: bool,
    is_fixed: bool,
    has_lower_limit: bool,
    has_upper_limit: bool,
    lower_limit: f64,
    upper_limit: f64,
}

impl MinuitParameter {
    /// Free parameter with no bounds.
    pub fn new(num: usize, name: impl Into<String>, value: f64, error: f64) -> Self {
        Self {
            num,
            name: name.into(),
            value,
            error,
            is_const: false,
            is_fixed: false,
            has_lower_limit: false,
            has_upper_limit: false,
            lower_limit: 0.0,
            upper_limit: 0.0,
        }
    }

    /// Parameter with lower bound only.
    pub fn with_lower_limit(
        num: usize,
        name: impl Into<String>,
        value: f64,
        error: f64,
        lower: f64,
    ) -> Self {
        Self {
            num,
            name: name.into(),
            value,
            error,
            is_const: false,
            is_fixed: false,
            has_lower_limit: true,
            has_upper_limit: false,
            lower_limit: lower,
            upper_limit: 0.0,
        }
    }

    /// Parameter with upper bound only.
    pub fn with_upper_limit(
        num: usize,
        name: impl Into<String>,
        value: f64,
        error: f64,
        upper: f64,
    ) -> Self {
        Self {
            num,
            name: name.into(),
            value,
            error,
            is_const: false,
            is_fixed: false,
            has_lower_limit: false,
            has_upper_limit: true,
            lower_limit: 0.0,
            upper_limit: upper,
        }
    }

    /// Parameter with both bounds.
    pub fn with_limits(
        num: usize,
        name: impl Into<String>,
        value: f64,
        error: f64,
        lower: f64,
        upper: f64,
    ) -> Self {
        Self {
            num,
            name: name.into(),
            value,
            error,
            is_const: false,
            is_fixed: false,
            has_lower_limit: true,
            has_upper_limit: true,
            lower_limit: lower,
            upper_limit: upper,
        }
    }

    /// Constant parameter (fixed, never released).
    pub fn constant(num: usize, name: impl Into<String>, value: f64) -> Self {
        Self {
            num,
            name: name.into(),
            value,
            error: 0.0,
            is_const: true,
            is_fixed: true,
            has_lower_limit: false,
            has_upper_limit: false,
            lower_limit: 0.0,
            upper_limit: 0.0,
        }
    }

    pub fn number(&self) -> usize {
        self.num
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn set_name(&mut self, name: impl Into<String>) {
        self.name = name.into();
    }

    pub fn value(&self) -> f64 {
        self.value
    }

    pub fn error(&self) -> f64 {
        self.error
    }

    pub fn set_value(&mut self, val: f64) {
        self.value = val;
    }

    pub fn set_error(&mut self, err: f64) {
        self.error = err;
    }

    // --- Limits ---

    pub fn set_limits(&mut self, lower: f64, upper: f64) {
        assert!(lower < upper, "lower limit must be less than upper limit");
        self.has_lower_limit = true;
        self.has_upper_limit = true;
        self.lower_limit = lower;
        self.upper_limit = upper;
    }

    pub fn set_lower_limit(&mut self, lower: f64) {
        self.has_lower_limit = true;
        self.lower_limit = lower;
    }

    pub fn set_upper_limit(&mut self, upper: f64) {
        self.has_upper_limit = true;
        self.upper_limit = upper;
    }

    pub fn remove_limits(&mut self) {
        self.has_lower_limit = false;
        self.has_upper_limit = false;
    }

    pub fn has_lower_limit(&self) -> bool {
        self.has_lower_limit
    }

    pub fn has_upper_limit(&self) -> bool {
        self.has_upper_limit
    }

    pub fn has_limits(&self) -> bool {
        self.has_lower_limit && self.has_upper_limit
    }

    pub fn lower_limit(&self) -> f64 {
        self.lower_limit
    }

    pub fn upper_limit(&self) -> f64 {
        self.upper_limit
    }

    // --- Fixed/Const ---

    pub fn fix(&mut self) {
        self.is_fixed = true;
    }

    pub fn release(&mut self) {
        if !self.is_const {
            self.is_fixed = false;
        }
    }

    pub fn is_fixed(&self) -> bool {
        self.is_fixed
    }

    pub fn is_const(&self) -> bool {
        self.is_const
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_parameter() {
        let p = MinuitParameter::new(0, "x", 1.0, 0.1);
        assert_eq!(p.number(), 0);
        assert_eq!(p.name(), "x");
        assert!((p.value() - 1.0).abs() < 1e-15);
        assert!(!p.is_fixed());
        assert!(!p.has_limits());
    }

    #[test]
    fn bounded_parameter() {
        let p = MinuitParameter::with_limits(0, "x", 5.0, 0.1, 0.0, 10.0);
        assert!(p.has_limits());
        assert!(p.has_lower_limit());
        assert!(p.has_upper_limit());
        assert!((p.lower_limit() - 0.0).abs() < 1e-15);
        assert!((p.upper_limit() - 10.0).abs() < 1e-15);
    }

    #[test]
    fn fix_release() {
        let mut p = MinuitParameter::new(0, "x", 1.0, 0.1);
        p.fix();
        assert!(p.is_fixed());
        p.release();
        assert!(!p.is_fixed());
    }

    #[test]
    fn const_cannot_release() {
        let mut p = MinuitParameter::constant(0, "x", 42.0);
        assert!(p.is_const());
        assert!(p.is_fixed());
        p.release();
        assert!(p.is_fixed()); // still fixed
    }

    #[test]
    fn set_name_updates_parameter_name() {
        let mut p = MinuitParameter::new(0, "x", 1.0, 0.1);
        p.set_name("alpha");
        assert_eq!(p.name(), "alpha");
    }
}
