/// Machine precision constants for numerical differentiation.
///
/// Uses `f64::EPSILON` (2^-52 ≈ 2.22e-16) as the machine epsilon basis.
#[derive(Debug, Clone, Copy)]
pub struct MnMachinePrecision {
    eps: f64,
    eps2: f64,
}

impl MnMachinePrecision {
    /// Create a new precision object with default f64 epsilon.
    pub fn new() -> Self {
        let eps = f64::EPSILON;
        Self {
            eps,
            eps2: 2.0 * eps.sqrt(),
        }
    }

    /// Machine epsilon (~2.22e-16 for f64).
    pub fn eps(&self) -> f64 {
        self.eps
    }

    /// 2 * sqrt(eps), used as default step size for gradient calculations.
    pub fn eps2(&self) -> f64 {
        self.eps2
    }

    /// Override machine epsilon (for testing or non-standard arithmetic).
    pub fn set_precision(&mut self, eps: f64) {
        self.eps = eps;
        self.eps2 = 2.0 * eps.sqrt();
    }

    pub fn compute_precision(&mut self) {
        self.set_precision(f64::EPSILON);
    }
}

impl Default for MnMachinePrecision {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_precision() {
        let p = MnMachinePrecision::new();
        assert!((p.eps() - f64::EPSILON).abs() < 1e-30);
        assert!((p.eps2() - 2.0 * f64::EPSILON.sqrt()).abs() < 1e-20);
    }

    #[test]
    fn custom_precision() {
        let mut p = MnMachinePrecision::new();
        p.set_precision(1e-10);
        assert!((p.eps() - 1e-10).abs() < 1e-25);
        assert!((p.eps2() - 2.0e-5).abs() < 1e-15);
    }
}
