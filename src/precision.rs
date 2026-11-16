/// Machine precision of the numerical steps.
///
/// `eps` is four times `f64::EPSILON`. `eps2` is `2 * sqrt(eps)`.
#[derive(Debug, Clone, Copy)]
pub struct MnMachinePrecision {
    eps: f64,
    eps2: f64,
}

impl MnMachinePrecision {
    /// Return the precision ROOT stores for double arithmetic.
    pub fn new() -> Self {
        // MnMachinePrecision.cxx:26-27. Four times machine epsilon.
        let eps = 4.0 * f64::EPSILON;
        Self {
            eps,
            eps2: 2.0 * eps.sqrt(),
        }
    }

    /// Return the stored machine epsilon.
    pub fn eps(&self) -> f64 {
        self.eps
    }

    /// Return `2 * sqrt(eps)`.
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
        let eps = 4.0 * f64::EPSILON;
        assert!((p.eps() - eps).abs() < 1e-30);
        assert!((p.eps2() - 2.0 * eps.sqrt()).abs() < 1e-20);
    }

    #[test]
    fn custom_precision() {
        let mut p = MnMachinePrecision::new();
        p.set_precision(1e-10);
        assert!((p.eps() - 1e-10).abs() < 1e-25);
        assert!((p.eps2() - 2.0e-5).abs() < 1e-15);
    }
}
