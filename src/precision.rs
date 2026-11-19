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
    #[must_use]
    pub fn new() -> Self {
        // MnMachinePrecision.cxx:26-27. Four times machine epsilon.
        let eps = 4.0 * f64::EPSILON;
        Self {
            eps,
            eps2: 2.0 * eps.sqrt(),
        }
    }

    /// Return the stored machine epsilon.
    #[must_use]
    pub fn eps(&self) -> f64 {
        self.eps
    }

    /// Return `2 * sqrt(eps)`.
    #[must_use]
    pub fn eps2(&self) -> f64 {
        self.eps2
    }

    /// Override machine epsilon (for testing or non-standard arithmetic).
    pub fn set_precision(&mut self, eps: f64) {
        self.eps = eps;
        self.eps2 = 2.0 * eps.sqrt();
    }

    /// Reset to the value ROOT's `ComputePrecision` finds for double arithmetic.
    ///
    /// The loop of MnMachinePrecision.cxx:29-56 ends at `8 * epstry` with
    /// `epstry = 2^-53`. That is `4 * f64::EPSILON`, the value of `new()`.
    pub fn compute_precision(&mut self) {
        *self = Self::new();
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
        // Both sides are the same expression. 4 * EPSILON = 2^-50 is a power of
        // two, so its square root 2^-25 and the factor 2 are exact and both
        // differences are 0. The bound 1e-30 is about 5 ulp of eps (ulp 2.0e-31).
        assert!((p.eps() - eps).abs() < 1e-30);
        // The bound 1e-20 is about 750 ulp of eps2 (ulp 1.3e-23). It leaves
        // room for a few ulp if `sqrt` or the product rounded differently.
        assert!((p.eps2() - 2.0 * eps.sqrt()).abs() < 1e-20);
    }

    #[test]
    fn custom_precision() {
        let mut p = MnMachinePrecision::new();
        p.set_precision(1e-10);
        assert!((p.eps() - 1e-10).abs() < 1e-25);
        assert!((p.eps2() - 2.0e-5).abs() < 1e-15);
    }

    /// ROOT's `ComputePrecision` ends at `8 * epstry` with `epstry = 2^-53`,
    /// which is `4 * f64::EPSILON` (MnMachinePrecision.cxx:29-56). The result
    /// is the value of the constructor, not `f64::EPSILON`.
    #[test]
    fn compute_precision_gives_root_value() {
        let mut p = MnMachinePrecision::new();
        p.set_precision(1e-10);
        p.compute_precision();
        let eps = 4.0 * f64::EPSILON;
        // Same expression as the constructor, so the values are bit-equal.
        assert_eq!(p.eps(), eps);
        assert_eq!(p.eps2(), 2.0 * eps.sqrt());
    }
}
