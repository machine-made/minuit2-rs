use super::ParameterTransform;
use crate::precision::MnMachinePrecision;

/// Square-root transform for lower-bounded parameters.
///
/// Maps `x in [lower, +inf)` to an unconstrained internal coordinate with
/// `x = lower - 1 + sqrt(y^2 + 1)`.
/// # Examples
///
/// ```
/// use minuit2::transform::{SqrtLowTransform, ParameterTransform};
/// use minuit2::MnMachinePrecision;
/// let t = SqrtLowTransform;
/// let internal = t.ext2int(2.0, 0.0, 0.0, &MnMachinePrecision::new());
/// assert!((t.int2ext(internal, 0.0, 0.0) - 2.0).abs() < 1e-12);
/// ```
pub struct SqrtLowTransform;

impl SqrtLowTransform {
    /// Return the derivative of the internal coordinate with respect to the external value.
    #[must_use]
    pub fn dext2int(&self, value: f64, upper: f64, lower: f64, prec: &MnMachinePrecision) -> f64 {
        let int = self.ext2int(value, upper, lower, prec);
        let d = self.dint2ext(int, upper, lower);
        if d.abs() > prec.eps2() { 1.0 / d } else { 0.0 }
    }
}

impl ParameterTransform for SqrtLowTransform {
    fn int2ext(&self, value: f64, _upper: f64, lower: f64) -> f64 {
        lower - 1.0 + (value * value + 1.0).sqrt()
    }

    fn ext2int(&self, value: f64, _upper: f64, lower: f64, prec: &MnMachinePrecision) -> f64 {
        let yy = value - lower + 1.0;
        let yy2 = yy * yy - 1.0;
        if yy2 < prec.eps2() {
            // Too close to the bound — return 0
            0.0
        } else {
            yy2.sqrt()
        }
    }

    fn dint2ext(&self, value: f64, _upper: f64, _lower: f64) -> f64 {
        value / (value * value + 1.0).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() {
        let t = SqrtLowTransform;
        let prec = MnMachinePrecision::new();
        let lo = 2.0;

        for &ext in &[3.0, 5.0, 100.0] {
            let int = t.ext2int(ext, f64::INFINITY, lo, &prec);
            let back = t.int2ext(int, f64::INFINITY, lo);
            assert!(
                (back - ext).abs() < 1e-10,
                "roundtrip failed for {ext}: got {back}"
            );
        }
    }

    #[test]
    fn at_bound() {
        let t = SqrtLowTransform;
        let prec = MnMachinePrecision::new();
        let int = t.ext2int(2.0, f64::INFINITY, 2.0, &prec);
        assert!((int).abs() < 1e-10, "at bound should be ~0, got {int}");
    }

    #[test]
    fn int_zero_maps_to_bound() {
        let t = SqrtLowTransform;
        let ext = t.int2ext(0.0, f64::INFINITY, 5.0);
        // int=0: lo - 1 + sqrt(0+1) = lo - 1 + 1 = lo
        assert!((ext - 5.0).abs() < 1e-15);
    }
}
