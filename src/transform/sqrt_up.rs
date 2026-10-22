use super::ParameterTransform;
use crate::precision::MnMachinePrecision;

/// Square-root transform for upper-bounded parameters.
///
/// Maps `x in (-inf, upper]` to an unconstrained internal coordinate with
/// `x = upper + 1 - sqrt(y^2 + 1)`.
pub struct SqrtUpTransform;

impl SqrtUpTransform {
    pub fn dext2int(&self, value: f64, upper: f64, lower: f64, prec: &MnMachinePrecision) -> f64 {
        let int = self.ext2int(value, upper, lower, prec);
        let d = self.dint2ext(int, upper, lower);
        if d.abs() > prec.eps2() { 1.0 / d } else { 0.0 }
    }
}

impl ParameterTransform for SqrtUpTransform {
    fn int2ext(&self, value: f64, upper: f64, _lower: f64) -> f64 {
        upper + 1.0 - (value * value + 1.0).sqrt()
    }

    fn ext2int(&self, value: f64, upper: f64, _lower: f64, prec: &MnMachinePrecision) -> f64 {
        let yy = upper - value + 1.0;
        let yy2 = yy * yy - 1.0;
        if yy2 < prec.eps2() { 0.0 } else { yy2.sqrt() }
    }

    fn dint2ext(&self, value: f64, _upper: f64, _lower: f64) -> f64 {
        -value / (value * value + 1.0).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() {
        let t = SqrtUpTransform;
        let prec = MnMachinePrecision::new();
        let hi = 10.0;

        for &ext in &[9.0, 5.0, -100.0] {
            let int = t.ext2int(ext, hi, f64::NEG_INFINITY, &prec);
            let back = t.int2ext(int, hi, f64::NEG_INFINITY);
            assert!(
                (back - ext).abs() < 1e-10,
                "roundtrip failed for {ext}: got {back}"
            );
        }
    }

    #[test]
    fn at_bound() {
        let t = SqrtUpTransform;
        let prec = MnMachinePrecision::new();
        let int = t.ext2int(10.0, 10.0, f64::NEG_INFINITY, &prec);
        assert!((int).abs() < 1e-10, "at bound should be ~0, got {int}");
    }

    #[test]
    fn int_zero_maps_to_bound() {
        let t = SqrtUpTransform;
        let ext = t.int2ext(0.0, 5.0, f64::NEG_INFINITY);
        // int=0: hi + 1 - sqrt(0+1) = hi + 1 - 1 = hi
        assert!((ext - 5.0).abs() < 1e-15);
    }
}
