use super::ParameterTransform;
use crate::precision::MnMachinePrecision;

/// Square-root transform for upper-bounded parameters.
///
/// Maps `x in (-inf, upper]` to an unconstrained internal coordinate with
/// `x = upper + 1 - sqrt(y^2 + 1)`.
pub struct SqrtUpTransform;

impl SqrtUpTransform {
    #[must_use]
    pub fn dext2int(&self, value: f64, upper: f64, lower: f64, prec: &MnMachinePrecision) -> f64 {
        let int = self.ext2int(value, upper, lower, prec);
        let d = self.dint2ext(int, upper, lower);
        if d.abs() > prec.eps2() { 1.0 / d } else { 0.0 }
    }
}

impl ParameterTransform for SqrtUpTransform {
    fn int2ext(&self, value: f64, upper: f64, _lower: f64) -> f64 {
        // ROOT SqrtUpParameterTransformation.cxx:22-26 evaluates this expression
        // in long double before returning an external double.
        // Keep rounding residuals through the same operations to avoid premature rounding.
        let square = value * value;
        let (radicand, addition_error) = sum_with_error(square, 1.0);
        let root = radicand.sqrt();
        let (offset, offset_error) = sum_with_error(upper, 1.0);
        let (external, external_error) = sum_with_error(offset, -root);
        if !external.is_finite() {
            return external;
        }
        let radicand_error = value.mul_add(value, -square) + addition_error;
        let root_error = ((-root).mul_add(root, radicand) + radicand_error) / (2.0 * root);
        external + (external_error + offset_error - root_error)
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

/// Return a sum and its rounding residual for finite inputs.
fn sum_with_error(left: f64, right: f64) -> (f64, f64) {
    let sum = left + right;
    let recovered_right = sum - left;
    let error = (left - (sum - recovered_right)) + (right - recovered_right);
    (sum, error)
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
