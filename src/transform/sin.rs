use std::f64::consts::FRAC_PI_2;

use super::ParameterTransform;
use crate::precision::MnMachinePrecision;

/// Sine transform for doubly-bounded parameters.
///
/// Maps a bounded external parameter `x in [lower, upper]` to an internal
/// angle `theta` with
/// `x = lower + (upper - lower) * (sin(theta) + 1) / 2`.
/// Near exact bounds, the inverse is clamped away from `+/- pi/2` so the
/// derivative remains numerically usable.
pub struct SinTransform;

impl SinTransform {
    pub fn dext2int(&self, value: f64, upper: f64, lower: f64, prec: &MnMachinePrecision) -> f64 {
        let int = self.ext2int(value, upper, lower, prec);
        let d = self.dint2ext(int, upper, lower);
        if d.abs() > prec.eps2() { 1.0 / d } else { 0.0 }
    }
}

impl ParameterTransform for SinTransform {
    fn int2ext(&self, angle: f64, upper: f64, lower: f64) -> f64 {
        let width = upper - lower;
        lower + 0.5 * width * (angle.sin() + 1.0)
    }

    fn ext2int(&self, value: f64, upper: f64, lower: f64, prec: &MnMachinePrecision) -> f64 {
        let width = upper - lower;
        let normalized = 2.0 * (value - lower) / width - 1.0;

        let boundary_margin = 8.0 * prec.eps2().sqrt();
        let max_angle = FRAC_PI_2 - boundary_margin;
        let min_angle = -FRAC_PI_2 + boundary_margin;

        match normalized {
            n if n <= -1.0 + boundary_margin => min_angle,
            n if n >= 1.0 - boundary_margin => max_angle,
            n => n.asin(),
        }
    }

    fn dint2ext(&self, angle: f64, upper: f64, lower: f64) -> f64 {
        0.5 * (upper - lower).abs() * angle.cos().abs()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() {
        let t = SinTransform;
        let prec = MnMachinePrecision::new();
        let (lo, hi) = (1.0, 10.0);

        for &ext in &[2.0, 5.5, 9.0] {
            let int = t.ext2int(ext, hi, lo, &prec);
            let back = t.int2ext(int, hi, lo);
            assert!((back - ext).abs() < 1e-12, "roundtrip failed for {ext}");
        }
    }

    #[test]
    fn midpoint() {
        let t = SinTransform;
        // internal = 0 should map to midpoint
        let ext = t.int2ext(0.0, 10.0, 0.0);
        assert!((ext - 5.0).abs() < 1e-15);
    }

    #[test]
    fn derivative_positive() {
        let t = SinTransform;
        let d = t.dint2ext(0.0, 10.0, 0.0);
        assert!(d > 0.0);
    }

    #[test]
    fn near_boundary() {
        let t = SinTransform;
        let prec = MnMachinePrecision::new();
        // Very close to upper bound
        let int = t.ext2int(9.9999999999, 10.0, 0.0, &prec);
        let ext = t.int2ext(int, 10.0, 0.0);
        assert!((ext - 10.0).abs() < 0.01);
    }

    #[test]
    fn exact_bounds_are_clamped_to_finite_internal_values() {
        let t = SinTransform;
        let prec = MnMachinePrecision::new();

        let lower = t.ext2int(0.0, 10.0, 0.0, &prec);
        let upper = t.ext2int(10.0, 10.0, 0.0, &prec);

        assert!(lower > -FRAC_PI_2);
        assert!(upper < FRAC_PI_2);
        assert!(lower.is_finite());
        assert!(upper.is_finite());
    }

    #[test]
    fn dext2int_matches_reciprocal_midpoint_derivative() {
        let t = SinTransform;
        let prec = MnMachinePrecision::new();

        let reciprocal = t.dext2int(5.0, 10.0, 0.0, &prec);

        assert!((reciprocal - 0.2).abs() < 1e-15);
    }
}
