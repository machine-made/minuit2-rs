//! Property-based metamorphic tests for the Migrad minimizer.
//!
//! Complements the fixed cases in `tests/metamorphic.rs` by sampling MANY
//! randomized objectives/start points with `proptest` and asserting oracle-free
//! invariants: relations between inputs and outputs that must hold regardless of
//! the (unknown) exact numeric answer.
//!
//! Each case is kept cheap (2-3 free parameters) so the suite stays fast.

use minuit2::MnMigrad;
use proptest::prelude::*;

/// Absolute tolerance for recovered parameter values. Migrad's default EDM
/// tolerance is `0.1 * up * 0.002`, so the argmin of a well-conditioned bowl is
/// pinned to roughly 1e-4; we use a slightly looser bound to stay honest across
/// the randomized range.
const PARAM_TOL: f64 = 1e-3;

fn assert_close(got: f64, want: f64, tol: f64, label: &str) -> Result<(), TestCaseError> {
    prop_assert!(
        (got - want).abs() <= tol,
        "{label}: expected {want}, got {got}, diff={}",
        (got - want).abs()
    );
    Ok(())
}

// proptest configuration: enough samples to exercise the space without making
// the suite slow. Each case runs at most a few hundred FCN evaluations.
proptest! {
    #![proptest_config(ProptestConfig {
        cases: 64,
        max_shrink_iters: 256,
        ..ProptestConfig::default()
    })]

    /// Translation invariance: a quadratic bowl `sum w_i * (p_i - t_i)^2` with a
    /// random true minimum `t` (within a bounded range) must be recovered by
    /// Migrad regardless of where `t` sits. The start point is offset from the
    /// truth so the minimizer has to actually move.
    #[test]
    fn migrad_recovers_random_quadratic_minimum(
        t0 in -5.0f64..5.0,
        t1 in -5.0f64..5.0,
        w0 in 0.5f64..4.0,
        w1 in 0.5f64..4.0,
    ) {
        let objective = move |p: &[f64]| {
            w0 * (p[0] - t0).powi(2) + w1 * (p[1] - t1).powi(2)
        };

        let min = MnMigrad::new()
            .add("x", t0 + 2.0, 1.0)
            .add("y", t1 - 2.0, 1.0)
            .minimize(&objective);

        prop_assert!(min.is_valid(), "minimization did not converge");
        let params = min.params();
        assert_close(params[0], t0, PARAM_TOL, "recovered x")?;
        assert_close(params[1], t1, PARAM_TOL, "recovered y")?;
        // At the minimum of a perfect bowl, fval -> 0.
        assert_close(min.fval(), 0.0, PARAM_TOL, "recovered fval")?;
    }

    /// Scaling invariance: multiplying the objective by a positive constant `c`
    /// must not move the recovered minimizer; only `fval` scales by `c`.
    #[test]
    fn migrad_argmin_invariant_under_positive_scaling(
        t0 in -4.0f64..4.0,
        t1 in -4.0f64..4.0,
        w0 in 0.5f64..3.0,
        w1 in 0.5f64..3.0,
        c in 0.25f64..8.0,
    ) {
        let base = move |p: &[f64]| {
            w0 * (p[0] - t0).powi(2) + w1 * (p[1] - t1).powi(2) + 1.0
        };
        let scaled = move |p: &[f64]| c * base(p);

        let base_min = MnMigrad::new()
            .add("x", t0 + 1.5, 1.0)
            .add("y", t1 - 1.5, 1.0)
            .minimize(&base);
        let scaled_min = MnMigrad::new()
            .add("x", t0 + 1.5, 1.0)
            .add("y", t1 - 1.5, 1.0)
            .minimize(&scaled);

        prop_assert!(base_min.is_valid(), "base minimization did not converge");
        prop_assert!(scaled_min.is_valid(), "scaled minimization did not converge");

        let bp = base_min.params();
        let sp = scaled_min.params();
        assert_close(sp[0], bp[0], PARAM_TOL, "scaled argmin x")?;
        assert_close(sp[1], bp[1], PARAM_TOL, "scaled argmin y")?;
        // fval relationship: scaled fval == c * base fval. The bowl bottom is 1,
        // so base fval ~ 1 and scaled fval ~ c; allow tolerance proportional to c.
        assert_close(scaled_min.fval(), c * base_min.fval(), PARAM_TOL * (1.0 + c), "scaled fval")?;
    }

    /// Permutation invariance: reordering the parameters of a separable quadratic
    /// (and the matching builder order) must yield the same recovered values,
    /// permuted the same way. Catches index/name mixups in result mapping.
    #[test]
    fn migrad_commutes_with_parameter_permutation(
        t0 in -4.0f64..4.0,
        t1 in -4.0f64..4.0,
        t2 in -4.0f64..4.0,
        w0 in 0.5f64..3.0,
        w1 in 0.5f64..3.0,
        w2 in 0.5f64..3.0,
    ) {
        // Original order: p[0], p[1], p[2].
        let original = move |p: &[f64]| {
            w0 * (p[0] - t0).powi(2) + w1 * (p[1] - t1).powi(2) + w2 * (p[2] - t2).powi(2)
        };
        // Permuted order (1,2,0): the FCN's p[0]<-old p1, p[1]<-old p2, p[2]<-old p0.
        let permuted = move |p: &[f64]| {
            w1 * (p[0] - t1).powi(2) + w2 * (p[1] - t2).powi(2) + w0 * (p[2] - t0).powi(2)
        };

        let original_min = MnMigrad::new()
            .add("a", t0 - 2.0, 1.0)
            .add("b", t1 + 2.0, 1.0)
            .add("c", t2 - 2.0, 1.0)
            .minimize(&original);
        let permuted_min = MnMigrad::new()
            .add("b", t1 + 2.0, 1.0)
            .add("c", t2 - 2.0, 1.0)
            .add("a", t0 - 2.0, 1.0)
            .minimize(&permuted);

        prop_assert!(original_min.is_valid(), "original minimization did not converge");
        prop_assert!(permuted_min.is_valid(), "permuted minimization did not converge");

        let op = original_min.params();
        let pp = permuted_min.params();
        // Undo the (1,2,0) permutation: original index 0 is permuted slot 2, etc.
        let permuted_back = [pp[2], pp[0], pp[1]];
        assert_close(permuted_back[0], op[0], PARAM_TOL, "permuted-back a")?;
        assert_close(permuted_back[1], op[1], PARAM_TOL, "permuted-back b")?;
        assert_close(permuted_back[2], op[2], PARAM_TOL, "permuted-back c")?;
    }

    /// Starting-point robustness: for a fixed convex quadratic, Migrad must
    /// converge to the same minimum from random valid starting points.
    #[test]
    fn migrad_converges_to_same_minimum_from_random_starts(
        s0 in -8.0f64..8.0,
        s1 in -8.0f64..8.0,
    ) {
        // Fixed convex bowl with a known minimizer at (2.0, -1.0).
        let objective = |p: &[f64]| {
            2.0 * (p[0] - 2.0).powi(2) + 0.75 * (p[1] + 1.0).powi(2)
        };

        let min = MnMigrad::new()
            .add("x", s0, 1.0)
            .add("y", s1, 1.0)
            .minimize(&objective);

        prop_assert!(min.is_valid(), "minimization did not converge from random start");
        let params = min.params();
        assert_close(params[0], 2.0, PARAM_TOL, "x from random start")?;
        assert_close(params[1], -1.0, PARAM_TOL, "y from random start")?;
    }
}
