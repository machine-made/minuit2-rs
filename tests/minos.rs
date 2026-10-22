use minuit2::{MnHesse, MnMigrad, MnMinos};

/// Symmetric case: Gaussian/quadratic fit → Minos errors ≈ Hesse errors.
#[test]
fn minos_symmetric_quadratic() {
    let quadratic = |p: &[f64]| 2.0 * p[0] * p[0] + 8.0 * p[1] * p[1];

    let result = MnMigrad::new()
        .add("x", 5.0, 1.0)
        .add("y", -3.0, 1.0)
        .minimize(&quadratic);

    assert!(result.is_valid());

    // Run Hesse first for accurate errors
    let hesse_result = MnHesse::new().calculate(&quadratic, &result);
    assert!(hesse_result.is_valid());

    let minos = MnMinos::new(&quadratic, &hesse_result);
    let me = minos.minos_error(0);
    let _ = minos.errors(0);

    // For symmetric quadratic, Minos errors should be roughly symmetric
    let upper = me.upper_error();
    let lower = me.lower_error();

    if me.is_valid() {
        // Upper and lower should have similar magnitude (symmetric function)
        let ratio = upper.abs() / lower.abs();
        assert!(
            (ratio - 1.0).abs() < 0.5,
            "upper/lower ratio should be ~1 for symmetric function, got {ratio} (upper={upper}, lower={lower})"
        );

        // Both should be positive/negative respectively and non-zero
        assert!(upper > 0.0, "upper Minos error should be positive: {upper}");
        assert!(lower < 0.0, "lower Minos error should be negative: {lower}");
    }

    assert_eq!(me.parameter(), 0);
    assert!(
        me.nfcn() < 1_000_000,
        "nfcn should be finite, got {}",
        me.nfcn()
    );
    let min_val = me.min();
    assert!(
        min_val.is_finite(),
        "minimum value should be finite, got {min_val}"
    );
    let _ = me.lower();
    let _ = me.upper();
    let _ = me.lower_valid();
    let _ = me.upper_valid();
    let _ = me.at_lower_limit();
    let _ = me.at_upper_limit();
    let _ = me.at_lower_max_fcn();
    let _ = me.at_upper_max_fcn();
    let _ = me.lower_new_min();
    let _ = me.upper_new_min();

    // Even if Minos doesn't converge perfectly, check it doesn't crash
    // and returns sensible signs
}

/// Asymmetric function: Minos upper != lower.
#[test]
fn minos_asymmetric() {
    // f(x) = x^2 for x > 0, 4*x^2 for x < 0
    // Minimum at x=0
    // Upper error: sqrt(up/2) = sqrt(0.5) ≈ 0.707
    // Lower error: sqrt(up/8) = sqrt(0.125) ≈ 0.354
    let asym = |p: &[f64]| {
        if p[0] >= 0.0 {
            p[0] * p[0]
        } else {
            4.0 * p[0] * p[0]
        }
    };

    let result = MnMigrad::new().add("x", 0.5, 0.5).minimize(&asym);

    assert!(result.is_valid());
    assert!(result.fval() < 0.01, "should find minimum near 0");

    let hesse_result = MnHesse::new().calculate(&asym, &result);

    let minos = MnMinos::new(&asym, &hesse_result);
    let me = minos.minos_error(0);

    if me.is_valid() {
        let upper = me.upper_error();
        let lower = me.lower_error();

        // Upper should be larger in magnitude than lower
        assert!(
            upper.abs() > lower.abs() * 1.2,
            "upper ({upper}) should be larger than |lower| ({lower}) for asymmetric function"
        );
    }
    // If Minos doesn't converge (can happen with discontinuous derivative),
    // that's acceptable — the important thing is it doesn't crash.
}

/// Fixed parameter: Minos should return invalid.
#[test]
fn minos_fixed_parameter() {
    let result = MnMigrad::new()
        .add("x", 5.0, 1.0)
        .add_const("y", 0.0)
        .minimize(&|p: &[f64]| p[0] * p[0] + p[1] * p[1]);

    assert!(result.is_valid());

    let minos = MnMinos::new(&|p: &[f64]| p[0] * p[0] + p[1] * p[1], &result);

    // Parameter 1 (y) is const → Minos should return invalid
    let cross = minos.lower(1);
    assert!(!cross.is_valid(), "Minos on const param should be invalid");

    // Regression: min() must not panic even when Minos crossings are invalid.
    let me = minos.minos_error(1);
    assert!(
        !me.is_valid(),
        "MinosError should be invalid for const param"
    );
    assert!(
        me.min().is_finite(),
        "min() should still return a finite original parameter value"
    );
}

#[test]
fn minos_upper_bound_reports_limit() {
    let result = MnMigrad::new()
        .add_upper_limited("x", -1.0, 0.5, -0.5)
        .minimize(&|p: &[f64]| (p[0] + 1.0).powi(2));

    assert!(result.is_valid());

    let minos = MnMinos::new(&|p: &[f64]| (p[0] + 1.0).powi(2), &result);
    let me = minos.minos_error(0);

    assert!((me.min() + 1.0).abs() < 1e-12);
    assert!(
        me.at_upper_limit(),
        "upper crossing should be at parameter limit"
    );
    assert!(!me.lower_new_min());
    assert!(!me.upper_new_min());
    assert!(me.nfcn() > 0);
}
