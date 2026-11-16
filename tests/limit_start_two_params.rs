//! Check Migrad when a bounded parameter starts on its limit.

use minuit2::MnMigrad;
use minuit2::minimum::FunctionMinimum;

fn quadratic2(p: &[f64]) -> f64 {
    let (x, y) = (p[0], p[1]);
    (x - 1.0).powi(2) + 4.0 * (y + 2.0).powi(2) + 0.3 * x * y
}

fn assert_closed_form(r: &FunctionMinimum) {
    let f = r.fval();
    let x = r.params()[0];
    let y = r.params()[1];
    // The failure this file guards against ended with valid=true on a limit, so
    // the status is part of the check.
    assert!(
        r.is_valid(),
        "result is not valid: fval {f}, x {x}, y {y}, edm {}, nfcn {}",
        r.edm(),
        r.nfcn()
    );
    // 1e-4 on fval: the cases below end within 1.4e-6 of the closed form and
    // iminuit 2.33.0 (scripts/limit_start_iminuit.py) within 2.9e-5, so 1e-4 is
    // 3 times the iminuit worst case. The failure this file guards against ends
    // at fval 3.6775, 4.4 above the minimum.
    assert!(
        (f - (-0.69956)).abs() < 1e-4,
        "fval {f}, x {x}, y {y}, edm {}, nfcn {}, expected -0.69956 within 1e-4",
        r.edm(),
        r.nfcn()
    );
    // -1.5: the failing runs ended with y on its limit at -1; the minimum is at
    // y = -2.04903, so -1.5 separates the two.
    assert!(y < -1.5, "y {y} stayed on the limit, expected y below -1.5");
    // 1e-2 on x and y: the cases below end within 1.2e-3 in x and 3.8e-4 in y
    // of the closed form, and iminuit 2.33.0 within 2.4e-3 in x and 2.5e-3 in y,
    // so 1e-2 is 4 times the iminuit worst case.
    assert!(
        (y - (-2.04903)).abs() < 1e-2,
        "y {y}, expected -2.04903 within 1e-2"
    );
    assert!(
        (x - 1.3074).abs() < 1e-2,
        "x {x}, expected 1.3074 within 1e-2"
    );
}

#[test]
fn both_limits_y_starts_on_upper() {
    let r = MnMigrad::new()
        .add_limited("x", 0.4, 0.1, 0.0, 2.0)
        .add_limited("y", -1.0, 0.1, -3.0, -1.0)
        .tolerance(0.1)
        .minimize(&quadratic2);
    assert_closed_form(&r);
}

#[test]
fn both_limits_y_starts_on_upper_step_0_5() {
    let r = MnMigrad::new()
        .add_limited("x", 0.4, 0.1, 0.0, 2.0)
        .add_limited("y", -1.0, 0.5, -3.0, -1.0)
        .tolerance(0.1)
        .minimize(&quadratic2);
    assert_closed_form(&r);
}

#[test]
fn free_x_y_starts_on_upper() {
    let r = MnMigrad::new()
        .add("x", 0.4, 0.1)
        .add_limited("y", -1.0, 0.1, -3.0, -1.0)
        .tolerance(0.1)
        .minimize(&quadratic2);
    assert_closed_form(&r);
}

#[test]
fn lower_limited_x_y_starts_on_upper() {
    let r = MnMigrad::new()
        .add_lower_limited("x", 0.4, 0.1, 0.0)
        .add_limited("y", -1.0, 0.1, -3.0, -1.0)
        .tolerance(0.1)
        .minimize(&quadratic2);
    assert_closed_form(&r);
}

#[test]
fn free_x_y_upper_limit_only() {
    let r = MnMigrad::new()
        .add("x", 0.4, 0.1)
        .add_upper_limited("y", -1.0, 0.1, -1.0)
        .tolerance(0.1)
        .minimize(&quadratic2);
    assert_closed_form(&r);
}

#[test]
fn both_on_one_sided_limits() {
    let r = MnMigrad::new()
        .add_lower_limited("x", 0.0, 0.1, 0.0)
        .add_upper_limited("y", -1.0, 0.1, -1.0)
        .tolerance(0.1)
        .minimize(&quadratic2);
    assert_closed_form(&r);
}

#[test]
fn both_on_double_limits() {
    let r = MnMigrad::new()
        .add_limited("x", 0.0, 0.1, 0.0, 2.0)
        .add_limited("y", -1.0, 0.1, -3.0, -1.0)
        .tolerance(0.1)
        .minimize(&quadratic2);
    assert_closed_form(&r);
}
