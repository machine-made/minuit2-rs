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
    assert!(
        (f - (-0.69956)).abs() < 1e-4,
        "fval {f}, x {x}, y {y}, edm {}, nfcn {}, expected -0.69956 within 1e-4",
        r.edm(),
        r.nfcn()
    );
    assert!(y < -1.5, "y {y} stayed on the limit, expected y below -1.5");
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
