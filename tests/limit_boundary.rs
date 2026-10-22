//! Regression: a free parameter started EXACTLY on a bound must still escape
//! toward an interior minimum (matches ROOT/iminuit). Before the negative-g2
//! seed escape, the symmetric bound transform gave a zero internal gradient,
//! so migrad reported convergence (EDM=0) without moving off the bound.

use minuit2::MnMigrad;

fn quad(p: &[f64]) -> f64 {
    (p[0] - 1.0).powi(2) // interior minimum at a = 1.0
}

#[test]
fn lower_limited_started_at_lower_bound_escapes() {
    let r = MnMigrad::new()
        .add_lower_limited("a", 0.0, 0.1, 0.0) // start == lower bound
        .minimize(&quad);
    assert!(
        (r.params()[0] - 1.0).abs() < 1e-2,
        "started at lower bound, should reach 1.0, got {}",
        r.params()[0]
    );
}

#[test]
fn upper_limited_started_at_upper_bound_escapes() {
    let r = MnMigrad::new()
        .add_upper_limited("a", 5.0, 0.1, 5.0) // start == upper bound
        .minimize(&quad);
    assert!(
        (r.params()[0] - 1.0).abs() < 1e-2,
        "started at upper bound, should reach 1.0, got {}",
        r.params()[0]
    );
}

#[test]
fn doubly_limited_started_at_lower_bound_escapes() {
    let r = MnMigrad::new()
        .add_limited("a", 0.0, 0.1, 0.0, 5.0) // start == lower bound of [0,5]
        .minimize(&quad);
    assert!(
        (r.params()[0] - 1.0).abs() < 1e-2,
        "started at lower bound of [0,5], should reach 1.0, got {}",
        r.params()[0]
    );
}
