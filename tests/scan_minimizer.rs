mod common;

use minuit2::MnScanMinimizer;

#[test]
fn quadratic_shared_smoke() {
    let result = MnScanMinimizer::new()
        .add("x", 5.0, 2.5)
        .add("y", -3.0, 1.5)
        .minimize(&|p: &[f64]| p[0] * p[0] + p[1] * p[1]);

    assert!(
        result.is_valid(),
        "SCAn should produce a valid FunctionMinimum"
    );
    common::assert_function_minimum_display(&result, &["x", "y"]);
    assert!(
        result.fval() < result.seed().fval(),
        "SCAn should improve the shared quadratic: start fval={}, final fval={}",
        result.seed().fval(),
        result.fval()
    );
    assert!(
        result.fval() <= 1.0e-24,
        "SCAn grid should hit the quadratic minimum exactly with this setup, got fval={}",
        result.fval()
    );
}

#[test]
fn bounded_parameter_case_respects_limits() {
    let result = MnScanMinimizer::new()
        .add_limited("x", 2.0, 1.5, 0.0, 4.0)
        .minimize(&|p: &[f64]| (p[0] - 3.0).powi(2));

    assert!(result.is_valid(), "bounded SCAn result should be valid");
    let x = result.params()[0];
    assert!((0.0..=4.0).contains(&x), "x={x} should stay inside [0, 4]");
    assert!((x - 3.0).abs() <= 0.05, "x should scan near 3.0, got {x}");
    assert!(
        result.fval() <= 0.0026,
        "bounded scan fval should be small, got {}",
        result.fval()
    );
}

#[test]
fn bad_start_improves_before_migrad() {
    let result = MnScanMinimizer::new()
        .add("x", -50.0, 30.0)
        .add("y", 80.0, 40.0)
        .with_strategy(1)
        .max_fcn(500)
        .minimize(&|p: &[f64]| (p[0] - 10.0).powi(2) + (p[1] + 5.0).powi(2));

    assert!(result.is_valid(), "bad-start SCAn should remain valid");
    assert!(
        result.fval() < result.seed().fval(),
        "bad-start SCAn should improve fval: start={}, final={}",
        result.seed().fval(),
        result.fval()
    );
    assert!(
        result.fval() < 1000.0,
        "bad-start SCAn should make a substantial pre-pass improvement, got fval={}",
        result.fval()
    );
}
