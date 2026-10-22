use minuit2::{MnMigrad, MnMinimize};

#[test]
fn minimize_converges_on_quadratic() {
    let quadratic = |p: &[f64]| p[0] * p[0] + 2.0 * p[1] * p[1];

    let result = MnMinimize::new()
        .add("x", 4.0, 1.0)
        .add("y", -3.0, 1.0)
        .tolerance(0.1)
        .minimize(&quadratic);

    assert!(result.is_valid(), "MnMinimize should converge");
    let params = result.params();
    assert!(
        params[0].abs() < 1e-3,
        "x should be near 0, got {}",
        params[0]
    );
    assert!(
        params[1].abs() < 1e-3,
        "y should be near 0, got {}",
        params[1]
    );

    let migrad = MnMigrad::new()
        .add("x", 4.0, 1.0)
        .add("y", -3.0, 1.0)
        .tolerance(0.1)
        .minimize(&quadratic);
    assert!(migrad.is_valid());

    // Hybrid should be in the same quality regime as pure Migrad on smooth bowls.
    assert!(
        result.fval() < migrad.fval() + 1e-8,
        "MnMinimize fval={} should be comparable to MnMigrad fval={}",
        result.fval(),
        migrad.fval()
    );
}

#[test]
fn minimize_respects_fixed_and_limited_parameters() {
    let constrained = |p: &[f64]| (p[0] - 3.0).powi(2) + (p[1] + 2.0).powi(2);

    let result = MnMinimize::new()
        .add_limited("x", 0.5, 0.5, 0.0, 2.0)
        .add("y", -3.0, 0.5)
        .fix(1) // keep y fixed
        .max_fcn(600)
        .tolerance(0.5)
        .minimize(&constrained);

    assert!(result.is_valid());
    let params = result.params();

    // x optimum is 3 but constrained to [0, 2], so it should end at upper bound.
    assert!(
        (params[0] - 2.0).abs() < 1e-3,
        "x should be at upper bound, got {}",
        params[0]
    );
    // y is fixed to the initial value.
    assert!(
        (params[1] + 3.0).abs() < 1e-12,
        "y should stay fixed at -3, got {}",
        params[1]
    );
}
