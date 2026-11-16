//! Seed g2 from the second derivative of the function.

use minuit2::fcn::{FCN, FCNGradient};
use minuit2::migrad::seed::MigradSeedGenerator;
use minuit2::minimum::seed::MinimumSeed;
use minuit2::{MinuitParameter, MnStrategy, MnUserTransformation};

/// f = 3 x^2 + 5 y^2. The diagonal second derivative is 6 and 10.
struct QuadraticG2;

impl FCN for QuadraticG2 {
    fn value(&self, p: &[f64]) -> f64 {
        3.0 * p[0] * p[0] + 5.0 * p[1] * p[1]
    }

    fn has_g2(&self) -> bool {
        true
    }

    fn g2(&self, _p: &[f64]) -> Vec<f64> {
        vec![6.0, 10.0]
    }
}

impl FCNGradient for QuadraticG2 {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        vec![6.0 * p[0], 10.0 * p[1]]
    }
}

/// f = -3 x^2 + 5 y^2. The diagonal second derivative is -6 and 10.
struct ConcaveX;

impl FCN for ConcaveX {
    fn value(&self, p: &[f64]) -> f64 {
        -3.0 * p[0] * p[0] + 5.0 * p[1] * p[1]
    }

    fn has_g2(&self) -> bool {
        true
    }

    fn g2(&self, _p: &[f64]) -> Vec<f64> {
        vec![-6.0, 10.0]
    }
}

impl FCNGradient for ConcaveX {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        vec![-6.0 * p[0], 10.0 * p[1]]
    }
}

fn assert_seed_g2(seed: &MinimumSeed, trafo: &MnUserTransformation, g2_ext: &[f64]) {
    let g2 = seed.gradient().g2();
    let internal = seed.parameters().vec();
    for i in 0..g2.len() {
        let ext = trafo.ext_of_int(i);
        let dxdi = trafo.dint2ext(ext, internal[i]);
        let expected = dxdi * dxdi * g2_ext[ext];
        assert!(
            (g2[i] - expected).abs() < 1e-10,
            "g2[{i}] = {}, dxdi^2 * g2_ext = {expected}",
            g2[i]
        );
    }
}

/// Check the seed g2 for a positive diagonal, with and without a limit.
#[test]
fn seed_g2_equals_dxdi_squared_times_user_g2() {
    let free = MnUserTransformation::new(vec![
        MinuitParameter::new(0, "x", 1.2, 0.1),
        MinuitParameter::new(1, "y", -0.7, 0.1),
    ]);
    let free_seed =
        MigradSeedGenerator::generate_with_gradient(&QuadraticG2, &free, &MnStrategy::new(1));
    let free_start = free.initial_internal_values();
    assert_seed_g2(&free_seed, &free, &[6.0, 10.0]);
    assert!(
        (free_seed.parameters().vec()[0] - free_start[0]).abs() < 1e-15,
        "a positive g2 must leave x at the start"
    );
    assert!((free_seed.gradient().g2()[0] - 6.0).abs() < 1e-12);
    assert!((free_seed.gradient().g2()[1] - 10.0).abs() < 1e-12);

    let (lo, up, x0) = (0.5, 3.0, 1.5);
    let limited = MnUserTransformation::new(vec![
        MinuitParameter::with_limits(0, "x", x0, 0.1, lo, up),
        MinuitParameter::new(1, "y", -0.7, 0.1),
    ]);
    let limited_seed =
        MigradSeedGenerator::generate_with_gradient(&QuadraticG2, &limited, &MnStrategy::new(1));
    let limited_start = limited.initial_internal_values();
    assert_seed_g2(&limited_seed, &limited, &[6.0, 10.0]);
    assert!(
        (limited_seed.parameters().vec()[0] - limited_start[0]).abs() < 1e-12,
        "a positive g2 must leave x at the start"
    );
    let x_ext = limited.int2ext(0, limited_seed.parameters().vec()[0]);
    let closed = (x_ext - lo) * (up - x_ext) * 6.0;
    assert!(
        (limited_seed.gradient().g2()[0] - closed).abs() < 1e-10,
        "g2[0] = {}, (x-lo)*(up-x)*g2_ext = {closed}",
        limited_seed.gradient().g2()[0]
    );
    assert!((limited_seed.gradient().g2()[1] - 10.0).abs() < 1e-12);
}

/// Check that the negative-curvature escape keeps the g2 of the function.
#[test]
fn seed_escape_keeps_user_g2() {
    let free = MnUserTransformation::new(vec![
        MinuitParameter::new(0, "x", 0.4, 0.1),
        MinuitParameter::new(1, "y", 0.3, 0.1),
    ]);
    let free_seed =
        MigradSeedGenerator::generate_with_gradient(&ConcaveX, &free, &MnStrategy::new(1));
    assert_seed_g2(&free_seed, &free, &[-6.0, 10.0]);
    let free_start = free.initial_internal_values();
    assert!(
        (free_seed.parameters().vec()[0] - free_start[0]).abs() > 1e-3,
        "a negative g2 must move x, x = {}",
        free_seed.parameters().vec()[0]
    );

    let (lo, up, x0) = (0.0, 2.0, 0.4);
    let limited = MnUserTransformation::new(vec![
        MinuitParameter::with_limits(0, "x", x0, 0.1, lo, up),
        MinuitParameter::new(1, "y", 0.3, 0.1),
    ]);
    let limited_seed =
        MigradSeedGenerator::generate_with_gradient(&ConcaveX, &limited, &MnStrategy::new(1));
    assert_seed_g2(&limited_seed, &limited, &[-6.0, 10.0]);
    let x_ext = limited.int2ext(0, limited_seed.parameters().vec()[0]);
    let dxdi2 = (x_ext - lo) * (up - x_ext);
    let closed = dxdi2 * -6.0;
    assert!(
        (limited_seed.gradient().g2()[0] - closed).abs() < 1e-9,
        "g2[0] = {}, (x-lo)*(up-x)*g2_ext = {closed}",
        limited_seed.gradient().g2()[0]
    );
}
