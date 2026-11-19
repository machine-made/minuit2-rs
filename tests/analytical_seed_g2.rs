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
        // dxdi squared times g2 rounds as a product of two doubles.
        // The bound 1e-10 sits well above one ulp of that product.
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
    // The seed copies the start. The bound 1e-15 is a few ulp of a value near 1.
    assert!(
        (free_seed.parameters().vec()[0] - free_start[0]).abs() < 1e-15,
        "a positive g2 must leave x at the start"
    );
    // An unbounded parameter has dxdi = 1, so g2 is the supplied value.
    // The bound 1e-12 is about 2000 ulp of 6 (ulp 8.9e-16).
    assert!((free_seed.gradient().g2()[0] - 6.0).abs() < 1e-12);
    // Same identity transform. The bound 1e-12 is about 1000 ulp of 10 (ulp 1.8e-15).
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
    // The seed copies the start. The bound 1e-12 is a few thousand ulp of an angle near 1.
    assert!(
        (limited_seed.parameters().vec()[0] - limited_start[0]).abs() < 1e-12,
        "a positive g2 must leave x at the start"
    );
    let x_ext = limited.int2ext(0, limited_seed.parameters().vec()[0]);
    let closed = (x_ext - lo) * (up - x_ext) * 6.0;
    // (x-lo)*(up-x)*g2 rounds as a product of two doubles.
    // The bound 1e-10 sits above the sine-map rounding of that product.
    assert!(
        (limited_seed.gradient().g2()[0] - closed).abs() < 1e-10,
        "g2[0] = {}, (x-lo)*(up-x)*g2_ext = {closed}",
        limited_seed.gradient().g2()[0]
    );
    // y is unbounded, so dxdi is 1 and g2 is the supplied 10.
    // The bound 1e-12 is about 1000 ulp of 10 (ulp 1.8e-15).
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
    // The analytical probe has length 1, so a real move is far above 1e-3.
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
    // (x-lo)*(up-x)*g2 after the move rounds as a product of two doubles.
    // The bound 1e-9 covers that product at the escaped point.
    assert!(
        (limited_seed.gradient().g2()[0] - closed).abs() < 1e-9,
        "g2[0] = {}, (x-lo)*(up-x)*g2_ext = {closed}",
        limited_seed.gradient().g2()[0]
    );
}

/// Second-derivative supply. The g2 and the Hessian diagonal differ.
#[derive(Clone, Copy)]
enum Supply {
    G2,
    Hessian,
    Both,
    Neither,
}

/// f = 2 x^2 + 4 y^2. The g2 is 6 and 10. The Hessian diagonal is 4 and 8.
struct SplitSupply {
    supply: Supply,
}

impl SplitSupply {
    fn user_g2() -> [f64; 2] {
        [6.0, 10.0]
    }

    fn hessian_diag() -> [f64; 2] {
        [4.0, 8.0]
    }
}

impl FCN for SplitSupply {
    fn value(&self, p: &[f64]) -> f64 {
        2.0 * p[0] * p[0] + 4.0 * p[1] * p[1]
    }

    fn has_g2(&self) -> bool {
        matches!(self.supply, Supply::G2 | Supply::Both)
    }

    fn g2(&self, _p: &[f64]) -> Vec<f64> {
        Self::user_g2().to_vec()
    }

    fn has_hessian(&self) -> bool {
        matches!(self.supply, Supply::Hessian | Supply::Both)
    }

    /// Packed upper triangle: h00, h01, h11.
    fn hessian(&self, _p: &[f64]) -> Vec<f64> {
        let d = Self::hessian_diag();
        vec![d[0], 1.5, d[1]]
    }
}

impl FCNGradient for SplitSupply {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        vec![4.0 * p[0], 8.0 * p[1]]
    }
}

fn trafo_pair(limited: bool) -> MnUserTransformation {
    let x = if limited {
        MinuitParameter::with_limits(0, "x", 1.5, 0.1, 0.5, 3.0)
    } else {
        MinuitParameter::new(0, "x", 1.2, 0.1)
    };
    MnUserTransformation::new(vec![x, MinuitParameter::new(1, "y", -0.7, 0.1)])
}

/// Step-size g2, 2 / dirin^2, from first principles and not from the
/// production arithmetic. A free parameter has dirin = error. A parameter
/// with two limits uses the sine transform, int = asin(2 (x - a) / (b - a) - 1),
/// and dirin is half the sum of the two internal steps of plus and minus error.
#[expect(
    clippy::manual_midpoint,
    reason = "keeps the operation order of the ROOT reference"
)]
fn step_size_g2(trafo: &MnUserTransformation, int_i: usize) -> f64 {
    let p = &trafo.parameters()[trafo.ext_of_int(int_i)];
    let (value, error) = (p.value(), p.error());
    let dirin = if p.has_lower_limit() && p.has_upper_limit() {
        let (a, b) = (p.lower_limit(), p.upper_limit());
        let internal = |x: f64| (2.0 * (x - a) / (b - a) - 1.0).asin();
        0.5 * ((internal(value + error) - internal(value)).abs()
            + (internal(value - error) - internal(value)).abs())
    } else {
        error
    };
    2.0 / (dirin * dirin)
}

fn assert_matches_external(seed: &MinimumSeed, trafo: &MnUserTransformation, ext_diag: &[f64]) {
    let g2 = seed.gradient().g2();
    let internal = seed.parameters().vec();
    let start = trafo.initial_internal_values();
    for i in 0..g2.len() {
        // A positive diagonal must leave the parameter at the start.
        // The seed copies that value, so the comparison is exact.
        assert_eq!(internal[i], start[i], "parameter {i} moved");
        let ext = trafo.ext_of_int(i);
        let dxdi = trafo.dint2ext(ext, internal[i]);
        let expected = dxdi * dxdi * ext_diag[ext];
        // Same product the seed stores. Equality is exact for that product.
        assert_eq!(g2[i], expected, "g2[{i}]");
    }
}

fn assert_step_size(seed: &MinimumSeed, trafo: &MnUserTransformation) {
    let g2 = seed.gradient().g2();
    let internal = seed.parameters().vec();
    let start = trafo.initial_internal_values();
    for i in 0..g2.len() {
        assert_eq!(internal[i], start[i], "parameter {i} moved");
        let expected = step_size_g2(trafo, i);
        // The oracle orders the arithmetic differently from the seed, so the two
        // differ by a few units of the last place. 1e-12 relative is far above that
        // and far below any change of the formula.
        assert!(
            (g2[i] - expected).abs() <= 1e-12 * expected,
            "g2[{i}] = {}, expected {expected}",
            g2[i]
        );
        // Reciprocation preserves the 1e-12 relative rounding allowance above.
        // The diagonal of the seed covariance is 1/g2, not the identity value 1.
        let diag = seed.error().matrix()[(i, i)];
        assert!(
            (diag - 1.0 / expected).abs() <= 1e-12 / expected,
            "covariance diagonal {diag}, expected {}",
            1.0 / expected
        );
    }
}

/// Check g2 only, Hessian only, both, and neither, with and without a limit.
#[test]
fn seed_g2_source_for_each_supply() {
    for limited in [false, true] {
        let trafo = trafo_pair(limited);
        let g2_only = MigradSeedGenerator::generate_with_gradient(
            &SplitSupply { supply: Supply::G2 },
            &trafo,
            &MnStrategy::new(1),
        );
        assert_matches_external(&g2_only, &trafo, &SplitSupply::user_g2());

        let hessian_only = MigradSeedGenerator::generate_with_gradient(
            &SplitSupply {
                supply: Supply::Hessian,
            },
            &trafo,
            &MnStrategy::new(1),
        );
        assert_matches_external(&hessian_only, &trafo, &SplitSupply::hessian_diag());

        let both = MigradSeedGenerator::generate_with_gradient(
            &SplitSupply {
                supply: Supply::Both,
            },
            &trafo,
            &MnStrategy::new(1),
        );
        assert_matches_external(&both, &trafo, &SplitSupply::hessian_diag());

        let neither = MigradSeedGenerator::generate_with_gradient(
            &SplitSupply {
                supply: Supply::Neither,
            },
            &trafo,
            &MnStrategy::new(1),
        );
        assert_step_size(&neither, &trafo);
    }
}

/// f = (x^2 - 1)^2 + 4 y^2. The g2 is 6 and 10. The Hessian diagonal is -4 and 8.
struct BothNegativeHessian;

impl FCN for BothNegativeHessian {
    fn value(&self, p: &[f64]) -> f64 {
        (p[0] * p[0] - 1.0).powi(2) + 4.0 * p[1] * p[1]
    }

    fn has_g2(&self) -> bool {
        true
    }

    fn g2(&self, _p: &[f64]) -> Vec<f64> {
        vec![6.0, 10.0]
    }

    fn has_hessian(&self) -> bool {
        true
    }

    fn hessian(&self, _p: &[f64]) -> Vec<f64> {
        vec![-4.0, 0.0, 8.0]
    }
}

impl FCNGradient for BothNegativeHessian {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        vec![4.0 * p[0] * (p[0] * p[0] - 1.0), 8.0 * p[1]]
    }
}

/// Check that the escape uses the g2 after a step when both exist.
#[test]
fn escape_after_step_uses_user_g2_when_both_exist() {
    for limited in [false, true] {
        let x = if limited {
            MinuitParameter::with_limits(0, "x", 0.4, 0.1, 0.0, 2.0)
        } else {
            MinuitParameter::new(0, "x", 0.1, 0.1)
        };
        let trafo = MnUserTransformation::new(vec![x, MinuitParameter::new(1, "y", 0.3, 0.1)]);
        let seed = MigradSeedGenerator::generate_with_gradient(
            &BothNegativeHessian,
            &trafo,
            &MnStrategy::new(1),
        );
        let start = trafo.initial_internal_values();
        let x_end = seed.parameters().vec()[0];
        // The analytical probe has length 1, so a real move is far above 1e-3.
        assert!(
            (x_end - start[0]).abs() > 1e-3,
            "limited={limited}: the escape did not move x, x = {x_end}"
        );
        assert_seed_g2(&seed, &trafo, &[6.0, 10.0]);
        let dxdi = trafo.dint2ext(0, x_end);
        let hessian_at_end = dxdi * dxdi * -4.0;
        // The user g2 is 6 and the Hessian diagonal is -4. The gap is far above 1.
        assert!(
            (seed.gradient().g2()[0] - hessian_at_end).abs() > 1.0,
            "limited={limited}: g2 stayed on the Hessian diagonal"
        );
    }
}

/// A NaN or an infinity in a supplied g2.
struct NonFiniteG2 {
    bad: f64,
}

impl FCN for NonFiniteG2 {
    fn value(&self, p: &[f64]) -> f64 {
        p[0] * p[0] + p[1] * p[1]
    }

    fn has_g2(&self) -> bool {
        true
    }

    fn g2(&self, _p: &[f64]) -> Vec<f64> {
        vec![self.bad, 10.0]
    }
}

impl FCNGradient for NonFiniteG2 {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        vec![2.0 * p[0], 2.0 * p[1]]
    }
}

/// Check that a non-finite supplied g2 keeps the step-size estimate.
#[test]
fn nonfinite_supplied_g2_keeps_step_size() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let trafo = MnUserTransformation::new(vec![
            MinuitParameter::new(0, "x", 1.0, 0.1),
            MinuitParameter::new(1, "y", 0.5, 0.1),
        ]);
        let seed = MigradSeedGenerator::generate_with_gradient(
            &NonFiniteG2 { bad },
            &trafo,
            &MnStrategy::new(1),
        );
        assert_step_size(&seed, &trafo);
    }
}

/// External parameter 0 is fixed. Parameter 1 is variable.
struct FixedSelect {
    g2: Vec<f64>,
}

impl FCN for FixedSelect {
    fn value(&self, p: &[f64]) -> f64 {
        p[1] * p[1]
    }

    fn has_g2(&self) -> bool {
        true
    }

    fn g2(&self, _p: &[f64]) -> Vec<f64> {
        self.g2.clone()
    }
}

impl FCNGradient for FixedSelect {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        vec![0.0, 2.0 * p[1]]
    }
}

fn fixed_trafo() -> MnUserTransformation {
    let mut x = MinuitParameter::new(0, "x", 1.0, 0.1);
    x.fix();
    MnUserTransformation::new(vec![x, MinuitParameter::new(1, "y", 0.2, 0.1)])
}

/// Check a fixed parameter: g2 selects the variable entry, and a short vector does not.
#[test]
fn fixed_parameter_selects_variable_g2() {
    let trafo = fixed_trafo();
    let full = MigradSeedGenerator::generate_with_gradient(
        &FixedSelect {
            g2: vec![100.0, 6.0],
        },
        &trafo,
        &MnStrategy::new(1),
    );
    // One variable parameter. dxdi is 1, so the seed g2 is the external entry 6.
    assert_eq!(full.gradient().g2().len(), 1);
    assert_eq!(full.gradient().g2()[0], 6.0);

    let trafo = fixed_trafo();
    let short = MigradSeedGenerator::generate_with_gradient(
        &FixedSelect { g2: vec![100.0] },
        &trafo,
        &MnStrategy::new(1),
    );
    assert_step_size(&short, &trafo);
}

/// f = (x^2 - 1)^2. The g2 is 12 x^2 - 4.
struct MovingG2;

impl FCN for MovingG2 {
    fn value(&self, p: &[f64]) -> f64 {
        (p[0] * p[0] - 1.0).powi(2)
    }

    fn has_g2(&self) -> bool {
        true
    }

    fn g2(&self, p: &[f64]) -> Vec<f64> {
        vec![12.0 * p[0] * p[0] - 4.0]
    }
}

impl FCNGradient for MovingG2 {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        vec![4.0 * p[0] * (p[0] * p[0] - 1.0)]
    }
}

/// Check that a position-dependent g2 is evaluated again after the escape.
#[test]
fn position_dependent_g2_follows_the_moved_point() {
    let x0 = 0.1;
    let trafo = MnUserTransformation::new(vec![MinuitParameter::new(0, "x", x0, 0.1)]);
    let seed = MigradSeedGenerator::generate_with_gradient(&MovingG2, &trafo, &MnStrategy::new(1));
    let x_end = seed.parameters().vec()[0];
    // The analytical probe has length 1, so a real move is far above 1e-3.
    assert!(
        (x_end - x0).abs() > 1e-3,
        "the escape did not move x, x = {x_end}"
    );
    let at_end = 12.0 * x_end * x_end - 4.0;
    let at_start = 12.0 * x0 * x0 - 4.0;
    // Both sides are 12 x^2 - 4 at the same x. The bound 1e-10 covers that product.
    assert!(
        (seed.gradient().g2()[0] - at_end).abs() < 1e-10,
        "g2 = {}, g2(x) = {at_end}",
        seed.gradient().g2()[0]
    );
    // g2 = 12 x^2 - 4 changes by more than 1 between x = 0.1 and a point near 1.
    assert!(
        (seed.gradient().g2()[0] - at_start).abs() > 1.0,
        "g2 stayed at the initial value {at_start}"
    );
}

/// f = 3 x^2 + 2 xy + 2 y^2 has H = [[6,2],[2,4]].
/// H^-1 = [[4,-2],[-2,6]] / 20; the seed stores H^-1, before user scaling.
struct CoupledQuadratic;

impl FCN for CoupledQuadratic {
    fn value(&self, p: &[f64]) -> f64 {
        3.0 * p[0] * p[0] + 2.0 * p[0] * p[1] + 2.0 * p[1] * p[1]
    }
}

impl FCNGradient for CoupledQuadratic {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        vec![6.0 * p[0] + 2.0 * p[1], 2.0 * p[0] + 4.0 * p[1]]
    }
}

#[test]
fn strategy2_analytical_seed_contains_the_full_inverse_hessian() {
    let trafo = MnUserTransformation::new(vec![
        MinuitParameter::new(0, "x", 0.0, 0.1),
        MinuitParameter::new(1, "y", 0.0, 0.1),
    ]);
    let seed =
        MigradSeedGenerator::generate_with_gradient(&CoupledQuadratic, &trafo, &MnStrategy::new(2));
    assert!(seed.error().is_valid());
    // Successful Hesse constructs MinimumError with the literal dcovar 0.0.
    assert_eq!(seed.error().dcovar(), 0.0);
    let expected = [[0.2, -0.1], [-0.1, 0.3]];
    for (i, row) in expected.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            // Quadratic differences have no truncation error; 1e-9 covers cancellation and inversion rounding.
            assert!(
                (seed.error().matrix()[(i, j)] - value).abs() < 1e-9,
                "seed covariance ({i},{j}) = {}, expected {value}",
                seed.error().matrix()[(i, j)]
            );
        }
    }
}
