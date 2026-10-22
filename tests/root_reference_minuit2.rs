use minuit2::{FCN, FCNGradient, MnHesse, MnMigrad, MnMinimize, MnSimplex};
use std::cell::Cell;

/// Regression behavior checked against ROOT Minuit2 `math/minuit2/test/testMinuit2.cxx`:
/// Hessian indexing must follow floating parameters when fixed parameters exist.
#[test]
fn root_hessian_external_indexing_numeric() {
    let f = |p: &[f64]| {
        let (x, y, z) = (p[0], p[1], p[2]);
        x * x + 10.0 * y * y + 100.0 * z * z + 2.0 * x * y + 4.0 * x * z + 8.0 * y * z
    };

    let minimum = MnMigrad::new()
        .add("x", 1.0, 0.1)
        .add("y", 2.0, 0.1)
        .add("z", 3.0, 0.1)
        .fix(0) // break ext/int identity mapping
        .minimize(&f);
    assert!(minimum.is_valid());

    let hesse = MnHesse::new().calculate(&f, &minimum);
    assert!(hesse.is_valid());
    let state = hesse.user_state();
    let cov = state.covariance().expect("Hesse should provide covariance");

    // Only floating parameters should appear in covariance (y, z).
    assert_eq!(cov.nrow(), 2);

    // For fixed x, Hessian(y,z) = [[20, 8], [8, 200]].
    // ROOT Minuit2 user covariance convention is 2 * Hessian^-1.
    // So expected covariance = 2 * [[200, -8], [-8, 20]] / 3936.
    let c00 = 400.0 / 3936.0;
    let c01 = -16.0 / 3936.0;
    let c11 = 40.0 / 3936.0;

    assert!(
        (cov.get(0, 0) - c00).abs() < 5e-4,
        "cov(0,0) expected ~{c00}, got {}",
        cov.get(0, 0)
    );
    assert!(
        (cov.get(0, 1) - c01).abs() < 5e-4,
        "cov(0,1) expected ~{c01}, got {}",
        cov.get(0, 1)
    );
    assert!(
        (cov.get(1, 1) - c11).abs() < 5e-4,
        "cov(1,1) expected ~{c11}, got {}",
        cov.get(1, 1)
    );
}

struct QuadraticGrad;

impl FCN for QuadraticGrad {
    fn value(&self, p: &[f64]) -> f64 {
        let (x, y, z) = (p[0], p[1], p[2]);
        x * x + 10.0 * y * y + 100.0 * z * z + 2.0 * x * y + 4.0 * x * z + 8.0 * y * z
    }
}

impl FCNGradient for QuadraticGrad {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        let (x, y, z) = (p[0], p[1], p[2]);
        vec![
            2.0 * x + 2.0 * y + 4.0 * z,
            2.0 * x + 20.0 * y + 8.0 * z,
            4.0 * x + 8.0 * y + 200.0 * z,
        ]
    }
}

/// Same regression objective as above, but minimizing with analytical gradient path.
#[test]
fn root_hessian_external_indexing_with_gradient() {
    let f = QuadraticGrad;
    let minimum = MnMigrad::new()
        .add("x", 1.0, 0.1)
        .add("y", 2.0, 0.1)
        .add("z", 3.0, 0.1)
        .fix(0)
        .minimize_grad(&f);
    assert!(minimum.is_valid());

    let hesse = MnHesse::new().calculate(&f, &minimum);
    let cov = hesse
        .user_state()
        .covariance()
        .expect("Hesse should provide covariance");
    assert_eq!(cov.nrow(), 2);

    // Same expected 2x2 block for (y,z), with ROOT 2*H^-1 convention.
    assert!((cov.get(0, 0) - 400.0 / 3936.0).abs() < 5e-4);
    assert!((cov.get(0, 1) - (-16.0 / 3936.0)).abs() < 5e-4);
    assert!((cov.get(1, 1) - 40.0 / 3936.0).abs() < 5e-4);
}

struct Quadratic2;

impl FCN for Quadratic2 {
    fn value(&self, p: &[f64]) -> f64 {
        let (x, y) = (p[0], p[1]);
        let dx = x - 1.0;
        let dy = y + 2.0;
        dx * dx + 4.0 * dy * dy + 0.3 * x * y
    }
}

struct Rosenbrock2;

impl FCN for Rosenbrock2 {
    fn value(&self, p: &[f64]) -> f64 {
        let (x, y) = (p[0], p[1]);
        let t1 = y - x * x;
        let t2 = 1.0 - x;
        100.0 * t1 * t1 + t2 * t2
    }
}

struct QuadraticNoG2 {
    g2_called: Cell<usize>,
}

impl FCN for QuadraticNoG2 {
    fn value(&self, p: &[f64]) -> f64 {
        let (x, y) = (p[0], p[1]);
        let dx = x - 1.0;
        let dy = y + 2.0;
        dx * dx + dy * dy
    }

    fn has_hessian(&self) -> bool {
        true
    }

    fn hessian(&self, _p: &[f64]) -> Vec<f64> {
        // Packed upper triangle for 2x2 identity-scaled Hessian.
        vec![2.0, 0.0, 2.0]
    }

    fn has_g2(&self) -> bool {
        false
    }

    fn g2(&self, _p: &[f64]) -> Vec<f64> {
        self.g2_called.set(self.g2_called.get() + 1);
        vec![2.0, 2.0]
    }
}

impl FCNGradient for QuadraticNoG2 {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        let (x, y) = (p[0], p[1]);
        vec![2.0 * (x - 1.0), 2.0 * (y + 2.0)]
    }
}

/// Regression for ROOT v6-36-08 parity:
/// Simplex uses tolerance directly as EDM target (scaled by Up),
/// not tolerance*1e-3.
#[test]
fn root_simplex_reference_quadratic2() {
    let f = Quadratic2;
    let min = MnSimplex::new()
        .add("x", 0.4, 0.1)
        .add("y", -1.0, 0.1)
        .tolerance(0.1)
        .minimize(&f);

    assert!(min.is_valid());
    let params = min.params();
    assert!((params[0] - 1.15).abs() < 5e-3, "x mismatch: {}", params[0]);
    assert!((params[1] + 2.05).abs() < 5e-3, "y mismatch: {}", params[1]);
    assert!(
        (min.fval() - (-0.67475)).abs() < 5e-4,
        "fval mismatch: {}",
        min.fval()
    );
    assert!(
        (min.edm() - 0.0334887).abs() < 5e-4,
        "edm mismatch: {}",
        min.edm()
    );
}

/// Regression for ROOT v6-36-08 parity:
/// MnMinimize is not "Simplex then Migrad" by default.
/// It runs Migrad first and falls back only if needed, yielding a covariant minimum here.
#[test]
fn root_minimize_reference_rosenbrock2() {
    let f = Rosenbrock2;
    let min = MnMinimize::new()
        .add("x", 0.0, 0.1)
        .add("y", 0.0, 0.1)
        .tolerance(0.1)
        .minimize(&f);

    assert!(min.is_valid());
    assert!(
        min.user_state().has_covariance(),
        "MnMinimize should end with covariance for this workload"
    );
    assert!(
        min.nfcn() > 50,
        "unexpectedly low call count: {}",
        min.nfcn()
    );

    let params = min.params();
    assert!(
        (params[0] - 0.99415).abs() < 0.02,
        "x mismatch: {}",
        params[0]
    );
    assert!(
        (params[1] - 0.98827).abs() < 0.02,
        "y mismatch: {}",
        params[1]
    );
}

/// Regression for the rosenbrock2 Migrad seed NFCN path against ROOT v6-36-08.
/// ROOT's Hesse-verification loop appends the post-Hesse state for this
/// strategy-1 fit; the Rust path now includes the same verification calls.
#[test]
fn root_migrad_rosenbrock2_nfcn_stays_at_root_parity() {
    let f = Rosenbrock2;
    let min = MnMigrad::new()
        .add("x", 0.0, 0.1)
        .add("y", 0.0, 0.1)
        .tolerance(0.1)
        .minimize(&f);

    assert!(min.is_valid());
    assert_eq!(min.nfcn(), 152);
}

#[test]
fn root_migrad_continuation_budget_does_not_report_call_limit_when_converged() {
    let f = Rosenbrock2;
    let maxfcn = 140;
    let min = MnMigrad::new()
        .add("x", 0.0, 0.1)
        .add("y", 0.0, 0.1)
        .max_fcn(maxfcn)
        .tolerance(0.1)
        .minimize(&f);

    assert!(
        min.nfcn() > maxfcn && min.nfcn() <= (1.3 * maxfcn as f64) as usize,
        "regression must exercise the continuation window: nfcn={} maxfcn={}",
        min.nfcn(),
        maxfcn
    );
    assert!(
        min.is_valid(),
        "continuation-converged fit should be valid: nfcn={} call_limit={} above_edm={} edm={:.17e}",
        min.nfcn(),
        min.reached_call_limit(),
        min.is_above_max_edm(),
        min.edm()
    );
    assert!(
        !min.reached_call_limit(),
        "continuation-converged fit must not report call limit"
    );
}

/// Reference behavior for `NoG2CallsWhenFCHasNoG2`:
/// if FCN reports no G2 support, G2 must not be called.
#[test]
fn root_no_g2_calls_when_fcn_has_no_g2() {
    let f = QuadraticNoG2 {
        g2_called: Cell::new(0),
    };
    let min = MnMigrad::new()
        .add("x", 0.4, 0.1)
        .add("y", -1.0, 0.1)
        .tolerance(0.1)
        .minimize_grad(&f);

    assert!(min.is_valid());
    assert_eq!(
        f.g2_called.get(),
        0,
        "G2() must not be called when has_g2() reports false"
    );
}
