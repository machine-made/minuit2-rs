use std::cell::Cell;

use minuit2::hesse::calculator;
use minuit2::minimum::error::MinimumError;
use minuit2::minimum::gradient::FunctionGradient;
use minuit2::minimum::parameters::MinimumParameters;
use minuit2::minimum::state::MinimumState;
use minuit2::mn_fcn::MnFcn;
use minuit2::parameter::MinuitParameter;
use minuit2::strategy::MnStrategy;
use minuit2::user_transformation::MnUserTransformation;
use minuit2::{FCN, FCNGradient, MnHesse, MnMigrad};
use nalgebra::{DMatrix, DVector};

struct Quadratic {
    calls: Cell<usize>,
    hessian_calls: Cell<usize>,
    gradient: bool,
    hessian: Option<Vec<f64>>,
}

impl Quadratic {
    fn new(hessian: Option<Vec<f64>>) -> Self {
        Self {
            calls: Cell::new(0),
            hessian_calls: Cell::new(0),
            gradient: true,
            hessian,
        }
    }
}

impl FCN for Quadratic {
    fn value(&self, p: &[f64]) -> f64 {
        self.calls.set(self.calls.get() + 1);
        (p[0] - 1.0).powi(2) + (p[1] + 2.0).powi(2)
    }

    fn error_def(&self) -> f64 {
        0.5
    }

    fn has_gradient(&self) -> bool {
        self.gradient
    }

    fn has_hessian(&self) -> bool {
        self.hessian.is_some()
    }

    fn hessian(&self, p: &[f64]) -> Vec<f64> {
        self.hessian_calls.set(self.hessian_calls.get() + 1);
        assert_eq!(p.len(), 2);
        self.hessian.clone().unwrap()
    }
}

impl FCNGradient for Quadratic {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        vec![2.0 * (p[0] - 1.0), 2.0 * (p[1] + 2.0)]
    }
}

fn migrad() -> MnMigrad {
    MnMigrad::new().add("x", 0.4, 0.1).add("y", -1.0, 0.1)
}

/// Relative comparison for results of floating-point arithmetic.
///
/// The factor 1e-12 is about 4500 times the rounding unit 2.2e-16 of `f64`.
/// A printed run showed exact equality in every comparison except one, which differed by 1.3e-16.
/// The margin allows another math library or compiler to change `sin`, `cos` and `asin` by a few units in the last place.
/// Most matrices here have entries of order one. A wrong factor, index or sign changes a result in its first digits.
///
/// The floor `f64::EPSILON` matters only when the expected value is zero.
/// Without it the bound would be zero. With it the bound is 1e-12 * 2.2e-16, about 2.2e-28 absolute.
/// The expected zeros here are off-diagonal covariances of diagonal Hessians.
/// They come out as exact zero, so the floor only has to accept a signed zero.
/// Do not use this function for expected values below 2.2e-16, where the floor sets the bound.
fn assert_relative(got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-12 * want.abs().max(f64::EPSILON),
        "got {got:.17e}, expected {want:.17e}"
    );
}

/// Relative comparison without the floor, for values far from one, such as 1e-200 or 1e300.
///
/// The factor 1e-12 has the same reason as in `assert_relative`.
/// A scaled inverse differs from the exact one by a few units in the last place, about 1e-15.
fn assert_relative_at_scale(got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-12 * want.abs(),
        "got {got:.17e}, expected {want:.17e}"
    );
}

/// Compare a 2x2 matrix with `unit * pattern`. Every entry gets the bound 1e-12 * `unit`.
///
/// The entries of `pattern` are at most one. An expected zero is therefore compared at the size of the diagonal.
/// A coupled inverse can leave a rounding residue of that size in an off-diagonal entry.
fn assert_scaled_matrix(got: &DMatrix<f64>, unit: f64, pattern: [[f64; 2]; 2]) {
    for (i, row) in pattern.iter().enumerate() {
        for (j, entry) in row.iter().enumerate() {
            let want = unit * entry;
            assert!(
                (got[(i, j)] - want).abs() <= 1e-12 * unit,
                "entry ({i}, {j}): got {:.17e}, expected {want:.17e}",
                got[(i, j)]
            );
        }
    }
}

#[test]
fn hesse_makes_no_calls_with_an_analytical_hessian() {
    let f = Quadratic::new(Some(vec![2.0, 0.0, 2.0]));
    let minimum = migrad().minimize_grad(&f);
    assert!(minimum.is_valid());
    assert_eq!(minimum.nfcn(), f.calls.get());
    let before = f.calls.get();
    let hessian_before = f.hessian_calls.get();
    let hesse = MnHesse::new().calculate(&f, &minimum);
    let extra = f.calls.get() - before;
    assert_eq!(extra, 0, "Hesse made {extra} extra function calls");
    assert_eq!(f.hessian_calls.get(), hessian_before + 1);
    assert_eq!(hesse.nfcn(), minimum.nfcn());
    assert!(hesse.is_valid());
    assert!(hesse.state().error().is_accurate());
    let cov = hesse.state().error().matrix();
    for i in 0..2 {
        for j in 0..2 {
            assert_relative(cov[(i, j)], if i == j { 0.5 } else { 0.0 });
        }
    }
    let errors = MnHesse::new().calculate_errors(&f, &minimum);
    assert_eq!(f.calls.get(), before);
    assert_eq!(errors.nfcn(), minimum.nfcn());
    assert_relative(errors.covariance().unwrap().get(0, 0), 0.5);
}

#[test]
fn strategy_two_seed_and_verification_use_analytical_hessians() {
    let f = Quadratic::new(Some(vec![2.0, 0.0, 2.0]));
    let minimum = migrad().with_strategy(2).minimize_grad(&f);
    assert!(minimum.is_valid());
    assert_eq!(minimum.nfcn(), 2);
    assert_eq!(minimum.nfcn(), f.calls.get());
    assert!(f.hessian_calls.get() >= 3);
    assert_relative(minimum.seed().error().matrix()[(0, 0)], 0.5);
}

#[test]
fn limited_parameter_covariance_matches_stored_root() {
    let f = Quadratic::new(Some(vec![2.0, 0.0, 2.0]));
    let minimum = MnMigrad::new()
        .add_limited("x", 0.4, 0.1, -3.0, 4.0)
        .add("y", -1.0, 0.1)
        .minimize_grad(&f);
    assert!(minimum.is_valid());
    let before = f.calls.get();
    let hesse = MnHesse::new().calculate(&f, &minimum);
    let extra = f.calls.get() - before;
    assert_eq!(extra, 0, "Hesse made {extra} extra function calls");
    assert_eq!(hesse.nfcn(), minimum.nfcn());
    assert!(hesse.is_valid());
    let root: serde_json::Value = serde_json::from_str(include_str!(
        "../reports/verification/raw/ref/quadratic2_no_g2_migrad.json"
    ))
    .unwrap();
    // ROOT's Jacobians cancel in the external covariance for this constant Hessian.
    // Scale ROOT's Up = 1 covariance to this function's Up = 0.5.
    let cov = hesse.user_state().covariance().unwrap();
    for i in 0..2 {
        for j in 0..2 {
            let want = root["covariance"][i][j].as_f64().unwrap() * f.error_def();
            assert_relative(cov.get(i, j), want);
        }
    }
    let trafo = minimum.seed().trafo();
    let jacobian = trafo.dint2ext(0, minimum.state().parameters().vec()[0]);
    assert_relative(hesse.state().gradient().g2()[0], 2.0 * jacobian.powi(2));
    assert_relative(
        hesse.state().error().matrix()[(0, 0)],
        0.5 / jacobian.powi(2),
    );
}

#[test]
fn gradient_without_hessian_keeps_numerical_hesse() {
    let f = Quadratic::new(None);
    let minimum = migrad().minimize_grad(&f);
    let before = f.calls.get();
    let hesse = MnHesse::new().calculate(&f, &minimum);
    assert!(f.calls.get() > before);
    assert_eq!(hesse.nfcn(), f.calls.get());
    assert_eq!(f.hessian_calls.get(), 0);
    assert!(hesse.is_valid());
}

#[test]
fn numerical_gradient_keeps_numerical_hesse_even_with_a_hessian() {
    let f = Quadratic::new(Some(vec![2.0, 0.0, 2.0]));
    let minimum = migrad().minimize(&f);
    let before = f.calls.get();
    let hesse = MnHesse::new().calculate(&f, &minimum);
    assert!(f.calls.get() > before);
    assert_eq!(hesse.nfcn(), f.calls.get());
    assert_eq!(f.hessian_calls.get(), 0);
}

#[test]
fn function_without_gradient_flag_keeps_numerical_hesse() {
    let mut f = Quadratic::new(Some(vec![2.0, 0.0, 2.0]));
    f.gradient = false;
    let minimum = migrad().minimize_grad(&f);
    let before = f.calls.get();
    let hesse = MnHesse::new().calculate(&f, &minimum);
    assert!(f.calls.get() > before);
    assert_eq!(hesse.nfcn(), f.calls.get());
}

fn calculate_at_nonzero_gradient(hessian: Vec<f64>) -> calculator::HesseResult {
    let f = Quadratic::new(Some(hessian));
    let minimum = migrad().with_strategy(0).minimize_grad(&f);
    let state = MinimumState::new(
        MinimumParameters::new(DVector::from_vec(vec![0.4, -1.0]), 1.36),
        MinimumError::new(DMatrix::identity(2, 2), 1.0),
        FunctionGradient::analytical(DVector::from_vec(vec![-1.2, 2.0])),
        42.0,
        0,
    );
    let mn_fcn = MnFcn::new(&f, minimum.seed().trafo());
    let result = calculator::calculate(
        &mn_fcn,
        &state,
        minimum.seed().trafo(),
        &MnStrategy::default(),
        0,
    );
    assert_eq!(mn_fcn.num_of_calls(), 0);
    assert_eq!(result.state.nfcn(), 0);
    result
}

#[test]
fn correlated_hessian_updates_g2_and_edm() {
    let result = calculate_at_nonzero_gradient(vec![4.0, 1.0, 2.0]);
    let state = result.state;
    let cov = state.error().matrix();
    for (i, j, want) in [(0, 0, 2.0 / 7.0), (0, 1, -1.0 / 7.0), (1, 1, 4.0 / 7.0)] {
        assert_relative(cov[(i, j)], want);
    }
    // Exact: the Jacobian of an unrestricted parameter is one, so g2 is the Hessian diagonal.
    assert_eq!(state.gradient().g2().as_slice(), &[4.0, 2.0]);
    assert!(state.gradient().is_analytical());
    assert_relative(state.edm(), 0.5 * (2.88 + 4.8 + 16.0) / 7.0);
    assert!(state.error().is_accurate());
}

#[test]
fn missing_or_malformed_hessian_returns_failed_state() {
    for hessian in [
        vec![],
        vec![2.0],
        vec![2.0, 0.0, 2.0, 1.0],
        vec![2.0, f64::NAN, 2.0],
    ] {
        let result = calculate_at_nonzero_gradient(hessian);
        assert!(result.hesse_failed);
        assert!(result.state.error().hesse_failed());
        assert!(!result.state.error().is_valid());
        // Exact: a failed Hesse hands back the old EDM unchanged, so any recomputation would show.
        assert_eq!(result.state.edm(), 42.0);
        assert!(result.state.gradient().is_analytical());
    }
}

#[test]
fn inversion_failure_returns_original_reciprocal_diagonal() {
    // The exact inverse has entries near 5e310, which `f64` cannot hold, so the inversion fails.
    // The smallest eigenvalue of the normalised matrix is 1e-6, above the threshold of about 1.2e-7.
    // The positive-definiteness step then leaves the matrix unchanged.
    let tiny: f64 = 1e-305;
    let correlation: f64 = 0.999_999;
    // The exact diagonal of the inverse is 1 / (tiny * (1 - correlation^2)), about 10^310.7.
    // The largest `f64` is about 10^308.25, so the failure is genuine and not an effect of the scaling.
    let inverse_log10 = -tiny.log10() - (1.0 - correlation * correlation).log10();
    assert!(inverse_log10 > f64::MAX.log10());
    let result = calculate_at_nonzero_gradient(vec![tiny, correlation * tiny, tiny]);
    assert!(result.invert_failed);
    assert!(!result.made_pos_def);
    assert!(result.state.error().invert_failed());
    assert!(!result.state.error().is_valid());
    // Exact: the failed state keeps the old EDM and the Hessian diagonal as g2.
    assert_eq!(result.state.edm(), 42.0);
    assert_eq!(result.state.gradient().g2().as_slice(), &[tiny, tiny]);
    // Exact: ROOT returns 1 / g2 from one division per entry (MnHesse.cxx:145), so the bits must agree.
    // The off-diagonal entries are exactly zero.
    let matrix = result.state.error().matrix();
    assert_eq!(matrix[(0, 0)], 1.0 / tiny);
    assert_eq!(matrix[(1, 1)], 1.0 / tiny);
    assert_eq!(matrix[(0, 1)], 0.0);
}

#[test]
fn inverse_that_overflows_is_a_failure_not_an_infinite_covariance() {
    // The determinant of diag(1e-320, 1) is 1e-320, which `f64` can hold, but the inverse entry 1e320 it gives is infinite.
    let result = calculate_at_nonzero_gradient(vec![1e-320, 0.0, 1.0]);
    assert!(
        result.invert_failed,
        "covariance {:e}, accurate {}",
        result.state.error().matrix()[(0, 0)],
        result.state.error().is_accurate()
    );
    assert!(!result.state.error().is_valid());
    // Exact: the failed state keeps the old EDM.
    assert_eq!(result.state.edm(), 42.0);
}

/// Constant packed Hessian for any number of unrestricted parameters. The objective must never be evaluated.
struct Packed(Vec<f64>);

impl FCN for Packed {
    fn value(&self, _p: &[f64]) -> f64 {
        panic!("analytical Hesse must not evaluate the objective");
    }
    fn has_gradient(&self) -> bool {
        true
    }
    fn has_hessian(&self) -> bool {
        true
    }
    fn hessian(&self, _p: &[f64]) -> Vec<f64> {
        self.0.clone()
    }
}

/// Run the analytical Hesse for a diagonal Hessian with `n` equal entries and a gradient of ones.
fn calculate_with_diagonal(n: usize, entry: f64) -> calculator::HesseResult {
    let mut packed = vec![0.0; n * (n + 1) / 2];
    for i in 0..n {
        packed[i * (i + 1) / 2 + i] = entry;
    }
    let trafo = MnUserTransformation::new(
        (0..n)
            .map(|i| MinuitParameter::new(i, format!("p{i}"), 0.0, 0.1))
            .collect(),
    );
    let state = MinimumState::new(
        MinimumParameters::new(DVector::zeros(n), 0.0),
        MinimumError::new(DMatrix::identity(n, n), 1.0),
        FunctionGradient::analytical(DVector::from_element(n, 1.0)),
        42.0,
        0,
    );
    let function = Packed(packed);
    let fcn = MnFcn::new(&function, &trafo);
    calculator::calculate(&fcn, &state, &trafo, &MnStrategy::default(), 0)
}

#[test]
#[expect(
    clippy::float_cmp,
    reason = "the Hessian diagonal must be copied exactly"
)]
#[expect(
    clippy::cast_precision_loss,
    reason = "the dimension is at most five, so the conversion is exact"
)]
fn diagonal_hessians_of_any_size_invert_at_extreme_scales() {
    // The determinant of an n x n diagonal matrix is entry^n, so it overflows or underflows earlier for larger n.
    // Each case would give zero, NaN or a failed inversion without the scaling.
    for (n, entry) in [
        (2, 1e200),
        (3, 1e110),
        (4, 1e90),
        (4, 1e-90),
        (5, 1e70),
        (3, 1e-110),
    ] {
        let result = calculate_with_diagonal(n, entry);
        assert!(!result.invert_failed, "n = {n}, entry = {entry:e}");
        assert!(!result.made_pos_def, "n = {n}, entry = {entry:e}");
        assert!(
            result.state.error().is_accurate(),
            "n = {n}, entry = {entry:e}"
        );
        // Exact: the Jacobians of unrestricted parameters are one, so g2 is the Hessian diagonal.
        assert!(result.state.gradient().g2().iter().all(|g2| *g2 == entry));
        let covariance = result.state.error().matrix();
        for i in 0..n {
            for j in 0..n {
                let want = if i == j { 1.0 / entry } else { 0.0 };
                assert!(
                    (covariance[(i, j)] - want).abs() <= 1e-12 / entry,
                    "n = {n}, entry = {entry:e}, ({i}, {j}): got {:e}",
                    covariance[(i, j)]
                );
            }
        }
        // The gradient is ones, so the EDM is 0.5 * n / entry.
        assert_relative_at_scale(result.state.edm(), 0.5 * n as f64 / entry);
    }
}

#[test]
fn tiny_hessian_inverts_as_in_root() {
    // ROOT normalises the diagonal before inverting (MnMatrix.cxx:33-41), so a tiny Hessian is not singular.
    let result = calculate_at_nonzero_gradient(vec![1e-200, 0.0, 1e-200]);
    assert!(!result.invert_failed);
    assert!(!result.made_pos_def);
    assert!(result.state.error().is_accurate());
    // Exact: the Jacobians of unrestricted parameters are one, so g2 is the Hessian diagonal.
    assert_eq!(result.state.gradient().g2().as_slice(), &[1e-200, 1e-200]);
    assert_scaled_matrix(
        result.state.error().matrix(),
        1e200,
        [[1.0, 0.0], [0.0, 1.0]],
    );
    // The gradient is (-1.2, 2.0), so the EDM is 0.5 * (1.44 + 4.0) * 1e200.
    assert_relative_at_scale(result.state.edm(), 2.72e200);
}

#[test]
fn huge_hessian_inverts_without_overflow() {
    // The determinant of diag(1e200, 1e200) is 1e400, which overflows `f64`.
    let result = calculate_at_nonzero_gradient(vec![1e200, 0.0, 1e200]);
    assert!(!result.invert_failed);
    assert!(!result.made_pos_def);
    assert!(result.state.error().is_accurate());
    assert_eq!(result.state.gradient().g2().as_slice(), &[1e200, 1e200]);
    assert_scaled_matrix(
        result.state.error().matrix(),
        1e-200,
        [[1.0, 0.0], [0.0, 1.0]],
    );
    // The gradient is (-1.2, 2.0), so the EDM is 0.5 * (1.44 + 4.0) * 1e-200.
    assert_relative_at_scale(result.state.edm(), 2.72e-200);
}

#[test]
fn huge_correlated_hessian_keeps_its_inverse_pattern() {
    // 1e200 times the matrix of `correlated_hessian_updates_g2_and_edm`, whose inverse is (1/7) [[2, -1], [-1, 4]].
    let result = calculate_at_nonzero_gradient(vec![4e200, 1e200, 2e200]);
    assert!(!result.invert_failed);
    assert!(result.state.error().is_accurate());
    assert_scaled_matrix(
        result.state.error().matrix(),
        1e-200,
        [[2.0 / 7.0, -1.0 / 7.0], [-1.0 / 7.0, 4.0 / 7.0]],
    );
    assert_relative_at_scale(result.state.edm(), 0.5 * (2.88 + 4.8 + 16.0) / 7.0 * 1e-200);
}

#[test]
fn mixed_scale_hessian_keeps_the_unscaled_inverse() {
    // The determinant of diag(1e300, 1e-300) is one, so the plain inversion works.
    // Scaling by the largest entry would underflow the small one. The scaled attempt must fall back to the plain one.
    let result = calculate_at_nonzero_gradient(vec![1e300, 0.0, 1e-300]);
    assert!(!result.invert_failed);
    let covariance = result.state.error().matrix();
    assert_relative_at_scale(covariance[(0, 0)], 1e-300);
    assert_relative_at_scale(covariance[(1, 1)], 1e300);
    // Exact: the inverse of a diagonal matrix has exact zeros off the diagonal.
    assert_eq!(covariance[(0, 1)], 0.0);
}

#[test]
fn indefinite_hessian_is_made_positive_definite() {
    let result = calculate_at_nonzero_gradient(vec![1.0, 2.0, 1.0]);
    assert!(result.made_pos_def);
    assert!(result.state.error().is_made_pos_def());
    // Exact: the status constant for a corrected matrix, not a computed value.
    assert_eq!(result.state.error().dcovar(), 1.0);
    assert!(!result.state.error().is_accurate());
    assert!(
        result
            .state
            .error()
            .matrix()
            .clone()
            .symmetric_eigen()
            .eigenvalues
            .min()
            > 0.0
    );
    // Exact: g2 is saved before the matrix is made positive definite, so it keeps the supplied diagonal.
    assert_eq!(result.state.gradient().g2().as_slice(), &[1.0, 1.0]);
    assert_relative(
        result.state.error().matrix()[(0, 0)],
        2.003 / (2.003_f64.powi(2) - 4.0),
    );
}

#[test]
fn positive_nearly_singular_hessian_is_not_modified() {
    let off_diagonal = 0.999_999_5;
    let result = calculate_at_nonzero_gradient(vec![1.0, off_diagonal, 1.0]);
    assert!(!result.made_pos_def);
    assert!(result.state.error().is_accurate());
    assert_eq!(result.state.error().dcovar(), 0.0);
    let inverse_diagonal = 1.0 / (1.0 - off_diagonal.powi(2));
    assert_relative(result.state.error().matrix()[(0, 0)], inverse_diagonal);
}

#[test]
fn diagonal_shift_alone_preserves_root_accurate_status() {
    let result = calculate_at_nonzero_gradient(vec![-2.0, 0.0, 2.0]);
    assert!(!result.made_pos_def);
    assert!(result.state.error().is_accurate());
    assert_eq!(result.state.error().dcovar(), 0.0);
    let shift = 2.5 + minuit2::precision::MnMachinePrecision::new().eps2();
    assert_relative(result.state.error().matrix()[(0, 0)], 1.0 / (-2.0 + shift));
    // Exact: the saved g2 is not shifted, only the matrix that is inverted.
    assert_eq!(result.state.gradient().g2().as_slice(), &[-2.0, 2.0]);
}

/// One parameter with a constant Hessian entry. The objective must never be evaluated.
struct Scalar(f64);

impl FCN for Scalar {
    fn value(&self, _p: &[f64]) -> f64 {
        panic!("analytical Hesse must not evaluate the objective");
    }
    fn has_gradient(&self) -> bool {
        true
    }
    fn has_hessian(&self) -> bool {
        true
    }
    fn hessian(&self, _p: &[f64]) -> Vec<f64> {
        vec![self.0]
    }
}

#[test]
fn scalar_hessian_uses_root_precision_threshold() {
    let eps = minuit2::precision::MnMachinePrecision::new().eps();
    for (curvature, inverse, corrected) in [
        (-2.0, 1.0, true),
        (1e-200, 1.0, true),
        (2.0, 0.5, false),
        (eps, 1.0 / eps, false),
    ] {
        let trafo = MnUserTransformation::new(vec![MinuitParameter::new(0, "x", 0.0, 0.1)]);
        let state = MinimumState::new(
            MinimumParameters::new(DVector::zeros(1), 0.0),
            MinimumError::new(DMatrix::identity(1, 1), 1.0),
            FunctionGradient::analytical(DVector::from_vec(vec![3.0])),
            42.0,
            0,
        );
        let scalar = Scalar(curvature);
        let fcn = MnFcn::new(&scalar, &trafo);
        let result = calculator::calculate(&fcn, &state, &trafo, &MnStrategy::default(), 0);
        assert_eq!(result.made_pos_def, corrected);
        // Exact: the Jacobian is one, so g2 is the supplied curvature.
        assert_eq!(result.state.gradient().g2()[0], curvature);
        assert_relative(result.state.error().matrix()[(0, 0)], inverse);
        assert_relative(result.state.edm(), 4.5 * inverse);
    }
}

#[test]
fn overflow_after_the_jacobians_returns_the_failed_state() {
    // At the midpoint of the limits the Jacobian is 0.5 * 2e155 = 1e155.
    // Then 1e155 * 2.0 * 1e155 overflows to infinity, although the supplied Hessian is finite.
    let trafo = MnUserTransformation::new(vec![MinuitParameter::with_limits(
        0, "x", 0.0, 0.1, -1e155, 1e155,
    )]);
    let internal = DVector::from_vec(trafo.initial_internal_values());
    assert_relative_at_scale(trafo.dint2ext(0, internal[0]), 1e155);
    let state = MinimumState::new(
        MinimumParameters::new(internal.clone(), 0.0),
        MinimumError::new(DMatrix::identity(1, 1), 1.0),
        FunctionGradient::analytical(DVector::from_vec(vec![3.0])),
        42.0,
        7,
    );
    let scalar = Scalar(2.0);
    let fcn = MnFcn::new(&scalar, &trafo);
    let result = calculator::calculate(&fcn, &state, &trafo, &MnStrategy::default(), 0);
    assert!(
        result.hesse_failed,
        "covariance {:e}, accurate {}",
        result.state.error().matrix()[(0, 0)],
        result.state.error().is_accurate()
    );
    assert!(!result.invert_failed);
    assert!(result.state.error().hesse_failed());
    assert!(!result.state.error().is_valid());
    // Exact: a failed Hesse returns the input parameters, gradient, EDM and call count unchanged.
    assert_eq!(result.state.parameters().vec(), &internal);
    assert_eq!(result.state.gradient().grad().as_slice(), &[3.0]);
    assert!(result.state.gradient().is_analytical());
    assert_eq!(result.state.edm(), 42.0);
    assert_eq!(result.state.nfcn(), 7);
    assert!(
        result
            .state
            .error()
            .matrix()
            .iter()
            .all(|entry| *entry == 0.0)
    );
}

#[test]
fn fixed_parameters_and_limits_preserve_packed_hessian_indexing() {
    struct Coupled;
    impl FCN for Coupled {
        fn value(&self, p: &[f64]) -> f64 {
            let (a, b, c) = (p[0] - 3.0, p[1] - 1.0, p[2] + 2.0);
            4.0 * a * a + 2.0 * a * b + a * c + 2.0 * b * b + 0.5 * b * c + c * c
        }
        fn has_gradient(&self) -> bool {
            true
        }
        fn has_hessian(&self) -> bool {
            true
        }
        fn hessian(&self, p: &[f64]) -> Vec<f64> {
            assert_eq!(p.len(), 3);
            assert_eq!(p[0], 3.0);
            vec![8.0, 2.0, 4.0, 1.0, 0.5, 2.0]
        }
    }
    impl FCNGradient for Coupled {
        fn gradient(&self, p: &[f64]) -> Vec<f64> {
            let (a, b, c) = (p[0] - 3.0, p[1] - 1.0, p[2] + 2.0);
            vec![
                8.0 * a + 2.0 * b + c,
                2.0 * a + 4.0 * b + 0.5 * c,
                a + 0.5 * b + 2.0 * c,
            ]
        }
    }
    let minimum = MnMigrad::new()
        .add("fixed", 3.0, 0.1)
        .add_limited("y", 0.4, 0.1, -3.0, 4.0)
        .add_upper_limited("z", -1.0, 0.1, 0.0)
        .fix(0)
        .minimize_grad(&Coupled);
    let hesse = MnHesse::new().calculate(&Coupled, &minimum);
    assert!(hesse.is_valid());
    assert_eq!(hesse.nfcn(), minimum.nfcn());
    let cov = hesse.user_state().covariance().unwrap();
    for (i, j, inverse) in [(0, 0, 2.0 / 7.75), (0, 1, -0.5 / 7.75), (1, 1, 4.0 / 7.75)] {
        assert_relative(cov.get(i, j), 2.0 * inverse);
    }
}

#[test]
fn internal_gradient_space_still_uses_external_hessian_and_jacobians() {
    struct InternalGradient;
    impl FCN for InternalGradient {
        fn value(&self, _p: &[f64]) -> f64 {
            panic!("analytical Hesse must not evaluate the objective");
        }
        fn has_gradient(&self) -> bool {
            true
        }
        fn has_hessian(&self) -> bool {
            true
        }
        fn grad_parameter_space(&self) -> minuit2::fcn::GradientParameterSpace {
            minuit2::fcn::GradientParameterSpace::Internal
        }
        fn hessian(&self, p: &[f64]) -> Vec<f64> {
            // The limited coordinate passes through `asin` and `sin` and comes back 2.2e-16 off, one unit in the last place.
            // The bound 1e-14 is about 45 times that.
            // The unrestricted coordinate is exact. Passing an internal value instead would move the first coordinate by about 1.
            assert!((p[0] - 1.3).abs() < 1e-14);
            assert!((p[1] + 2.5).abs() < 1e-14);
            vec![4.0, 1.0, 2.0]
        }
    }
    let trafo = MnUserTransformation::new(vec![
        MinuitParameter::with_limits(0, "x", 1.3, 0.1, -3.0, 4.0),
        MinuitParameter::new(1, "y", -2.5, 0.1),
    ]);
    let internal = DVector::from_vec(trafo.initial_internal_values());
    let jacobian = trafo.dint2ext(0, internal[0]);
    let state = MinimumState::new(
        MinimumParameters::new(internal, 1.0),
        MinimumError::new(DMatrix::identity(2, 2), 1.0),
        FunctionGradient::analytical(DVector::from_vec(vec![1.0, 2.0])),
        42.0,
        0,
    );
    let fcn = MnFcn::new(&InternalGradient, &trafo);
    let result = calculator::calculate(&fcn, &state, &trafo, &MnStrategy::default(), 0);
    assert_eq!(result.state.nfcn(), 0);
    assert_relative(result.state.gradient().g2()[0], 4.0 * jacobian.powi(2));
    let cov = result.state.error().matrix();
    assert_relative(cov[(0, 0)], 2.0 / (7.0 * jacobian.powi(2)));
    assert_relative(cov[(0, 1)], -1.0 / (7.0 * jacobian));
    assert_relative(cov[(1, 1)], 4.0 / 7.0);
}
