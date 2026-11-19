//! Fixed-seed checks of the real-valued Lean statements in floating point.

use minuit2::migrad::builder::VariableMetricBuilder;
use minuit2::minimum::error::MinimumError;
use minuit2::minimum::gradient::FunctionGradient;
use minuit2::minimum::parameters::MinimumParameters;
use nalgebra::{DMatrix, DVector, SymmetricEigen};
use proptest::prelude::*;
use proptest::test_runner::{Config, RngAlgorithm, RngSeed, TestRunner};
use std::cell::Cell;

const CASES: u32 = 2000;
const SEED: u64 = 0x4446_505f_504f_5344;

// A uniform grid avoids special float values while retaining entries in [-1, 1].
fn entry() -> impl Strategy<Value = f64> {
    (-1_000_000i32..=1_000_000).prop_map(|x| f64::from(x) / 1_000_000.0)
}

#[derive(Clone, Debug)]
struct Case {
    v: DMatrix<f64>,
    dx: DVector<f64>,
    dg: DVector<f64>,
}

#[expect(
    clippy::cast_precision_loss,
    reason = "the dimension is at most six, so the conversion is exact"
)]
fn cases(positive: bool) -> impl Strategy<Value = Case> {
    (2usize..=6)
        .prop_flat_map(|n| {
            (
                Just(n),
                prop::collection::vec(entry(), n * n),
                prop::collection::vec(entry(), n),
                prop::collection::vec(entry(), n),
            )
        })
        .prop_map(|(n, a, dx, dg)| {
            let a = DMatrix::from_row_slice(n, n, &a);
            Case {
                v: &a * a.transpose() + DMatrix::identity(n, n) * n as f64,
                dx: DVector::from_vec(dx),
                dg: DVector::from_vec(dg),
            }
        })
        .prop_filter("curvature angle and nonzero vectors", move |case| {
            let scale = case.dx.norm() * case.dg.norm();
            let signed_dot = case.dx.dot(&case.dg) * if positive { 1.0 } else { -1.0 };
            scale > 0.0 && signed_dot >= 0.1 * scale
        })
}

fn runner() -> TestRunner {
    TestRunner::new(Config {
        cases: CASES,
        rng_seed: RngSeed::Fixed(SEED),
        rng_algorithm: RngAlgorithm::ChaCha,
        max_global_rejects: 20_000,
        failure_persistence: None,
        ..Config::default()
    })
}

fn update(case: &Case) -> DMatrix<f64> {
    let n = case.dx.len();
    let error = MinimumError::new(case.v.clone(), 0.0);
    let p_old = MinimumParameters::new(DVector::zeros(n), 0.0);
    let p_new = MinimumParameters::new(case.dx.clone(), 0.0);
    let g_old = FunctionGradient::analytical(DVector::zeros(n));
    let g_new = FunctionGradient::analytical(case.dg.clone());
    VariableMetricBuilder::update(&error, &p_new, &p_old, &g_new, &g_old).0
}

#[test]
fn posdef_kept_for_positive_delgam() {
    let dimensions = Cell::new(0u8);
    let branches = Cell::new(0u8);
    let count = Cell::new(0u32);
    runner()
        .run(&cases(true), |case| {
            let delgam = case.dx.dot(&case.dg);
            let gvg = case.dg.dot(&(&case.v * &case.dg));
            dimensions.set(dimensions.get() | (1 << (case.dx.len() - 2)));
            branches.set(branches.get() | if delgam > gvg { 1 } else { 2 });
            count.set(count.get() + 1);
            let updated = update(&case);
            let eigenvalues = SymmetricEigen::new(updated.clone()).eigenvalues;
            let smallest = eigenvalues.min();
            let largest = eigenvalues.max();
            // The floor 1e-12 is about 4500 times the f64 rounding unit. The generated matrices
            // reach a condition number of 704 in these cases; the cap 1e4 below makes a change of
            // the generator that leaves this range fail here, where the floor needs a new look.
            prop_assert!(
                smallest > 1e-12 * largest,
                "eigenvalue range: {smallest}..{largest}"
            );
            prop_assert!(
                largest / smallest < 1e4,
                "condition number: {}",
                largest / smallest
            );
            let residual = (&updated * &case.dg - &case.dx).norm();
            // The secant residual is rounding only: the worst case here is 5.4e-15 of |dx|,
            // and the bound leaves a factor of about 200.
            prop_assert!(
                residual <= 1e-12 * case.dx.norm(),
                "secant residual: {residual}"
            );
            Ok(())
        })
        .unwrap();
    assert_eq!(count.get(), CASES);
    assert_eq!(dimensions.get(), 0b1_1111);
    assert_eq!(branches.get(), 0b11);
}

#[test]
fn posdef_lost_for_negative_delgam() {
    let dimensions = Cell::new(0u8);
    let count = Cell::new(0u32);
    runner()
        .run(&cases(false), |case| {
            dimensions.set(dimensions.get() | (1 << (case.dx.len() - 2)));
            count.set(count.get() + 1);
            let updated = update(&case);
            let delgam = case.dx.dot(&case.dg);
            let quadratic = case.dg.dot(&(&updated * &case.dg));
            // In real arithmetic dg.V'dg equals delgam. The f64 value comes from chained
            // products that cancel, and delgam is at least 0.1 |dx| |dg|. The worst deviation
            // here is 2.2e-14 of delgam, so 1e-9 leaves a wide margin.
            prop_assert!(
                (quadratic - delgam).abs() <= 1e-9 * delgam.abs(),
                "quadratic value: {quadratic}, delgam: {delgam}"
            );
            let smallest = SymmetricEigen::new(updated).eigenvalues.min();
            prop_assert!(smallest < 0.0, "smallest eigenvalue: {smallest}");
            Ok(())
        })
        .unwrap();
    assert_eq!(count.get(), CASES);
    assert_eq!(dimensions.get(), 0b1_1111);
}
