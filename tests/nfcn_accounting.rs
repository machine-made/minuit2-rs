//! Tests reported objective call counts and running call budgets.
use std::cell::Cell;

use minuit2::migrad::seed::MigradSeedGenerator;
use minuit2::mn_fcn::MnFcn;
use minuit2::{
    FCN, FCNGradient, MinuitParameter, MnHesse, MnMigrad, MnStrategy, MnUserTransformation,
};

#[derive(Default)]
struct CountedQuadratic {
    calls: Cell<usize>,
}

impl FCN for CountedQuadratic {
    fn value(&self, p: &[f64]) -> f64 {
        self.calls.set(self.calls.get() + 1);
        (p[0] - 1.0).powi(2) + 4.0 * (p[1] + 2.0).powi(2) + 0.3 * p[0] * p[1]
    }

    fn has_g2(&self) -> bool {
        true
    }

    fn g2(&self, _p: &[f64]) -> Vec<f64> {
        vec![2.0, 8.0]
    }
}

impl FCNGradient for CountedQuadratic {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        vec![
            2.0 * (p[0] - 1.0) + 0.3 * p[1],
            8.0 * (p[1] + 2.0) + 0.3 * p[0],
        ]
    }
}

fn migrad(strategy: u32, max_fcn: usize) -> MnMigrad {
    MnMigrad::new()
        .add("x", 0.4, 0.1)
        .add("y", -1.0, 0.1)
        .tolerance(0.1)
        .with_strategy(strategy)
        .max_fcn(max_fcn)
}

#[test]
fn reported_nfcn_equals_user_calls_after_hesse() {
    let f = CountedQuadratic::default();
    let mut minimum = migrad(1, 420).minimize(&f);
    assert_eq!(minimum.nfcn(), f.calls.get());
    for _ in 0..2 {
        minimum = MnHesse::new().calculate(&f, &minimum);
        assert!(minimum.is_valid());
        assert_eq!(
            minimum.nfcn(),
            f.calls.get(),
            "reported nfcn {} differs from user-counted {} after Migrad then Hesse",
            minimum.nfcn(),
            f.calls.get()
        );
        assert_eq!(minimum.user_state().nfcn(), minimum.nfcn());
    }
}

#[test]
fn strategy2_analytical_seed_counts_user_calls() {
    let f = CountedQuadratic::default();
    let trafo = MnUserTransformation::new(vec![
        MinuitParameter::new(0, "x", 0.4, 0.1),
        MinuitParameter::new(1, "y", -1.0, 0.1),
    ]);
    let seed = MigradSeedGenerator::generate_with_gradient(&f, &trafo, &MnStrategy::new(2));
    assert_eq!(
        seed.nfcn(),
        f.calls.get(),
        "seed must include its Hesse calls"
    );
}

#[test]
fn strategy2_numerical_seed_counts_user_calls() {
    let f = CountedQuadratic::default();
    let trafo = MnUserTransformation::new(vec![
        MinuitParameter::new(0, "x", 0.4, 0.1),
        MinuitParameter::new(1, "y", -1.0, 0.1),
    ]);
    let mn_fcn = MnFcn::new(&f, &trafo);
    let seed = MigradSeedGenerator::generate(&mn_fcn, &trafo, &MnStrategy::new(2));
    assert_eq!(seed.nfcn(), f.calls.get());
}

#[test]
fn strategy2_analytical_minimum_counts_user_calls() {
    let f = CountedQuadratic::default();
    let minimum = migrad(2, 420).minimize_grad(&f);
    assert!(minimum.is_valid());
    assert_eq!(minimum.nfcn(), f.calls.get());
}

#[test]
fn strategy2_numerical_minimum_counts_user_calls() {
    let f = CountedQuadratic::default();
    let minimum = migrad(2, 420).minimize(&f);
    assert!(minimum.is_valid());
    assert_eq!(minimum.nfcn(), f.calls.get());
}

#[test]
fn strategy2_stops_before_iteration_when_seed_exhausts_budget() {
    // ROOT builds the whole seed, then checks max_fcn before any iteration
    // (ModularFunctionMinimizer.cxx:57-59,79-85). The seed can exceed the limit.
    for analytical in [false, true] {
        for max_fcn in [1, 5] {
            let f = CountedQuadratic::default();
            let minimum = if analytical {
                migrad(2, max_fcn).minimize_grad(&f)
            } else {
                migrad(2, max_fcn).minimize(&f)
            };
            assert!(
                minimum.reached_call_limit(),
                "analytical={analytical} max_fcn={max_fcn}"
            );
            assert!(!minimum.is_valid());
            assert!(!minimum.is_above_max_edm());
            assert_eq!(minimum.nfcn(), minimum.seed().nfcn());
            assert_eq!(minimum.nfcn(), f.calls.get());
        }
    }
}

#[test]
fn fixed_parameter_hesse_reports_user_calls() {
    struct FixedQuadratic(Cell<usize>);

    impl FCN for FixedQuadratic {
        fn value(&self, p: &[f64]) -> f64 {
            self.0.set(self.0.get() + 1);
            let (x, y, z) = (p[0], p[1], p[2]);
            x * x + 10.0 * y * y + 100.0 * z * z + 2.0 * x * y + 4.0 * x * z + 8.0 * y * z
        }
    }

    impl FCNGradient for FixedQuadratic {
        fn gradient(&self, p: &[f64]) -> Vec<f64> {
            let (x, y, z) = (p[0], p[1], p[2]);
            vec![
                2.0 * x + 2.0 * y + 4.0 * z,
                2.0 * x + 20.0 * y + 8.0 * z,
                4.0 * x + 8.0 * y + 200.0 * z,
            ]
        }
    }

    let f = FixedQuadratic(Cell::new(0));
    let minimum = MnMigrad::new()
        .add("x", 1.0, 0.1)
        .add("y", 2.0, 0.1)
        .add("z", 3.0, 0.1)
        .fix(0)
        .tolerance(0.1)
        .minimize_grad(&f);
    assert_eq!(minimum.nfcn(), f.0.get());
    let minimum = MnHesse::new().calculate(&f, &minimum);
    assert!(minimum.is_valid());
    assert_eq!(minimum.nfcn(), f.0.get());
}

#[test]
fn hesse_reports_exhausted_running_budget() {
    for remaining in [0, 2] {
        let f = CountedQuadratic::default();
        let minimum = migrad(1, 420).minimize(&f);
        let limit = minimum.nfcn() + remaining;
        let result = MnHesse::new().with_max_calls(limit).calculate(&f, &minimum);
        assert_eq!(result.nfcn(), f.calls.get());
        assert!(result.reached_call_limit());
        assert!(result.state().error().has_reached_call_limit());
        assert!(!result.is_valid());
        // ROOT checks after the center call and a complete diagonal cycle.
        assert_eq!(result.nfcn(), minimum.nfcn() + 3);
        assert!(result.nfcn() > limit);
    }
}

#[test]
fn failed_hesse_errors_keep_running_count() {
    let f = CountedQuadratic::default();
    let minimum = migrad(1, 420).minimize(&f);
    let result = MnHesse::new()
        .with_max_calls(minimum.nfcn() + 2)
        .calculate_errors(&f, &minimum);
    assert!(f.calls.get() > minimum.nfcn());
    assert_eq!(result.nfcn(), f.calls.get());
    assert_eq!(result.errors(), minimum.user_state().errors());
}

#[test]
fn hesse_call_limit_error_is_neither_accurate_nor_valid() {
    let f = CountedQuadratic::default();
    let minimum = migrad(1, 420).minimize(&f);
    assert!(minimum.state().error().is_accurate());
    let result = MnHesse::new()
        .with_max_calls(minimum.nfcn())
        .calculate(&f, &minimum);
    let error = result.state().error();
    assert!(!error.is_accurate());
    assert!(!error.is_valid());
    assert!(error.has_reached_call_limit());
}
