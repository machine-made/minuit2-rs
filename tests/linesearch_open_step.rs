//! Regression: after the line search opens the step from the best point, the
//! bounds of the next parabola step must move with that point.
//!
//! f(x) = x^4/4 - 10 x^3/3 - 23 x^2/2 - 12 x. From x = 0 with step 1 and
//! gdel -12 the search visits 5 and then 15. The parabola through 1, 5 and 15
//! has its vertex near 17.5. Calling MnLineSearch in ROOT Minuit2 v6-36-08
//! on these inputs ends at x = 10.390954 with f = -2191.655428.
//! Inputs and full-precision output: reports/verification/raw/trace/linesearch_open_step.root.json.

use minuit2::fcn::FCN;
use minuit2::linesearch::mn_linesearch;
use minuit2::minimum::parameters::MinimumParameters;
use minuit2::mn_fcn::MnFcn;
use minuit2::parameter::MinuitParameter;
use minuit2::precision::MnMachinePrecision;
use minuit2::user_transformation::MnUserTransformation;
use nalgebra::DVector;

struct Quartic;

impl FCN for Quartic {
    fn value(&self, p: &[f64]) -> f64 {
        let x = p[0];
        x.powi(4) / 4.0 - 10.0 * x.powi(3) / 3.0 - 23.0 * x * x / 2.0 - 12.0 * x
    }
}

#[test]
fn open_step_moves_the_parabola_bounds_with_the_best_point() {
    let trafo = MnUserTransformation::new(vec![MinuitParameter::new(0, "x", 0.0, 0.1)]);
    let fcn = MnFcn::new(&Quartic, &trafo);
    let start = MinimumParameters::new(DVector::from_vec(vec![0.0]), 0.0);
    let step = DVector::from_vec(vec![1.0]);
    // f'(0) = -12, so gdel = step . gradient = -12.
    let gdel = -12.0;
    let prec = MnMachinePrecision::new();

    let result = mn_linesearch(&fcn, &start, &step, gdel, &prec);

    // 1e-4 covers the six-decimal reference rounding (5e-7) and small arithmetic drift;
    // it still excludes the former endpoint 10.345059, which missed by 0.046.
    assert!(
        (result.x - 10.390954).abs() < 1e-4,
        "line search ended at x {}, f {}; the reference ends at x 10.390954, f -2191.655",
        result.x,
        result.y
    );
    // -2150 separates the stale endpoint's -2103.33 from ROOT's -2191.66.
    // Near ROOT's endpoint |f'| is about 209, so the 1e-4 x bound costs only 0.021 in f.
    assert!(
        result.y < -2150.0,
        "line search ended at f {}, x {}; the stale bound gives -2103.33",
        result.y,
        result.x
    );
}

#[test]
fn zero_coordinate_does_not_pin_the_minimum_step() {
    let trafo = MnUserTransformation::new(vec![
        MinuitParameter::new(0, "x", 0.0, 0.1),
        MinuitParameter::new(1, "y", 1e10, 0.1),
    ]);
    let objective = |p: &[f64]| (p[0] - 4.0).powi(2);
    let fcn = MnFcn::new(&objective, &trafo);
    let start = MinimumParameters::new(DVector::from_vec(vec![0.0, 1e10]), 16.0);
    let step = DVector::from_vec(vec![1.0, 1.0]);
    let result = mn_linesearch(&fcn, &start, &step, -8.0, &MnMachinePrecision::new());

    // ROOT MnLineSearch.cxx:78 resets a zero slamin to the next ratio, 1e10.
    // slamin = 1e10 * eps2 > 5, so only the unit trial can be accepted.
    // The accepted lambda and value are copied exactly from that trial.
    assert_eq!(result.x, 1.0);
    assert_eq!(result.y, 9.0);
    assert_eq!(fcn.num_of_calls(), 1);
}

// An oscillatory valley needs both shortening and three-point refinements.
#[test]
fn shortening_and_refinement_share_the_root_trial_budget() {
    let (a, b, c, x0) = (
        8.542660106212638,
        6140.499261869658,
        0.25126429994345456,
        6.080884527495888,
    );
    let objective = |p: &[f64]| (p[0] - a).powi(2) + b * (c * p[0]).sin();
    let derivative = 2.0 * (x0 - a) + b * c * (c * x0).cos();
    let trafo = MnUserTransformation::new(vec![MinuitParameter::new(0, "x", x0, 0.1)]);
    let fcn = MnFcn::new(&objective, &trafo);
    let start = MinimumParameters::new(DVector::from_vec(vec![x0]), objective(&[x0]));
    let step = DVector::from_vec(vec![-derivative]);
    mn_linesearch(
        &fcn,
        &start,
        &step,
        -derivative * derivative,
        &MnMachinePrecision::new(),
    );

    // ROOT MnLineSearch.cxx:68-70 starts niter at 1 and stops at 12: at most 11 new calls.
    assert!(
        fcn.num_of_calls() <= 11,
        "search made {} calls",
        fcn.num_of_calls()
    );
}
