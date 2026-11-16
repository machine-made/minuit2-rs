//! Regression: after the line search opens the step from the best point, the
//! bounds of the next parabola step must move with that point.
//!
//! f(x) = x^4/4 - 10 x^3/3 - 23 x^2/2 - 12 x. From x = 0 with step 1 and
//! gdel -12 the search visits 5 and then 15. The parabola through 1, 5 and 15
//! has its vertex near 17.5. The reference search ends at x = 10.390954
//! with f = -2191.655428.

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

    // The reference search ends at x = 10.390954. 1e-4 is the required gap.
    assert!(
        (result.x - 10.390954).abs() < 1e-4,
        "line search ended at x {}, f {}; the reference ends at x 10.390954, f -2191.655",
        result.x,
        result.y
    );
    // -2150: the stale bound gives -2103.33 and ROOT gives -2191.66. The slope
    // there is about -209 per unit of x, so 0.1 from ROOT's point costs at
    // most 21, which leaves f below -2170.
    assert!(
        result.y < -2150.0,
        "line search ended at f {}, x {}; the stale bound gives -2103.33",
        result.y,
        result.x
    );
}
