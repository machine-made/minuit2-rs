use minuit2::fcn::FCN;
use minuit2::migrad::seed::MigradSeedGenerator;
use minuit2::mn_fcn::MnFcn;
use minuit2::{MinuitParameter, MnStrategy, MnUserTransformation};

struct CoupledSaddle;

impl FCN for CoupledSaddle {
    fn value(&self, p: &[f64]) -> f64 {
        let x = p[0];
        let y = p[1];
        -x - y - 0.5 * x * x - 0.5 * y * y + 10.0 * x * x * y * y
    }
}

#[test]
fn negative_g2_seed_line_search_repairs_one_coordinate_before_recomputing_gradient() {
    let trafo = MnUserTransformation::new(vec![
        MinuitParameter::new(0, "x", 0.0, 0.1),
        MinuitParameter::new(1, "y", 0.0, 0.1),
    ]);
    let fcn = MnFcn::new(&CoupledSaddle, &trafo);
    let seed = MigradSeedGenerator::generate(&fcn, &trafo, &MnStrategy::default());

    assert!(
        seed.parameters().vec()[0] > 0.0,
        "first negative-G2 coordinate should move during the seed escape"
    );
    assert_eq!(
        seed.parameters().vec()[1],
        0.0,
        "ROOT repairs only one coordinate before recomputing the full gradient"
    );
}
