//! Simplex seed generator.
//!
//! Creates the initial `MinimumSeed` for the Simplex minimizer by evaluating
//! the FCN at the starting point and computing an initial gradient estimate.

use nalgebra::DVector;

use crate::gradient::InitialGradientCalculator;
use crate::minimum::error::MinimumError;
use crate::minimum::parameters::MinimumParameters;
use crate::minimum::seed::MinimumSeed;
use crate::minimum::state::MinimumState;
use crate::mn_fcn::MnFcn;
use crate::strategy::MnStrategy;
use crate::user_transformation::MnUserTransformation;

pub struct SimplexSeedGenerator;

impl SimplexSeedGenerator {
    pub fn generate(
        fcn: &MnFcn,
        trafo: &MnUserTransformation,
        strategy: &MnStrategy,
    ) -> MinimumSeed {
        let n = trafo.variable_parameters();

        // Get initial internal parameter values
        let int_values = trafo.initial_internal_values();
        let int_vec = DVector::from_vec(int_values.clone());

        // Evaluate FCN at starting point (MnFcn handles int→ext transform)
        let fval = fcn.call(&int_values);

        // Build initial parameters
        let params = MinimumParameters::new(int_vec, fval);

        // Compute initial gradient (heuristic, no FCN calls)
        let grad_calc = InitialGradientCalculator::new(*strategy);
        let gradient = grad_calc.compute(fcn, &params, trafo);

        // Build diagonal covariance from 1/g2
        let mut diag = nalgebra::DMatrix::zeros(n, n);
        let eps2 = trafo.precision().eps2();
        for i in 0..n {
            let g2i = gradient.g2()[i];
            diag[(i, i)] = if g2i.abs() > eps2 { 1.0 / g2i } else { 1.0 };
        }

        let error = MinimumError::new(diag, 1.0);

        // EDM = gradient^T * error_matrix * gradient
        let edm = {
            let g = gradient.grad();
            let e = error.matrix();
            let tmp = e * g;
            g.dot(&tmp)
        };

        let state = MinimumState::new(params, error, gradient, edm, fcn.num_of_calls());

        MinimumSeed::new(state, trafo.clone())
    }
}
