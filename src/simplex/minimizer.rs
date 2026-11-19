//! `SimplexMinimizer`: composes seed generator + builder.
//!
//! Internal minimizer component used by the public `MnSimplex` API.

use super::builder::SimplexBuilder;
use super::seed::SimplexSeedGenerator;
use crate::minimum::FunctionMinimum;
use crate::mn_fcn::MnFcn;
use crate::strategy::MnStrategy;
use crate::user_transformation::MnUserTransformation;

/// Compose seed generation and Nelder–Mead iterations.
/// # Examples
///
/// ```
/// use minuit2::{MinuitParameter, MnStrategy, MnUserTransformation};
/// use minuit2::mn_fcn::MnFcn;
/// let trafo = MnUserTransformation::new(vec![MinuitParameter::new(0, "x", 0.0, 0.1)]);
/// let f = |p: &[f64]| (p[0] - 1.0).powi(2);
/// let fcn = MnFcn::new(&f, &trafo);
/// let strategy = MnStrategy::default();
/// let minimum = minuit2::simplex::minimizer::SimplexMinimizer::minimize(&fcn, &trafo, &strategy, 1000, 1e-8);
/// assert!((minimum.params()[0] - 1.0).abs() < 1e-3);
/// ```
pub struct SimplexMinimizer;

impl SimplexMinimizer {
    /// Generate a seed and run Nelder–Mead with the objective error definition scaling the tolerance.
    pub fn minimize(
        fcn: &MnFcn,
        trafo: &MnUserTransformation,
        strategy: &MnStrategy,
        maxfcn: usize,
        tolerance: f64,
    ) -> FunctionMinimum {
        let up = fcn.error_def();

        // Generate seed (initial point + gradient)
        let seed = SimplexSeedGenerator::generate(fcn, trafo, strategy);

        if !seed.is_valid() {
            return FunctionMinimum::new(seed, Vec::new(), up);
        }

        // The builder EDM target is tolerance scaled by Up.
        // (`0.001` scaling is specific to Migrad's internal criterion, not Simplex.)
        let minedm = tolerance * up;

        // Run Nelder-Mead iteration
        let states = SimplexBuilder::minimum(fcn, &seed, maxfcn, minedm);

        // Check if we hit call limit
        let nfcn = fcn.num_of_calls();
        if nfcn >= maxfcn {
            FunctionMinimum::with_call_limit(seed, states, up)
        } else if let Some(last) = states.last() {
            if last.edm() > minedm {
                FunctionMinimum::above_max_edm(seed, states, up)
            } else {
                FunctionMinimum::new(seed, states, up)
            }
        } else {
            FunctionMinimum::new(seed, states, up)
        }
    }
}
