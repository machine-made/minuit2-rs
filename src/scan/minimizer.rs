//! `SCAn` minimizer: composes seed generator + brute-force builder.

use super::builder::ScanBuilder;
use super::seed::ScanSeedGenerator;
use crate::minimum::FunctionMinimum;
use crate::mn_fcn::MnFcn;
use crate::strategy::MnStrategy;
use crate::user_transformation::MnUserTransformation;

/// Compose seed generation and coordinate scanning.
/// # Examples
///
/// ```
/// use minuit2::{MinuitParameter, MnStrategy, MnUserTransformation};
/// use minuit2::mn_fcn::MnFcn;
/// let trafo = MnUserTransformation::new(vec![MinuitParameter::new(0, "x", 0.0, 0.1)]);
/// let f = |p: &[f64]| (p[0] - 1.0).powi(2);
/// let fcn = MnFcn::new(&f, &trafo);
/// let strategy = MnStrategy::default();
/// let minimum = minuit2::scan::minimizer::ScanMinimizer::minimize(&fcn, &trafo, &strategy, 1000);
/// assert!(minimum.fval() < 1.0);
/// ```
pub struct ScanMinimizer;

impl ScanMinimizer {
    /// Generate a seed, scan each variable coordinate, and return the resulting minimum.
    pub fn minimize(
        fcn: &MnFcn,
        trafo: &MnUserTransformation,
        strategy: &MnStrategy,
        maxfcn: usize,
    ) -> FunctionMinimum {
        let up = fcn.error_def();
        let seed = ScanSeedGenerator::generate(fcn, trafo, strategy);
        if !seed.is_valid() {
            return FunctionMinimum::new(seed, Vec::new(), up);
        }

        let states = ScanBuilder::minimum(fcn, &seed);
        if fcn.num_of_calls() >= maxfcn {
            FunctionMinimum::with_call_limit(seed, states, up)
        } else {
            FunctionMinimum::new(seed, states, up)
        }
    }
}
