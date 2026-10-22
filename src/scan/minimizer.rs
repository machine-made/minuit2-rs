//! SCAn minimizer: composes seed generator + brute-force builder.

use super::builder::ScanBuilder;
use super::seed::ScanSeedGenerator;
use crate::minimum::FunctionMinimum;
use crate::mn_fcn::MnFcn;
use crate::strategy::MnStrategy;
use crate::user_transformation::MnUserTransformation;

pub struct ScanMinimizer;

impl ScanMinimizer {
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
