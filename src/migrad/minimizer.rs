//! VariableMetricMinimizer: composes seed generator + builder.
//!
//! Orchestrates the Migrad minimization by generating the seed, then running
//! the `VariableMetricBuilder` loop.

use super::builder::VariableMetricBuilder;
use super::seed::MigradSeedGenerator;
use crate::fcn::FCNGradient;
use crate::minimum::FunctionMinimum;
use crate::mn_fcn::MnFcn;
use crate::strategy::MnStrategy;
use crate::user_transformation::MnUserTransformation;

pub struct VariableMetricMinimizer;

impl VariableMetricMinimizer {
    /// Minimize using numerical gradients (central differences).
    pub fn minimize(
        fcn: &MnFcn,
        trafo: &MnUserTransformation,
        strategy: &MnStrategy,
        maxfcn: usize,
        tolerance: f64,
    ) -> FunctionMinimum {
        let up = fcn.error_def();

        // Generate seed: FCN eval + numerical gradient + V₀
        let seed = MigradSeedGenerator::generate(fcn, trafo, strategy);

        if !seed.is_valid() {
            return FunctionMinimum::new(seed, Vec::new(), up);
        }

        // EDM tolerance: F77 Minuit compatibility factor
        let edmval = tolerance * up * 0.002;

        // Run variable-metric iteration
        let states = VariableMetricBuilder::minimum(fcn, &seed, strategy, maxfcn, edmval);

        // Check outcome. ROOT evaluates convergence after the Hesse-verified
        // continuation with the extended budget before reporting a call limit
        // (VariableMetricBuilder.cxx:177-198); a valid state converged inside
        // (maxfcn, 1.3*maxfcn] must therefore not be marked call-limited.
        let nfcn = fcn.num_of_calls();
        if let Some(last) = states.last() {
            if !last.error().is_valid() {
                FunctionMinimum::above_max_edm(seed, states, up)
            } else if last.edm() <= 10.0 * edmval {
                FunctionMinimum::new(seed, states, up)
            } else if nfcn >= maxfcn {
                FunctionMinimum::with_call_limit(seed, states, up)
            } else {
                FunctionMinimum::above_max_edm(seed, states, up)
            }
        } else if nfcn >= maxfcn {
            FunctionMinimum::with_call_limit(seed, states, up)
        } else {
            FunctionMinimum::new(seed, states, up)
        }
    }

    /// Minimize using analytical gradients provided by the user.
    pub fn minimize_with_gradient(
        fcn: &dyn FCNGradient,
        trafo: &MnUserTransformation,
        strategy: &MnStrategy,
        maxfcn: usize,
        tolerance: f64,
    ) -> FunctionMinimum {
        let up = fcn.error_def();

        // Generate seed: FCN eval + analytical gradient + V₀
        let seed = MigradSeedGenerator::generate_with_gradient(fcn, trafo, strategy);

        if !seed.is_valid() {
            return FunctionMinimum::new(seed, Vec::new(), up);
        }

        // EDM tolerance: F77 Minuit compatibility factor
        let edmval = tolerance * up * 0.002;

        // Create a temporary MnFcn for call counting during iteration
        let mn_fcn = MnFcn::new(fcn, trafo);

        // Run variable-metric iteration with analytical gradient calculator
        let states = VariableMetricBuilder::minimum_with_gradient(
            &mn_fcn, fcn, &seed, strategy, maxfcn, edmval,
        );

        // Check outcome; see numerical-gradient path above for the ROOT
        // continuation/call-limit ordering.
        let nfcn = mn_fcn.num_of_calls();
        if let Some(last) = states.last() {
            if !last.error().is_valid() {
                FunctionMinimum::above_max_edm(seed, states, up)
            } else if last.edm() <= 10.0 * edmval {
                FunctionMinimum::new(seed, states, up)
            } else if nfcn >= maxfcn {
                FunctionMinimum::with_call_limit(seed, states, up)
            } else {
                FunctionMinimum::above_max_edm(seed, states, up)
            }
        } else if nfcn >= maxfcn {
            FunctionMinimum::with_call_limit(seed, states, up)
        } else {
            FunctionMinimum::new(seed, states, up)
        }
    }
}
