//! Brute-force SCAn minimum builder.
//!
//! Port of ROOT Minuit2 v6-36-08 `math/minuit2/src/ScanBuilder.cxx`.

use nalgebra::DVector;

use crate::minimum::parameters::MinimumParameters;
use crate::minimum::seed::MinimumSeed;
use crate::minimum::state::MinimumState;
use crate::mn_fcn::MnFcn;

pub struct ScanBuilder;

impl ScanBuilder {
    /// Scan each variable parameter once and keep the best point found.
    pub fn minimum(fcn: &MnFcn, seed: &MinimumSeed) -> Vec<MinimumState> {
        let n = seed.n_variable_params();
        let mut x = seed.parameters().vec().as_slice().to_vec();
        let mut external = seed.trafo().transform(&x);
        let mut amin = seed.fval();
        let mut dirin = vec![0.0; n];

        for i in 0..n {
            let ext = seed.trafo().ext_of_int(i);
            if let Some((best_ext, best_fval)) = Self::scan_parameter(fcn, seed, ext, &mut external)
            {
                if best_fval < amin {
                    amin = best_fval;
                    external[ext] = best_ext;
                    x[i] = seed.trafo().ext2int(ext, best_ext);
                }
            }
            dirin[i] = (2.0 * fcn.up() * seed.error().matrix()[(i, i)]).sqrt();
        }

        let params =
            MinimumParameters::with_step(DVector::from_vec(x), DVector::from_vec(dirin), amin);
        vec![MinimumState::from_params_edm(
            params,
            0.0,
            fcn.num_of_calls(),
        )]
    }

    fn scan_parameter(
        fcn: &MnFcn,
        seed: &MinimumSeed,
        ext: usize,
        external: &mut [f64],
    ) -> Option<(f64, f64)> {
        let p = seed.trafo().parameter(ext);
        let val = p.value();
        let err = p.error();
        let mut low = val - 2.0 * err;
        let mut high = val + 2.0 * err;
        if p.has_lower_limit() {
            low = low.max(p.lower_limit());
        }
        if p.has_upper_limit() {
            high = high.min(p.upper_limit());
        }

        let maxsteps = 41usize;
        let step = (high - low) / (maxsteps - 1) as f64;
        let original = external[ext];
        let mut best = None;
        for k in 0..maxsteps {
            let trial = low + k as f64 * step;
            external[ext] = trial;
            let internal: Vec<f64> = (0..seed.n_variable_params())
                .map(|i| {
                    let e = seed.trafo().ext_of_int(i);
                    seed.trafo().ext2int(e, external[e])
                })
                .collect();
            let fval = fcn.call(&internal);
            if best.is_none_or(|(_, best_fval)| fval < best_fval) {
                best = Some((trial, fval));
            }
        }
        external[ext] = original;
        best
    }
}
