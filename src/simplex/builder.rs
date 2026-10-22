//! Nelder-Mead simplex iteration (Minuit variant).
//!
//! This is not textbook Nelder-Mead: it keeps the Minuit-style adaptive
//! extrapolation step and final centroid check used by the public behavior
//! tests and ROOT parity workloads.

use nalgebra::DVector;

use super::parameters::SimplexParameters;
use crate::minimum::parameters::MinimumParameters;
use crate::minimum::seed::MinimumSeed;
use crate::minimum::state::MinimumState;
use crate::mn_fcn::MnFcn;

pub struct SimplexBuilder;

impl SimplexBuilder {
    pub fn minimum(
        fcn: &MnFcn,
        seed: &MinimumSeed,
        maxfcn: usize,
        minedm: f64,
    ) -> Vec<MinimumState> {
        let n = seed.n_variable_params();
        let prec = seed.precision();

        let start = seed.parameters().vec().clone();
        let mut initial_steps: Vec<f64> =
            (0..n).map(|i| 10.0 * seed.gradient().gstep()[i]).collect();

        const REFLECTION: f64 = 1.0;
        const CONTRACTION: f64 = 0.5;
        const EXPANSION: f64 = 2.0;
        const MIN_EXTRAPOLATION: f64 = 4.0;
        const MAX_EXTRAPOLATION: f64 = 8.0;

        let reflected_weight = 1.0 + REFLECTION;
        let expanded_weight = 1.0 + REFLECTION * EXPANSION;

        // Build initial simplex: N+1 vertices
        let mut vertices: Vec<(f64, Vec<f64>)> = Vec::with_capacity(n + 1);
        vertices.push((seed.fval(), start.as_slice().to_vec()));

        let mut trial_vertex = start.as_slice().to_vec();
        for i in 0..n {
            let min_step = 8.0 * prec.eps2() * (trial_vertex[i].abs() + prec.eps2());
            if initial_steps[i] < min_step {
                initial_steps[i] = min_step;
            }
            trial_vertex[i] += initial_steps[i];
            let fval = fcn.call(&trial_vertex);
            vertices.push((fval, trial_vertex.clone()));
            trial_vertex[i] -= initial_steps[i];
        }

        let mut simplex = SimplexParameters::new(vertices);
        let mut previous_edm;
        loop {
            let worst_index = simplex.jhigh();
            let best_value = simplex.fval_best();
            previous_edm = simplex.edm();

            let centroid = Self::centroid_without(&simplex, worst_index, n);
            let worst_vertex = simplex.params()[worst_index].1.clone();
            let reflected = Self::combine(&centroid, reflected_weight, &worst_vertex, -REFLECTION);
            let reflected_value = fcn.call(&reflected);

            if reflected_value > best_value {
                if reflected_value < simplex.params()[worst_index].0 {
                    simplex.update(worst_index, reflected_value, reflected.clone());
                    if worst_index != simplex.jhigh() {
                        if !Self::should_stop(&simplex, previous_edm, minedm, fcn, maxfcn) {
                            continue;
                        }
                        break;
                    }
                }

                let current_worst = simplex.params()[simplex.jhigh()].1.clone();
                let contracted =
                    Self::combine(&current_worst, CONTRACTION, &centroid, 1.0 - CONTRACTION);
                let contracted_value = fcn.call(&contracted);

                if contracted_value > simplex.params()[simplex.jhigh()].0 {
                    break;
                }
                simplex.update(simplex.jhigh(), contracted_value, contracted);
            } else {
                let expanded = Self::combine(&reflected, EXPANSION, &centroid, 1.0 - EXPANSION);
                let expanded_value = fcn.call(&expanded);

                let reflected_delta =
                    (reflected_value - simplex.params()[worst_index].0) * expanded_weight;
                let expanded_delta =
                    (expanded_value - simplex.params()[worst_index].0) * reflected_weight;
                let extrapolation_denominator = reflected_delta - expanded_delta;

                if extrapolation_denominator.abs() < 1e-30 {
                    if expanded_value < simplex.fval_best() {
                        simplex.update(worst_index, expanded_value, expanded);
                    } else {
                        simplex.update(worst_index, reflected_value, reflected);
                    }
                } else {
                    let extrapolation = 0.5
                        * (expanded_weight * reflected_delta - reflected_weight * expanded_delta)
                        / extrapolation_denominator;

                    if extrapolation < MIN_EXTRAPOLATION {
                        if expanded_value < simplex.fval_best() {
                            simplex.update(worst_index, expanded_value, expanded);
                        } else {
                            simplex.update(worst_index, reflected_value, reflected);
                        }
                    } else {
                        let extrapolation = extrapolation.min(MAX_EXTRAPOLATION);
                        let extrapolated = Self::combine(
                            &centroid,
                            extrapolation,
                            &worst_vertex,
                            1.0 - extrapolation,
                        );
                        let extrapolated_value = fcn.call(&extrapolated);

                        if extrapolated_value < simplex.fval_best()
                            && extrapolated_value < expanded_value
                        {
                            simplex.update(worst_index, extrapolated_value, extrapolated);
                        } else if expanded_value < simplex.fval_best() {
                            simplex.update(worst_index, expanded_value, expanded);
                        } else if extrapolated_value > simplex.fval_best() {
                            if expanded_value < simplex.fval_best() {
                                simplex.update(worst_index, expanded_value, expanded);
                            } else {
                                simplex.update(worst_index, reflected_value, reflected);
                            }
                        } else {
                            simplex.update(worst_index, reflected_value, reflected);
                        }
                    }
                }
            }

            if Self::should_stop(&simplex, previous_edm, minedm, fcn, maxfcn) {
                break;
            }
        }

        let worst_index = simplex.jhigh();
        let centroid = Self::centroid_without(&simplex, worst_index, n);
        let centroid_value = fcn.call(&centroid);

        let (final_vec, final_fval) = if centroid_value < simplex.fval_best() {
            simplex.update(simplex.jhigh(), centroid_value, centroid.clone());
            (centroid, centroid_value)
        } else {
            (simplex.best().to_vec(), simplex.fval_best())
        };

        // Compute dirin from simplex spread, scaled by sqrt(up/edm)
        let edm = simplex.edm();
        let up = fcn.up();
        let scale = if edm > f64::MIN_POSITIVE {
            (up / edm).sqrt()
        } else {
            1.0
        };

        let mut dirin = vec![0.0; n];
        for i in 0..n {
            let mut lo = f64::MAX;
            let mut hi = f64::MIN;
            for (_, v) in simplex.params().iter() {
                if v[i] < lo {
                    lo = v[i];
                }
                if v[i] > hi {
                    hi = v[i];
                }
            }
            dirin[i] = (hi - lo) * scale;
        }

        let final_params = MinimumParameters::with_step(
            DVector::from_vec(final_vec),
            DVector::from_vec(dirin),
            final_fval,
        );

        let state = MinimumState::from_params_edm(final_params, edm, fcn.num_of_calls());
        vec![state]
    }

    fn centroid_without(simplex: &SimplexParameters, excluded: usize, n: usize) -> Vec<f64> {
        let weight = 1.0 / n as f64;
        let mut centroid = vec![0.0; n];
        for (i, (_, vertex)) in simplex.params().iter().enumerate() {
            if i == excluded {
                continue;
            }
            for j in 0..n {
                centroid[j] += weight * vertex[j];
            }
        }
        centroid
    }

    fn combine(left: &[f64], left_scale: f64, right: &[f64], right_scale: f64) -> Vec<f64> {
        left.iter()
            .zip(right)
            .map(|(l, r)| left_scale * l + right_scale * r)
            .collect()
    }

    fn should_stop(
        simplex: &SimplexParameters,
        edm_prev: f64,
        minedm: f64,
        fcn: &MnFcn,
        maxfcn: usize,
    ) -> bool {
        if fcn.num_of_calls() >= maxfcn {
            return true;
        }
        // Both current and previous EDM must be below threshold
        simplex.edm() <= minedm && edm_prev <= minedm
    }
}
