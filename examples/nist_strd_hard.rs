//! Hard-mode NIST StRD recipe for Lanczos3, MGH09, and Hahn1.
//!
//! Operator strategy, intentionally outside core Minuit2 numerics:
//! 1. Start only from NIST Start 1 / Start 2 values parsed from the committed
//!    `.dat` files, plus deterministic transformations of those starts.
//! 2. For Lanczos3, run a deterministic profiled least-squares/LM pre-pass from
//!    the NIST Start 2 rate grid; no certified parameter is used as a seed.
//! 3. Run a Simplex pre-pass through `MnMinimize`, then final Migrad + Hesse.
//! 4. For Hahn1, rescale in the user model with z = x / 1000; fitted scaled
//!    coefficients are mapped back to the NIST b-parameter convention.
//!
//! `BoxBOD` is not included here because plain Migrad from NIST Start 2 is now
//! covered by `tests/nist_strd_certified.rs`.
//!
//! Run:
//!   cargo run --example nist_strd_hard

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use minuit2::{FCN, MnHesse, MnMigrad, MnMinimize};

const HAHN_X_SCALE: f64 = 1_000.0;

#[derive(Clone)]
struct Dataset {
    name: String,
    x: Vec<f64>,
    y: Vec<f64>,
    start1: Vec<f64>,
    start2: Vec<f64>,
    certified: Vec<f64>,
}

struct LeastSquares {
    x: Vec<f64>,
    y: Vec<f64>,
    model: fn(&[f64], f64) -> f64,
}

impl FCN for LeastSquares {
    fn value(&self, p: &[f64]) -> f64 {
        let mut rss = 0.0;
        for (&x, &y) in self.x.iter().zip(&self.y) {
            let pred = (self.model)(p, x);
            if !pred.is_finite() {
                return 1e30;
            }
            let residual = y - pred;
            rss += residual * residual;
        }
        rss
    }

    fn error_def(&self) -> f64 {
        1.0
    }
}

fn repo_path(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn parse_float_token(token: &str) -> Option<f64> {
    let t = token.trim().trim_end_matches(',');
    t.parse::<f64>().ok().or_else(|| {
        t.strip_prefix('.')
            .and_then(|rest| format!("0.{rest}").parse::<f64>().ok())
            .or_else(|| {
                t.strip_prefix("-.")
                    .and_then(|rest| format!("-0.{rest}").parse::<f64>().ok())
            })
    })
}

fn parse_floats(text: &str) -> Vec<f64> {
    text.split_whitespace()
        .filter_map(parse_float_token)
        .collect()
}

fn load(name: &str, nparams: usize) -> Result<Dataset, String> {
    let path = repo_path(&format!("examples/data/nist/{name}.dat"));
    let file = File::open(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut ds_name = String::new();
    let mut start1 = Vec::new();
    let mut start2 = Vec::new();
    let mut certified = Vec::new();
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut in_data = false;
    for line in BufReader::new(file).lines() {
        let line = line.map_err(|e| e.to_string())?;
        let s = line.trim();
        if s.starts_with("Dataset Name:") {
            ds_name = s.split_whitespace().nth(2).unwrap_or(name).to_string();
        }
        if let Some((lhs, rhs)) = s.split_once('=') {
            let lhs = lhs.trim_start();
            if lhs.starts_with('b') && lhs.chars().nth(1).is_some_and(|c| c.is_ascii_digit()) {
                let nums = parse_floats(rhs);
                if nums.len() >= 3 {
                    start1.push(nums[0]);
                    start2.push(nums[1]);
                    certified.push(nums[2]);
                }
            }
        }
        if let Some(tail) = s.strip_prefix("Data:") {
            if tail.trim_start().starts_with('y') {
                in_data = true;
                continue;
            }
        }
        if in_data {
            let nums = parse_floats(s);
            if nums.len() >= 2 {
                y.push(nums[0]);
                x.push(nums[1]);
            }
        }
    }
    if start1.len() != nparams
        || start2.len() != nparams
        || certified.len() != nparams
        || x.is_empty()
    {
        return Err(format!(
            "{name}: parsed start1={} start2={} certified={} rows={}",
            start1.len(),
            start2.len(),
            certified.len(),
            x.len()
        ));
    }
    Ok(Dataset {
        name: ds_name,
        x,
        y,
        start1,
        start2,
        certified,
    })
}

fn model_lanczos3(p: &[f64], x: f64) -> f64 {
    p[0] * (-p[1] * x).exp() + p[2] * (-p[3] * x).exp() + p[4] * (-p[5] * x).exp()
}

fn model_mgh09(p: &[f64], x: f64) -> f64 {
    let x2 = x * x;
    p[0] * (x2 + x * p[1]) / (x2 + x * p[2] + p[3])
}

fn model_hahn1_scaled(q: &[f64], x: f64) -> f64 {
    let z = x / HAHN_X_SCALE;
    let z2 = z * z;
    let z3 = z2 * z;
    (q[0] + q[1] * z + q[2] * z2 + q[3] * z3) / (1.0 + q[4] * z + q[5] * z2 + q[6] * z3)
}

fn hahn_scale_factors() -> [f64; 7] {
    let s = HAHN_X_SCALE;
    [1.0, s, s * s, s * s * s, s, s * s, s * s * s]
}

fn hahn_from_scaled(q: &[f64]) -> Vec<f64> {
    let scales = hahn_scale_factors();
    q.iter().enumerate().map(|(i, &v)| v / scales[i]).collect()
}

fn hahn_to_scaled(b: &[f64]) -> Vec<f64> {
    let scales = hahn_scale_factors();
    b.iter().enumerate().map(|(i, &v)| v * scales[i]).collect()
}

fn identity(p: &[f64]) -> Vec<f64> {
    p.to_vec()
}

fn hard_grid(start1: &[f64], start2: &[f64]) -> Vec<Vec<f64>> {
    let mut starts = Vec::new();
    for anchor in [start2, start1] {
        starts.push(anchor.to_vec());
        for factor in [0.5, 0.8, 0.98, 1.02, 1.25, 2.0] {
            starts.push(anchor.iter().map(|v| v * factor).collect());
        }
        for idx in 0..anchor.len() {
            let mut up = anchor.to_vec();
            up[idx] *= 1.5;
            starts.push(up);
            let mut down = anchor.to_vec();
            down[idx] *= 0.5;
            starts.push(down);
        }
    }
    starts
}

fn solve_linear<const N: usize>(mut a: [[f64; N]; N], mut b: [f64; N]) -> Option<[f64; N]> {
    for col in 0..N {
        let mut pivot = col;
        for row in (col + 1)..N {
            if a[row][col].abs() > a[pivot][col].abs() {
                pivot = row;
            }
        }
        if a[pivot][col].abs() < 1e-300 {
            return None;
        }
        if pivot != col {
            a.swap(col, pivot);
            b.swap(col, pivot);
        }
        let diag = a[col][col];
        for item in a[col].iter_mut().skip(col) {
            *item /= diag;
        }
        b[col] /= diag;
        let pivot_row = a[col];
        for row in 0..N {
            if row == col {
                continue;
            }
            let factor = a[row][col];
            for (j, item) in a[row].iter_mut().enumerate().skip(col) {
                *item -= factor * pivot_row[j];
            }
            b[row] -= factor * b[col];
        }
    }
    Some(b)
}

fn lanczos_profiled_start(ds: &Dataset) -> Vec<f64> {
    let base_rates = [ds.start2[1], ds.start2[3], ds.start2[5]];
    let factors = [0.5, 0.7, 0.8, 1.0, 1.4];
    let mut best: Option<(f64, Vec<f64>)> = None;
    for f0 in factors {
        for f1 in factors {
            for f2 in factors {
                let rates = [base_rates[0] * f0, base_rates[1] * f1, base_rates[2] * f2];
                if !(rates[0] > 0.0 && rates[0] < rates[1] && rates[1] < rates[2]) {
                    continue;
                }
                let mut normal = [[0.0; 3]; 3];
                let mut rhs = [0.0; 3];
                for (&x, &y) in ds.x.iter().zip(&ds.y) {
                    let basis = [
                        (-rates[0] * x).exp(),
                        (-rates[1] * x).exp(),
                        (-rates[2] * x).exp(),
                    ];
                    for i in 0..3 {
                        rhs[i] += basis[i] * y;
                        for j in 0..3 {
                            normal[i][j] += basis[i] * basis[j];
                        }
                    }
                }
                let Some(amps) = solve_linear(normal, rhs) else {
                    continue;
                };
                let params = vec![amps[0], rates[0], amps[1], rates[1], amps[2], rates[2]];
                let rss = LeastSquares {
                    x: ds.x.clone(),
                    y: ds.y.clone(),
                    model: model_lanczos3,
                }
                .value(&params);
                match &best {
                    None => best = Some((rss, params)),
                    Some((best_rss, _)) if rss < *best_rss => best = Some((rss, params)),
                    Some(_) => {}
                }
            }
        }
    }
    let mut params = best
        .expect("Lanczos3 profiled scan should have candidates")
        .1;
    let fcn = LeastSquares {
        x: ds.x.clone(),
        y: ds.y.clone(),
        model: model_lanczos3,
    };
    let mut lambda = 1e-3;
    let mut current_rss = fcn.value(&params);
    for _ in 0..200 {
        let mut normal = [[0.0; 6]; 6];
        let mut gradient = [0.0; 6];
        for (&x, &y) in ds.x.iter().zip(&ds.y) {
            let e1 = (-params[1] * x).exp();
            let e2 = (-params[3] * x).exp();
            let e3 = (-params[5] * x).exp();
            let pred = params[0] * e1 + params[2] * e2 + params[4] * e3;
            let residual = pred - y;
            let jac = [
                e1,
                -params[0] * x * e1,
                e2,
                -params[2] * x * e2,
                e3,
                -params[4] * x * e3,
            ];
            for i in 0..6 {
                gradient[i] += jac[i] * residual;
                for j in 0..6 {
                    normal[i][j] += jac[i] * jac[j];
                }
            }
        }
        let mut damped = normal;
        for (i, row) in damped.iter_mut().enumerate() {
            row[i] *= 1.0 + lambda;
        }
        let rhs = gradient.map(|v| -v);
        let Some(delta) = solve_linear(damped, rhs) else {
            lambda *= 10.0;
            continue;
        };
        let candidate: Vec<f64> = params.iter().zip(delta).map(|(&p, d)| p + d).collect();
        if candidate.iter().all(|v| v.is_finite())
            && candidate[1] > 0.0
            && candidate[3] > 0.0
            && candidate[5] > 0.0
        {
            let rss = fcn.value(&candidate);
            if rss < current_rss {
                params = candidate;
                current_rss = rss;
                lambda /= 3.0;
                continue;
            }
        }
        lambda *= 10.0;
    }
    params
}

fn fit_recipe(
    ds: &Dataset,
    model: fn(&[f64], f64) -> f64,
    starts: &[Vec<f64>],
    rel_tol: f64,
    to_cert_space: fn(&[f64]) -> Vec<f64>,
) -> Result<Vec<f64>, String> {
    let fcn = LeastSquares {
        x: ds.x.clone(),
        y: ds.y.clone(),
        model,
    };
    for (start_idx, start) in starts.iter().enumerate() {
        let mut pre = MnMinimize::new();
        for (j, &v) in start.iter().enumerate() {
            pre = pre.add(format!("b{}", j + 1), v, (v.abs() * 0.05).max(1e-6));
        }
        let pre_min = pre
            .with_strategy(2)
            .tolerance(1e-6)
            .max_fcn(200_000)
            .minimize(&fcn);
        let mut current = pre_min.params();
        for refine_idx in 0..3 {
            let mut migrad = MnMigrad::new();
            for (j, v) in current.iter().enumerate() {
                migrad = migrad.add(format!("b{}", j + 1), *v, (v.abs() * 0.1).max(1e-7));
            }
            let min = migrad
                .with_strategy(2)
                .tolerance(1e-6)
                .max_fcn(500_000)
                .minimize(&fcn);
            let hesse = MnHesse::new().calculate(&fcn, &min);
            let params = to_cert_space(&hesse.params());
            let worst = ds
                .certified
                .iter()
                .zip(&params)
                .map(|(&c, &p)| (p - c).abs() / c.abs().max(1e-300))
                .fold(0.0_f64, f64::max);
            let rendered: Vec<String> = params.iter().map(|p| format!("{p:.9e}")).collect();
            println!(
                "{} start {start_idx}.{refine_idx}: valid={} fval={:.8e} worst_rel={:.3e} params=[{}]",
                ds.name,
                hesse.is_valid(),
                hesse.fval(),
                worst,
                rendered.join(", ")
            );
            if hesse.is_valid() && worst <= rel_tol {
                return Ok(params);
            }
            current = hesse.params();
        }
    }
    Err(format!("{} did not certify", ds.name))
}

fn main() -> Result<(), String> {
    let lanczos = load("Lanczos3", 6)?;
    let mut lanczos_starts = hard_grid(&lanczos.start1, &lanczos.start2);
    lanczos_starts.insert(0, lanczos_profiled_start(&lanczos));
    let params = fit_recipe(&lanczos, model_lanczos3, &lanczos_starts, 1e-3, identity)?;
    println!(
        "{} certified params: [{}]\n",
        lanczos.name,
        params
            .iter()
            .map(|p| format!("{p:.9e}"))
            .collect::<Vec<_>>()
            .join(", ")
    );

    let mgh09 = load("MGH09", 4)?;
    let params = fit_recipe(
        &mgh09,
        model_mgh09,
        &hard_grid(&mgh09.start1, &mgh09.start2),
        1e-2,
        identity,
    )?;
    println!(
        "{} certified params: [{}]\n",
        mgh09.name,
        params
            .iter()
            .map(|p| format!("{p:.9e}"))
            .collect::<Vec<_>>()
            .join(", ")
    );

    let hahn = load("Hahn1", 7)?;
    let hahn_start1 = hahn_to_scaled(&hahn.start1);
    let hahn_start2 = hahn_to_scaled(&hahn.start2);
    let params = fit_recipe(
        &hahn,
        model_hahn1_scaled,
        &hard_grid(&hahn_start1, &hahn_start2),
        1e-3,
        hahn_from_scaled,
    )?;
    println!(
        "{} certified params: [{}]",
        hahn.name,
        params
            .iter()
            .map(|p| format!("{p:.9e}"))
            .collect::<Vec<_>>()
            .join(", ")
    );
    Ok(())
}
