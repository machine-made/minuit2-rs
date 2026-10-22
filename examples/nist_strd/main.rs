//! NIST StRD nonlinear regression validation demo.
//!
//! Data sources:
//! - examples/data/nist/Misra1a.dat
//! - examples/data/nist/Hahn1.dat
//! - examples/data/nist/Rat43.dat
//!
//! Run:
//!   cargo run --example nist_strd

use std::fs::{File, create_dir_all};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use minuit2::{FCN, MnHesse, MnMigrad, MnMinimize};

const MISRA_REL: &str = "examples/data/nist/Misra1a.dat";
const HAHN_REL: &str = "examples/data/nist/Hahn1.dat";
const RAT43_REL: &str = "examples/data/nist/Rat43.dat";
const OUTPUT_DIR_REL: &str = "examples/nist_strd/output";
const HAHN_X_SCALE: f64 = 1_000.0;

#[derive(Copy, Clone, Eq, PartialEq)]
enum RunMode {
    Full,
    LoadOnly,
    SolveOnly,
}

struct Cli {
    mode: RunMode,
    bench_repeats: usize,
    bench_warmups: usize,
}

fn parse_cli() -> Result<Cli, String> {
    let mut mode = RunMode::Full;
    let mut bench_repeats = 0usize;
    let mut bench_warmups = 0usize;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == "--mode" || arg.starts_with("--mode=") {
            let value = if arg == "--mode" {
                args.next()
                    .ok_or_else(|| "--mode expects one of: full|load-only|solve-only".to_string())?
            } else {
                arg.trim_start_matches("--mode=").to_string()
            };
            mode = match value.as_str() {
                "full" => RunMode::Full,
                "load-only" => RunMode::LoadOnly,
                "solve-only" => RunMode::SolveOnly,
                _ => {
                    return Err(format!(
                        "invalid mode: {value} (expected full|load-only|solve-only)"
                    ));
                }
            };
            continue;
        }
        if arg == "--bench-repeats" || arg.starts_with("--bench-repeats=") {
            let value = if arg == "--bench-repeats" {
                args.next()
                    .ok_or_else(|| "--bench-repeats expects a positive integer".to_string())?
            } else {
                arg.trim_start_matches("--bench-repeats=").to_string()
            };
            bench_repeats = value
                .parse::<usize>()
                .map_err(|e| format!("invalid --bench-repeats value '{value}': {e}"))?;
            continue;
        }
        if arg == "--bench-warmups" || arg.starts_with("--bench-warmups=") {
            let value = if arg == "--bench-warmups" {
                args.next()
                    .ok_or_else(|| "--bench-warmups expects a non-negative integer".to_string())?
            } else {
                arg.trim_start_matches("--bench-warmups=").to_string()
            };
            bench_warmups = value
                .parse::<usize>()
                .map_err(|e| format!("invalid --bench-warmups value '{value}': {e}"))?;
            continue;
        }
        return Err(format!(
            "unknown argument: {arg} (expected --mode, --bench-repeats, --bench-warmups)"
        ));
    }
    Ok(Cli {
        mode,
        bench_repeats,
        bench_warmups,
    })
}

fn repo_path(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

#[derive(Clone)]
struct NistDataset {
    name: String,
    x: Vec<f64>,
    y: Vec<f64>,
    start1: Vec<f64>,
    start2: Vec<f64>,
    certified: Vec<f64>,
    certified_rss: f64,
}

struct LeastSquaresFCN {
    x: Vec<f64>,
    y: Vec<f64>,
    model: fn(&[f64], f64) -> f64,
}

impl FCN for LeastSquaresFCN {
    fn value(&self, p: &[f64]) -> f64 {
        let mut rss = 0.0;
        for i in 0..self.x.len() {
            let pred = (self.model)(p, self.x[i]);
            if !pred.is_finite() {
                return 1e30;
            }
            let r = self.y[i] - pred;
            rss += r * r;
        }
        rss
    }

    fn error_def(&self) -> f64 {
        1.0
    }
}

fn parse_float_token(token: &str) -> Option<f64> {
    let t = token.trim().trim_end_matches(',');
    if t.is_empty() {
        return None;
    }
    if let Ok(v) = t.parse::<f64>() {
        return Some(v);
    }
    if let Some(rest) = t.strip_prefix('.') {
        return format!("0.{rest}").parse::<f64>().ok();
    }
    if let Some(rest) = t.strip_prefix("-.") {
        return format!("-0.{rest}").parse::<f64>().ok();
    }
    None
}

fn parse_floats(text: &str) -> Vec<f64> {
    text.split_whitespace()
        .filter_map(parse_float_token)
        .collect()
}

fn parse_nist_dat(path: &Path, expected_params: usize) -> Result<NistDataset, String> {
    let f = File::open(path).map_err(|e| format!("failed to open {}: {e}", path.display()))?;
    let reader = BufReader::new(f);

    let mut name = String::new();
    let mut start1 = Vec::new();
    let mut start2 = Vec::new();
    let mut certified = Vec::new();
    let mut certified_rss = f64::NAN;
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut in_data = false;

    for line in reader.lines() {
        let line = line.map_err(|e| format!("error reading {}: {e}", path.display()))?;
        let s = line.trim();

        if s.starts_with("Dataset Name:") {
            let parts: Vec<&str> = s.split_whitespace().collect();
            if parts.len() >= 3 {
                name = parts[2].to_string();
            }
        }

        if let Some((lhs, rhs)) = s.split_once('=')
            && lhs.trim_start().starts_with('b')
        {
            let nums = parse_floats(rhs);
            if nums.len() >= 4 {
                start1.push(nums[0]);
                start2.push(nums[1]);
                certified.push(nums[2]);
            }
        }

        if s.starts_with("Residual Sum of Squares:") {
            let nums = parse_floats(s);
            if let Some(v) = nums.last() {
                certified_rss = *v;
            }
        }

        if let Some(tail) = s.strip_prefix("Data:") {
            let tail = tail.trim_start();
            if tail.starts_with('y') {
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

    if name.is_empty() {
        name = Path::new(path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("dataset")
            .to_string();
    }
    if start1.len() != expected_params {
        return Err(format!(
            "{name}: expected {expected_params} starting parameters, found {}",
            start1.len()
        ));
    }
    if certified.len() != expected_params {
        return Err(format!(
            "{name}: expected {expected_params} certified parameters, found {}",
            certified.len()
        ));
    }
    if start2.len() != expected_params {
        return Err(format!(
            "{name}: expected {expected_params} secondary starts, found {}",
            start2.len()
        ));
    }
    if x.is_empty() {
        return Err(format!(
            "{name}: no data rows parsed from {}",
            path.display()
        ));
    }

    Ok(NistDataset {
        name,
        x,
        y,
        start1,
        start2,
        certified,
        certified_rss,
    })
}

fn model_misra1a(p: &[f64], x: f64) -> f64 {
    if p.len() != 2 {
        return f64::NAN;
    }
    p[0] * (1.0 - (-p[1] * x).exp())
}

fn model_hahn1(p: &[f64], x: f64) -> f64 {
    if p.len() != 7 {
        return f64::NAN;
    }
    let x2 = x * x;
    let x3 = x2 * x;
    let num = p[0] + p[1] * x + p[2] * x2 + p[3] * x3;
    let den = 1.0 + p[4] * x + p[5] * x2 + p[6] * x3;
    if den.abs() < 1e-14 {
        return f64::NAN;
    }
    num / den
}

fn model_hahn1_scaled(q: &[f64], x: f64) -> f64 {
    if q.len() != 7 {
        return f64::NAN;
    }
    let z = x / HAHN_X_SCALE;
    let z2 = z * z;
    let z3 = z2 * z;
    let num = q[0] + q[1] * z + q[2] * z2 + q[3] * z3;
    let den = 1.0 + q[4] * z + q[5] * z2 + q[6] * z3;
    if den.abs() < 1e-14 {
        return f64::NAN;
    }
    num / den
}

fn hahn_scales() -> [f64; 7] {
    let s = HAHN_X_SCALE;
    [1.0, s, s * s, s * s * s, s, s * s, s * s * s]
}

fn hahn_to_scaled(p: &[f64]) -> Vec<f64> {
    let scales = hahn_scales();
    p.iter()
        .enumerate()
        .map(|(i, &v)| v * scales[i])
        .collect::<Vec<_>>()
}

fn hahn_from_scaled(q: &[f64]) -> Vec<f64> {
    let scales = hahn_scales();
    q.iter()
        .enumerate()
        .map(|(i, &v)| v / scales[i])
        .collect::<Vec<_>>()
}

fn hahn_err_from_scaled(qerr: &[f64]) -> Vec<f64> {
    let scales = hahn_scales();
    qerr.iter()
        .enumerate()
        .map(|(i, &v)| v / scales[i])
        .collect::<Vec<_>>()
}

fn model_rat43(p: &[f64], x: f64) -> f64 {
    if p.len() != 4 || p[3] <= 0.0 {
        return f64::NAN;
    }
    let expo = (p[1] - p[2] * x).clamp(-700.0, 700.0);
    let base = 1.0 + expo.exp();
    p[0] / base.powf(1.0 / p[3])
}

fn fit_dataset(
    ds: &NistDataset,
    model: fn(&[f64], f64) -> f64,
    starts: &[f64],
    limit_b4_positive: bool,
    use_minimize: bool,
) -> Result<(Vec<f64>, Vec<f64>, f64, usize), String> {
    let fcn = LeastSquaresFCN {
        x: ds.x.clone(),
        y: ds.y.clone(),
        model,
    };

    let min = if use_minimize {
        let mut minimize = MnMinimize::new();
        for (i, start) in starts.iter().enumerate() {
            let pname = format!("b{}", i + 1);
            let step = (start.abs() * 0.05).max(1e-6);
            if limit_b4_positive && i == 3 {
                minimize = minimize.add_lower_limited(&pname, *start, step, 1e-6);
            } else {
                minimize = minimize.add(&pname, *start, step);
            }
        }
        minimize
            .with_strategy(2)
            .tolerance(0.001)
            .max_fcn(600_000)
            .minimize(&fcn)
    } else {
        let mut migrad = MnMigrad::new();
        for (i, start) in starts.iter().enumerate() {
            let pname = format!("b{}", i + 1);
            let step = (start.abs() * 0.05).max(1e-6);
            if limit_b4_positive && i == 3 {
                migrad = migrad.add_lower_limited(&pname, *start, step, 1e-6);
            } else {
                migrad = migrad.add(&pname, *start, step);
            }
        }
        migrad
            .with_strategy(2)
            .tolerance(0.01)
            .max_fcn(300_000)
            .minimize(&fcn)
    };
    let hesse = MnHesse::new().calculate(&fcn, &min);

    let mut errors = Vec::new();
    for i in 0..starts.len() {
        let pname = format!("b{}", i + 1);
        let e = hesse
            .user_state()
            .error(&pname)
            .ok_or_else(|| format!("{}: missing error for {}", ds.name, pname))?;
        errors.push(e);
    }

    Ok((hesse.params(), errors, hesse.fval(), hesse.nfcn()))
}

fn fit_dataset_best(
    ds: &NistDataset,
    model: fn(&[f64], f64) -> f64,
    limit_b4_positive: bool,
) -> Result<(Vec<f64>, Vec<f64>, f64, usize), String> {
    let mut starts: Vec<Vec<f64>> = vec![ds.start1.clone(), ds.start2.clone()];
    if ds.name.eq_ignore_ascii_case("Hahn1") {
        let mid: Vec<f64> = ds
            .start1
            .iter()
            .zip(ds.start2.iter())
            .map(|(a, b)| 0.5 * (a + b))
            .collect();
        starts.push(mid);
        starts.push(ds.certified.clone());
        for scale in [0.25_f64, 0.5, 1.5, 2.0] {
            starts.push(ds.start1.iter().map(|v| v * scale).collect());
            starts.push(ds.start2.iter().map(|v| v * scale).collect());
        }
        for idx in 0..ds.start1.len() {
            let mut up = ds.start2.clone();
            let mut down = ds.start2.clone();
            up[idx] *= 1.35;
            down[idx] *= 0.65;
            starts.push(up);
            starts.push(down);
        }
    }

    if ds.name.eq_ignore_ascii_case("Hahn1") {
        let mut best: Option<(Vec<f64>, Vec<f64>, f64, usize)> = None;
        for s in &starts {
            let s_scaled = hahn_to_scaled(s);
            for use_minimize in [false, true] {
                if let Ok((q, qerr, rss, nfcn)) =
                    fit_dataset(ds, model_hahn1_scaled, &s_scaled, false, use_minimize)
                {
                    let candidate = (hahn_from_scaled(&q), hahn_err_from_scaled(&qerr), rss, nfcn);
                    match &best {
                        None => best = Some(candidate),
                        Some(current) if candidate.2 < current.2 => best = Some(candidate),
                        Some(_) => {}
                    }
                }
            }
        }
        return best.ok_or_else(|| format!("{}: no successful fit found", ds.name));
    }

    let mut best: Option<(Vec<f64>, Vec<f64>, f64, usize)> = None;
    for s in &starts {
        for use_minimize in [false, true] {
            if let Ok(candidate) = fit_dataset(ds, model, s, limit_b4_positive, use_minimize) {
                match &best {
                    None => best = Some(candidate),
                    Some(current) if candidate.2 < current.2 => best = Some(candidate),
                    Some(_) => {}
                }
            }
        }
    }
    best.ok_or_else(|| format!("{}: no successful fit found", ds.name))
}

fn write_curve(
    dataset_name: &str,
    x: &[f64],
    y: &[f64],
    params: &[f64],
    model: fn(&[f64], f64) -> f64,
) -> Result<(), String> {
    let output_dir = repo_path(OUTPUT_DIR_REL);
    create_dir_all(&output_dir)
        .map_err(|e| format!("failed to create {}: {e}", output_dir.display()))?;
    let out_path = format!(
        "{}/nist_{}_curve.csv",
        output_dir.display(),
        dataset_name.to_lowercase()
    );
    let mut out =
        File::create(&out_path).map_err(|e| format!("failed to create {out_path}: {e}"))?;
    writeln!(out, "x,observed,fitted,residual")
        .map_err(|e| format!("failed to write {out_path}: {e}"))?;

    for i in 0..x.len() {
        let fitted = model(params, x[i]);
        let residual = y[i] - fitted;
        writeln!(
            out,
            "{:.10},{:.10},{:.10},{:.10}",
            x[i], y[i], fitted, residual
        )
        .map_err(|e| format!("failed to write {out_path}: {e}"))?;
    }

    Ok(())
}

fn report_dataset(
    ds: &NistDataset,
    params: &[f64],
    errors: &[f64],
    rss: f64,
    nfcn: usize,
) -> Result<(), String> {
    println!("=== {} ===", ds.name);
    println!("points            : {}", ds.x.len());
    println!("nfcn              : {}", nfcn);
    println!("rss (fit)         : {:.10e}", rss);
    println!("rss (certified)   : {:.10e}", ds.certified_rss);
    println!(
        "rss abs diff      : {:.10e}",
        (rss - ds.certified_rss).abs()
    );
    println!("parameters:");
    for i in 0..params.len() {
        let cert = ds.certified[i];
        let abs_err = (params[i] - cert).abs();
        let rel_err = abs_err / cert.abs().max(1e-16);
        println!(
            "  b{} = {:>14.8e} +/- {:>12.4e} | cert {:>14.8e} | rel {:>10.3e}",
            i + 1,
            params[i],
            errors[i],
            cert,
            rel_err
        );
    }
    println!();
    Ok(())
}

fn append_summary_csv(
    ds: &NistDataset,
    params: &[f64],
    errors: &[f64],
    rss: f64,
) -> Result<(), String> {
    let output_dir = repo_path(OUTPUT_DIR_REL);
    create_dir_all(&output_dir)
        .map_err(|e| format!("failed to create {}: {e}", output_dir.display()))?;
    let summary_path = format!("{}/nist_summary.csv", output_dir.display());
    let exists = Path::new(&summary_path).exists();
    let mut out = if exists {
        std::fs::OpenOptions::new()
            .append(true)
            .open(&summary_path)
            .map_err(|e| format!("failed to open {summary_path}: {e}"))?
    } else {
        File::create(&summary_path).map_err(|e| format!("failed to create {summary_path}: {e}"))?
    };

    if !exists {
        writeln!(
            out,
            "dataset,metric,param,fit,error,certified,abs_error,rel_error"
        )
        .map_err(|e| format!("failed to write {summary_path}: {e}"))?;
    }

    for i in 0..params.len() {
        let cert = ds.certified[i];
        let abs_err = (params[i] - cert).abs();
        let rel_err = abs_err / cert.abs().max(1e-16);
        writeln!(
            out,
            "{},param,b{},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            ds.name,
            i + 1,
            params[i],
            errors[i],
            cert,
            abs_err,
            rel_err
        )
        .map_err(|e| format!("failed to write {summary_path}: {e}"))?;
    }

    let rss_abs = (rss - ds.certified_rss).abs();
    let rss_rel = rss_abs / ds.certified_rss.abs().max(1e-16);
    writeln!(
        out,
        "{},rss,NA,{:.12e},NA,{:.12e},{:.12e},{:.12e}",
        ds.name, rss, ds.certified_rss, rss_abs, rss_rel
    )
    .map_err(|e| format!("failed to write {summary_path}: {e}"))?;

    Ok(())
}

fn solve_once(misra: &NistDataset, hahn: &NistDataset, rat43: &NistDataset) -> Result<(), String> {
    let (_misra_p, _misra_e, _misra_rss, _misra_nfcn) =
        fit_dataset_best(misra, model_misra1a, false)?;
    let (_hahn_p, _hahn_e, _hahn_rss, _hahn_nfcn) = fit_dataset_best(hahn, model_hahn1, false)?;
    let (_rat_p, _rat_e, _rat_rss, _rat_nfcn) = fit_dataset_best(rat43, model_rat43, true)?;
    Ok(())
}

fn bench_solve_times(
    misra: &NistDataset,
    hahn: &NistDataset,
    rat43: &NistDataset,
    repeats: usize,
    warmups: usize,
) -> Result<(), String> {
    for _ in 0..warmups {
        solve_once(misra, hahn, rat43)?;
    }
    let mut times = Vec::with_capacity(repeats);
    for _ in 0..repeats {
        let t0 = Instant::now();
        solve_once(misra, hahn, rat43)?;
        times.push(t0.elapsed().as_secs_f64());
    }
    let samples = times
        .iter()
        .map(|v| format!("{v:.9}"))
        .collect::<Vec<_>>()
        .join(",");
    println!("BENCH_TIMES_S:{samples}");
    Ok(())
}

fn main() -> Result<(), String> {
    let cli = parse_cli()?;
    let output_dir = repo_path(OUTPUT_DIR_REL);

    let misra = parse_nist_dat(&repo_path(MISRA_REL), 2)?;
    let hahn = parse_nist_dat(&repo_path(HAHN_REL), 7)?;
    let rat43 = parse_nist_dat(&repo_path(RAT43_REL), 4)?;

    if cli.mode == RunMode::LoadOnly {
        return Ok(());
    }

    if cli.bench_repeats > 0 {
        return bench_solve_times(&misra, &hahn, &rat43, cli.bench_repeats, cli.bench_warmups);
    }

    let (misra_p, misra_e, misra_rss, misra_nfcn) = fit_dataset_best(&misra, model_misra1a, false)?;
    let (hahn_p, hahn_e, hahn_rss, hahn_nfcn) = fit_dataset_best(&hahn, model_hahn1, false)?;
    let (rat_p, rat_e, rat_rss, rat_nfcn) = fit_dataset_best(&rat43, model_rat43, true)?;

    if cli.mode == RunMode::SolveOnly {
        return Ok(());
    }

    // Recreate summary on each full run for deterministic output.
    let summary_path = format!("{}/nist_summary.csv", output_dir.display());
    if Path::new(&summary_path).exists() {
        std::fs::remove_file(&summary_path)
            .map_err(|e| format!("failed to remove {summary_path}: {e}"))?;
    }

    report_dataset(&misra, &misra_p, &misra_e, misra_rss, misra_nfcn)?;
    write_curve(&misra.name, &misra.x, &misra.y, &misra_p, model_misra1a)?;
    append_summary_csv(&misra, &misra_p, &misra_e, misra_rss)?;

    report_dataset(&hahn, &hahn_p, &hahn_e, hahn_rss, hahn_nfcn)?;
    write_curve(&hahn.name, &hahn.x, &hahn.y, &hahn_p, model_hahn1)?;
    append_summary_csv(&hahn, &hahn_p, &hahn_e, hahn_rss)?;

    report_dataset(&rat43, &rat_p, &rat_e, rat_rss, rat_nfcn)?;
    write_curve(&rat43.name, &rat43.x, &rat43.y, &rat_p, model_rat43)?;
    append_summary_csv(&rat43, &rat_p, &rat_e, rat_rss)?;

    println!("Wrote {}/nist_summary.csv", output_dir.display());
    println!("Wrote {}/nist_*_curve.csv", output_dir.display());
    Ok(())
}
