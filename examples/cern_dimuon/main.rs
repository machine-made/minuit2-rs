//! CERN dimuon mass peak fits (Z-region).
//!
//! Data sources:
//! - examples/data/cern/MuRun2010B_0.csv
//! - examples/data/cern/Zmumu.csv
//!
//! Run:
//!   cargo run --example cern_dimuon

use std::fs::{File, create_dir_all};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::time::Instant;

use minuit2::{FCN, MnHesse, MnMigrad};

const MURUN_PATH: &str = "examples/data/cern/MuRun2010B_0.csv";
const ZMUMU_PATH: &str = "examples/data/cern/Zmumu.csv";
const OUTPUT_MURUN: &str = "examples/cern_dimuon/output/cern_murun2010b0_jpsi_curve.csv";
const OUTPUT_ZMUMU: &str = "examples/cern_dimuon/output/cern_zmumu_zpeak_curve.csv";

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

struct HistChi2 {
    x: Vec<f64>,
    y: Vec<f64>,
    sigma: Vec<f64>,
}

impl HistChi2 {
    fn model(p: &[f64], x: f64) -> f64 {
        let amp = p[0];
        let mu = p[1];
        let sig = p[2];
        let c0 = p[3];
        let c1 = p[4];
        if sig <= 0.05 {
            return f64::NAN;
        }
        let z = (x - mu) / sig;
        let peak = amp * (-0.5 * z * z).exp();
        let bg = c0 + c1 * (x - 91.0);
        (peak + bg).max(1e-9)
    }
}

impl FCN for HistChi2 {
    fn value(&self, p: &[f64]) -> f64 {
        let mut chi2 = 0.0;
        for i in 0..self.x.len() {
            let pred = Self::model(p, self.x[i]);
            if !pred.is_finite() {
                return 1e30;
            }
            let r = (self.y[i] - pred) / self.sigma[i];
            chi2 += r * r;
        }
        chi2
    }

    fn error_def(&self) -> f64 {
        1.0
    }
}

fn parse_csv_header_indices(header: &str) -> Vec<String> {
    header.split(',').map(|h| h.trim().to_string()).collect()
}

fn read_masses_from_column(path: &str, column_name: &str) -> Result<Vec<f64>, String> {
    let f = File::open(path).map_err(|e| format!("failed to open {path}: {e}"))?;
    let mut reader = BufReader::new(f);

    let mut header = String::new();
    reader
        .read_line(&mut header)
        .map_err(|e| format!("failed reading header from {path}: {e}"))?;
    if header.trim().is_empty() {
        return Err(format!("empty header in {path}"));
    }

    let names = parse_csv_header_indices(&header);
    let col_idx = names
        .iter()
        .position(|n| n == column_name)
        .ok_or_else(|| format!("column {column_name} not found in {path}"))?;

    let mut masses = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("error reading {path}: {e}"))?;
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split(',').collect();
        if col_idx >= cols.len() {
            continue;
        }
        if let Ok(v) = cols[col_idx].trim().parse::<f64>()
            && v.is_finite()
        {
            masses.push(v);
        }
    }
    Ok(masses)
}

fn read_zmumu_masses_from_kinematics(path: &str) -> Result<Vec<f64>, String> {
    let f = File::open(path).map_err(|e| format!("failed to open {path}: {e}"))?;
    let mut reader = BufReader::new(f);

    let mut header = String::new();
    reader
        .read_line(&mut header)
        .map_err(|e| format!("failed reading header from {path}: {e}"))?;
    let names = parse_csv_header_indices(&header);

    let i_pt1 = names
        .iter()
        .position(|n| n == "pt1")
        .ok_or_else(|| format!("pt1 not found in {path}"))?;
    let i_eta1 = names
        .iter()
        .position(|n| n == "eta1")
        .ok_or_else(|| format!("eta1 not found in {path}"))?;
    let i_phi1 = names
        .iter()
        .position(|n| n == "phi1")
        .ok_or_else(|| format!("phi1 not found in {path}"))?;
    let i_pt2 = names
        .iter()
        .position(|n| n == "pt2")
        .ok_or_else(|| format!("pt2 not found in {path}"))?;
    let i_eta2 = names
        .iter()
        .position(|n| n == "eta2")
        .ok_or_else(|| format!("eta2 not found in {path}"))?;
    let i_phi2 = names
        .iter()
        .position(|n| n == "phi2")
        .ok_or_else(|| format!("phi2 not found in {path}"))?;

    let mut masses = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("error reading {path}: {e}"))?;
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split(',').collect();
        if cols.len() <= i_phi2 {
            continue;
        }
        let pt1 = match cols[i_pt1].trim().parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let eta1 = match cols[i_eta1].trim().parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let phi1 = match cols[i_phi1].trim().parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let pt2 = match cols[i_pt2].trim().parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let eta2 = match cols[i_eta2].trim().parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let phi2 = match cols[i_phi2].trim().parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };

        // Massless two-body approximation:
        // m^2 = 2 pT1 pT2 (cosh(eta1-eta2) - cos(phi1-phi2))
        let m2 = 2.0 * pt1 * pt2 * ((eta1 - eta2).cosh() - (phi1 - phi2).cos());
        if m2 > 0.0 && m2.is_finite() {
            masses.push(m2.sqrt());
        }
    }
    Ok(masses)
}

fn histogram(masses: &[f64], low: f64, high: f64, bins: usize) -> (Vec<f64>, Vec<f64>) {
    let mut counts = vec![0.0; bins];
    let width = (high - low) / bins as f64;
    for &m in masses {
        if m < low || m >= high {
            continue;
        }
        let mut idx = ((m - low) / width).floor() as usize;
        if idx >= bins {
            idx = bins - 1;
        }
        counts[idx] += 1.0;
    }
    let centers: Vec<f64> = (0..bins).map(|i| low + (i as f64 + 0.5) * width).collect();
    (centers, counts)
}

fn fit_histogram(
    label: &str,
    x: &[f64],
    y: &[f64],
    mu_start: f64,
    sigma_start: f64,
    output_path: &str,
    emit_output: bool,
) -> Result<(Vec<f64>, Vec<f64>), String> {
    let sigma: Vec<f64> = y.iter().map(|&c| c.max(1.0).sqrt()).collect();
    let max_count = y
        .iter()
        .copied()
        .fold(0.0_f64, |a, b| if a > b { a } else { b });
    let mean_bg = y.iter().sum::<f64>() / y.len().max(1) as f64;

    let fcn = HistChi2 {
        x: x.to_vec(),
        y: y.to_vec(),
        sigma: sigma.clone(),
    };

    let min = MnMigrad::new()
        .add_lower_limited("amp", max_count, (max_count * 0.1).max(1.0), 0.0)
        .add("mu", mu_start, (sigma_start * 0.2).max(0.01))
        .add_lower_limited("sigma", sigma_start, (sigma_start * 0.1).max(0.01), 0.05)
        .add_lower_limited("c0", mean_bg.max(1.0), 0.5, 0.0)
        .add("c1", 0.0, 0.05)
        .with_strategy(2)
        .tolerance(0.01)
        .max_fcn(100_000)
        .minimize(&fcn);
    let hesse = MnHesse::new().calculate(&fcn, &min);
    let hs = hesse.user_state();
    let p = hesse.params();
    if !hesse.is_valid() {
        return Err(format!("{label}: fit did not converge to a valid minimum"));
    }

    if emit_output {
        println!("=== {} ===", label);
        println!("bins                : {}", x.len());
        println!("valid               : {}", hesse.is_valid());
        println!("nfcn                : {}", hesse.nfcn());
        println!("chi2                : {:.6}", hesse.fval());
        println!("ndf                 : {}", x.len() as i64 - 5);
        println!(
            "chi2/ndf            : {:.6}",
            hesse.fval() / (x.len() as f64 - 5.0)
        );
        println!(
            "mu (GeV)            : {:.6} +/- {:.6}",
            hs.value("mu").unwrap_or(p[1]),
            hs.error("mu").unwrap_or(f64::NAN)
        );
        println!(
            "sigma (GeV)         : {:.6} +/- {:.6}",
            hs.value("sigma").unwrap_or(p[2]),
            hs.error("sigma").unwrap_or(f64::NAN)
        );
        println!();

        if let Some(parent) = Path::new(output_path).parent() {
            create_dir_all(parent).map_err(|e| format!("failed to create output dir: {e}"))?;
        }
        let mut out = File::create(output_path)
            .map_err(|e| format!("failed to create {}: {e}", output_path))?;
        writeln!(out, "bin_center,count,sigma,model,residual")
            .map_err(|e| format!("failed to write {}: {e}", output_path))?;
        for i in 0..x.len() {
            let m = HistChi2::model(&p, x[i]);
            let r = y[i] - m;
            writeln!(
                out,
                "{:.6},{:.6},{:.6},{:.6},{:.6}",
                x[i], y[i], sigma[i], m, r
            )
            .map_err(|e| format!("failed to write {}: {e}", output_path))?;
        }
    }

    Ok((p, sigma))
}

fn solve_once(
    x_murun: &[f64],
    y_murun: &[f64],
    x_zmumu: &[f64],
    y_zmumu: &[f64],
) -> Result<(), String> {
    let _ = fit_histogram(
        "MuRun2010B_0 (J/psi region)",
        x_murun,
        y_murun,
        3.10,
        0.12,
        OUTPUT_MURUN,
        false,
    )?;
    let _ = fit_histogram(
        "Zmumu (reconstructed mass)",
        x_zmumu,
        y_zmumu,
        91.0,
        2.5,
        OUTPUT_ZMUMU,
        false,
    )?;
    Ok(())
}

fn bench_solve_times(
    x_murun: &[f64],
    y_murun: &[f64],
    x_zmumu: &[f64],
    y_zmumu: &[f64],
    repeats: usize,
    warmups: usize,
) -> Result<(), String> {
    for _ in 0..warmups {
        solve_once(x_murun, y_murun, x_zmumu, y_zmumu)?;
    }
    let mut times = Vec::with_capacity(repeats);
    for _ in 0..repeats {
        let t0 = Instant::now();
        solve_once(x_murun, y_murun, x_zmumu, y_zmumu)?;
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

    // Fit 1: direct invariant mass column from MuRun2010B_0 (J/psi region)
    let masses_murun = read_masses_from_column(MURUN_PATH, "M")?;
    let jpsi_window_murun: Vec<f64> = masses_murun
        .iter()
        .copied()
        .filter(|m| *m >= 2.0 && *m <= 5.0)
        .collect();
    let (x_murun, y_murun) = histogram(&jpsi_window_murun, 2.0, 5.0, 60);

    // Fit 2: reconstructed mass from Zmumu kinematics.
    let masses_zmumu = read_zmumu_masses_from_kinematics(ZMUMU_PATH)?;
    let z_window_zmumu: Vec<f64> = masses_zmumu
        .iter()
        .copied()
        .filter(|m| *m >= 60.0 && *m <= 120.0)
        .collect();
    let (x_zmumu, y_zmumu) = histogram(&z_window_zmumu, 60.0, 120.0, 60);

    if cli.mode == RunMode::LoadOnly {
        return Ok(());
    }

    if cli.bench_repeats > 0 {
        return bench_solve_times(
            &x_murun,
            &y_murun,
            &x_zmumu,
            &y_zmumu,
            cli.bench_repeats,
            cli.bench_warmups,
        );
    }

    let emit_output = cli.mode == RunMode::Full;
    let _ = fit_histogram(
        "MuRun2010B_0 (J/psi region)",
        &x_murun,
        &y_murun,
        3.10,
        0.12,
        OUTPUT_MURUN,
        emit_output,
    )?;
    let _ = fit_histogram(
        "Zmumu (reconstructed mass)",
        &x_zmumu,
        &y_zmumu,
        91.0,
        2.5,
        OUTPUT_ZMUMU,
        emit_output,
    )?;

    if cli.mode == RunMode::Full {
        println!("Wrote {OUTPUT_MURUN}");
        println!("Wrote {OUTPUT_ZMUMU}");
    }
    Ok(())
}
