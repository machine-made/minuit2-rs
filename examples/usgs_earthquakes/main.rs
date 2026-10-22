//! USGS earthquake catalog: Gutenberg-Richter fit.
//!
//! Data source:
//! - examples/data/usgs/earthquakes_2025_m4p5.csv
//!
//! Run:
//!   cargo run --example usgs_earthquakes

use std::fs::{File, create_dir_all};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::time::Instant;

use minuit2::{FCN, MnHesse, MnMigrad};

const DATA_PATH: &str = "examples/data/usgs/earthquakes_2025_m4p5.csv";
const OUTPUT_PATH: &str = "examples/usgs_earthquakes/output/usgs_gutenberg_richter_curve.csv";

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

struct GRChi2 {
    m: Vec<f64>,
    log10_n: Vec<f64>,
    sigma_log10_n: Vec<f64>,
}

impl FCN for GRChi2 {
    fn value(&self, p: &[f64]) -> f64 {
        // log10 N(M>=m) = a - b m
        let a = p[0];
        let b = p[1];
        let mut chi2 = 0.0;
        for i in 0..self.m.len() {
            let pred = a - b * self.m[i];
            if !pred.is_finite() {
                return 1e30;
            }
            let r = (self.log10_n[i] - pred) / self.sigma_log10_n[i];
            chi2 += r * r;
        }
        chi2
    }

    fn error_def(&self) -> f64 {
        1.0
    }
}

fn parse_magnitudes(path: &str) -> Result<Vec<f64>, String> {
    let f = File::open(path).map_err(|e| format!("failed to open {path}: {e}"))?;
    let reader = BufReader::new(f);
    let mut mags = Vec::new();

    for (idx, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| format!("error reading {path}: {e}"))?;
        if idx == 0 || line.trim().is_empty() {
            continue;
        }
        // Mag is the 5th field in the USGS CSV schema.
        let cols: Vec<&str> = line.split(',').collect();
        if cols.len() < 5 {
            continue;
        }
        if let Ok(m) = cols[4].trim().parse::<f64>()
            && m.is_finite()
        {
            mags.push(m);
        }
    }

    if mags.is_empty() {
        return Err(format!("no magnitudes parsed from {path}"));
    }
    Ok(mags)
}

fn build_cumulative_curve(mags: &[f64], m_min: f64, m_max: f64, dm: f64) -> (Vec<f64>, Vec<f64>) {
    let mut thresholds = Vec::new();
    let mut counts = Vec::new();
    let mut m = m_min;
    while m <= m_max + 1e-12 {
        let n = mags.iter().filter(|&&x| x >= m).count() as f64;
        if n > 0.0 {
            thresholds.push(m);
            counts.push(n);
        }
        m += dm;
    }
    (thresholds, counts)
}

fn solve_once(fcn: &GRChi2) -> bool {
    let min = MnMigrad::new()
        .add("a", 5.0, 0.1)
        .add("b", 1.0, 0.05)
        .with_strategy(2)
        .tolerance(0.01)
        .max_fcn(20_000)
        .minimize(fcn);
    let hesse = MnHesse::new().calculate(fcn, &min);
    hesse.is_valid()
}

fn bench_solve_times(fcn: &GRChi2, repeats: usize, warmups: usize) -> Result<(), String> {
    for _ in 0..warmups {
        if !solve_once(fcn) {
            return Err("warmup fit failed".to_string());
        }
    }
    let mut times = Vec::with_capacity(repeats);
    for _ in 0..repeats {
        let t0 = Instant::now();
        if !solve_once(fcn) {
            return Err("measured fit failed".to_string());
        }
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
    let mags = parse_magnitudes(DATA_PATH)?;
    let m_min = 4.5;
    let m_max = mags
        .iter()
        .fold(f64::NEG_INFINITY, |a, &b| a.max(b))
        .floor();
    let (m_values, counts) = build_cumulative_curve(&mags, m_min, m_max, 0.1);

    if cli.mode == RunMode::LoadOnly {
        return Ok(());
    }

    let log10_n: Vec<f64> = counts.iter().map(|&n| n.log10()).collect();
    let sigma_log10_n: Vec<f64> = counts
        .iter()
        .map(|&n| 1.0 / (std::f64::consts::LN_10 * n.sqrt()))
        .collect();

    let fcn = GRChi2 {
        m: m_values.clone(),
        log10_n: log10_n.clone(),
        sigma_log10_n: sigma_log10_n.clone(),
    };

    if cli.bench_repeats > 0 {
        return bench_solve_times(&fcn, cli.bench_repeats, cli.bench_warmups);
    }

    let min = MnMigrad::new()
        .add("a", 5.0, 0.1)
        .add("b", 1.0, 0.05)
        .with_strategy(2)
        .tolerance(0.01)
        .max_fcn(20_000)
        .minimize(&fcn);
    let hesse = MnHesse::new().calculate(&fcn, &min);
    let hs = hesse.user_state();

    if cli.mode == RunMode::SolveOnly {
        return Ok(());
    }

    let a = hs
        .value("a")
        .ok_or_else(|| "missing parameter a".to_string())?;
    let b = hs
        .value("b")
        .ok_or_else(|| "missing parameter b".to_string())?;
    let a_err = hs.error("a").ok_or_else(|| "missing error a".to_string())?;
    let b_err = hs.error("b").ok_or_else(|| "missing error b".to_string())?;

    println!("=== USGS Gutenberg-Richter Fit ===");
    println!("events              : {}", mags.len());
    println!("fit points          : {}", m_values.len());
    println!("valid               : {}", hesse.is_valid());
    println!("nfcn                : {}", hesse.nfcn());
    println!("chi2                : {:.6}", hesse.fval());
    println!("ndf                 : {}", m_values.len() as i64 - 2);
    println!(
        "chi2/ndf            : {:.6}",
        hesse.fval() / (m_values.len() as f64 - 2.0)
    );
    println!("a                   : {:.8} +/- {:.6}", a, a_err);
    println!("b-value             : {:.8} +/- {:.6}", b, b_err);

    if let Some(parent) = Path::new(OUTPUT_PATH).parent() {
        create_dir_all(parent).map_err(|e| format!("failed to create output dir: {e}"))?;
    }
    let mut out =
        File::create(OUTPUT_PATH).map_err(|e| format!("failed to create {}: {e}", OUTPUT_PATH))?;
    writeln!(
        out,
        "magnitude_threshold,cumulative_count,log10_count,sigma_log10,pred_log10,residual"
    )
    .map_err(|e| format!("failed to write {}: {e}", OUTPUT_PATH))?;

    for i in 0..m_values.len() {
        let pred = a - b * m_values[i];
        let residual = log10_n[i] - pred;
        writeln!(
            out,
            "{:.3},{:.0},{:.8},{:.8},{:.8},{:.8}",
            m_values[i], counts[i], log10_n[i], sigma_log10_n[i], pred, residual
        )
        .map_err(|e| format!("failed to write {}: {e}", OUTPUT_PATH))?;
    }

    println!("Wrote {OUTPUT_PATH}");
    Ok(())
}
