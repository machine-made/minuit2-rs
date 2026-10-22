//! NOAA Mauna Loa CO2 harmonic trend fit.
//!
//! Data source:
//! - examples/data/noaa/co2_mm_mlo.csv
//!
//! Run:
//!   cargo run --example noaa_co2

use std::f64::consts::PI;
use std::fs::{File, create_dir_all};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::time::Instant;

use minuit2::{FCN, MnHesse, MnMigrad};

const DATA_PATH: &str = "examples/data/noaa/co2_mm_mlo.csv";
const OUTPUT_PATH: &str = "examples/noaa_co2/output/noaa_co2_curve.csv";

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

struct NoaaCo2Chi2 {
    t_years: Vec<f64>,
    y: Vec<f64>,
    sigma: Vec<f64>,
}

impl NoaaCo2Chi2 {
    fn model(p: &[f64], t: f64) -> f64 {
        let w1 = 2.0 * PI * t;
        let w2 = 4.0 * PI * t;
        p[0] + p[1] * t
            + p[2] * t * t
            + p[3] * w1.sin()
            + p[4] * w1.cos()
            + p[5] * w2.sin()
            + p[6] * w2.cos()
            + p[7] * t * w1.sin()
    }
}

impl FCN for NoaaCo2Chi2 {
    fn value(&self, p: &[f64]) -> f64 {
        let mut chi2 = 0.0;
        for i in 0..self.t_years.len() {
            let pred = Self::model(p, self.t_years[i]);
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

type NoaaSeries = (Vec<f64>, Vec<f64>, Vec<f64>);

fn parse_noaa_csv(path: &str) -> Result<NoaaSeries, String> {
    let f = File::open(path).map_err(|e| format!("failed to open {path}: {e}"))?;
    let reader = BufReader::new(f);

    let mut t = Vec::new();
    let mut y = Vec::new();
    let mut sigma = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| format!("error reading {path}: {e}"))?;
        let s = line.trim();
        if s.is_empty() || s.starts_with('#') || s.starts_with("year,month,decimal date") {
            continue;
        }

        let cols: Vec<&str> = s.split(',').collect();
        if cols.len() < 8 {
            continue;
        }

        let decimal_date = match cols[2].trim().parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let avg = match cols[3].trim().parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let unc = match cols[7].trim().parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };

        if unc <= 0.0 || !decimal_date.is_finite() || !avg.is_finite() {
            continue;
        }

        t.push(decimal_date);
        y.push(avg);
        sigma.push(unc.max(1e-6));
    }

    if t.is_empty() {
        return Err(format!("no usable rows parsed from {path}"));
    }

    Ok((t, y, sigma))
}

fn solve_once(fcn: &NoaaCo2Chi2, y0: f64) -> bool {
    let min = MnMigrad::new()
        .add("a0", y0, 0.5)
        .add("a1", 2.0, 0.2)
        .add("a2", 0.0, 0.01)
        .add("b1", 2.0, 0.2)
        .add("c1", 0.0, 0.2)
        .add("b2", 0.5, 0.1)
        .add("c2", 0.0, 0.1)
        .add("d1", 0.0, 0.01)
        .with_strategy(2)
        .tolerance(0.05)
        .max_fcn(100_000)
        .minimize(fcn);

    let hesse = MnHesse::new().calculate(fcn, &min);
    hesse.is_valid()
}

fn bench_solve_times(
    fcn: &NoaaCo2Chi2,
    y0: f64,
    repeats: usize,
    warmups: usize,
) -> Result<(), String> {
    for _ in 0..warmups {
        if !solve_once(fcn, y0) {
            return Err("warmup fit failed".to_string());
        }
    }
    let mut times = Vec::with_capacity(repeats);
    for _ in 0..repeats {
        let t0 = Instant::now();
        if !solve_once(fcn, y0) {
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
    let (t_abs, y, sigma) = parse_noaa_csv(DATA_PATH)?;
    let t0 = t_abs[0];
    let t_years: Vec<f64> = t_abs.iter().map(|&v| v - t0).collect();

    let fcn = NoaaCo2Chi2 {
        t_years: t_years.clone(),
        y: y.clone(),
        sigma: sigma.clone(),
    };

    if cli.bench_repeats > 0 {
        return bench_solve_times(&fcn, y[0], cli.bench_repeats, cli.bench_warmups);
    }

    if cli.mode == RunMode::LoadOnly {
        return Ok(());
    }

    // 8-parameter trend + harmonic model.
    let min = MnMigrad::new()
        .add("a0", y[0], 0.5)
        .add("a1", 2.0, 0.2)
        .add("a2", 0.0, 0.01)
        .add("b1", 2.0, 0.2)
        .add("c1", 0.0, 0.2)
        .add("b2", 0.5, 0.1)
        .add("c2", 0.0, 0.1)
        .add("d1", 0.0, 0.01)
        .with_strategy(2)
        .tolerance(0.05)
        .max_fcn(100_000)
        .minimize(&fcn);

    let hesse = MnHesse::new().calculate(&fcn, &min);
    let hs = hesse.user_state();

    if cli.mode == RunMode::SolveOnly {
        return Ok(());
    }

    println!("=== NOAA CO2 Harmonic Fit ===");
    println!("points      : {}", y.len());
    println!("valid       : {}", hesse.is_valid());
    println!("nfcn        : {}", hesse.nfcn());
    println!("chi2        : {:.6}", hesse.fval());
    println!("ndf         : {}", y.len() as i64 - 8);
    println!("chi2/ndf    : {:.6}", hesse.fval() / (y.len() as f64 - 8.0));
    println!();

    for name in ["a0", "a1", "a2", "b1", "c1", "b2", "c2", "d1"] {
        let v = hs
            .value(name)
            .ok_or_else(|| format!("missing fitted value for {name}"))?;
        let e = hs
            .error(name)
            .ok_or_else(|| format!("missing fitted error for {name}"))?;
        println!("{name:>2} = {v:>14.8} +/- {e:.6}");
    }

    let out_path = Path::new(OUTPUT_PATH);
    if let Some(parent) = out_path.parent() {
        create_dir_all(parent).map_err(|e| format!("failed to create output dir: {e}"))?;
    }

    let mut out =
        File::create(out_path).map_err(|e| format!("failed to create {}: {e}", OUTPUT_PATH))?;
    writeln!(out, "decimal_date,observed,uncertainty,fitted,residual")
        .map_err(|e| format!("failed to write {}: {e}", OUTPUT_PATH))?;

    let params = hesse.params();
    for i in 0..t_abs.len() {
        let fitted = NoaaCo2Chi2::model(&params, t_years[i]);
        let residual = y[i] - fitted;
        writeln!(
            out,
            "{:.6},{:.6},{:.6},{:.6},{:.6}",
            t_abs[i], y[i], sigma[i], fitted, residual
        )
        .map_err(|e| format!("failed to write {}: {e}", OUTPUT_PATH))?;
    }

    println!();
    println!("Wrote {OUTPUT_PATH}");
    Ok(())
}
