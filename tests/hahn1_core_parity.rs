use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use minuit2::{FCN, MnMigrad};

struct Hahn1Dataset {
    x: Vec<f64>,
    y: Vec<f64>,
    start2: Vec<f64>,
    certified: Vec<f64>,
}

struct Hahn1Fcn {
    x: Vec<f64>,
    y: Vec<f64>,
}

impl FCN for Hahn1Fcn {
    fn value(&self, p: &[f64]) -> f64 {
        let mut rss = 0.0;
        for (&x, &y) in self.x.iter().zip(self.y.iter()) {
            let x2 = x * x;
            let x3 = x2 * x;
            let den = 1.0 + p[4] * x + p[5] * x2 + p[6] * x3;
            if den.abs() < 1e-300 {
                return 1e30;
            }
            let pred = (p[0] + p[1] * x + p[2] * x2 + p[3] * x3) / den;
            if !pred.is_finite() {
                return 1e30;
            }
            let r = y - pred;
            rss += r * r;
        }
        rss
    }

    fn error_def(&self) -> f64 {
        1.0
    }
}

fn repo_path(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel)
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

fn load_hahn1(path: &Path) -> Hahn1Dataset {
    let f = File::open(path).unwrap_or_else(|e| panic!("failed to open {}: {e}", path.display()));
    let reader = BufReader::new(f);
    let mut start2 = Vec::new();
    let mut certified = Vec::new();
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut in_data = false;

    for line in reader.lines() {
        let line = line.unwrap_or_else(|e| panic!("error reading {}: {e}", path.display()));
        let s = line.trim();
        if let Some((lhs, rhs)) = s.split_once('=') {
            let lhs = lhs.trim_start();
            if lhs.starts_with('b') && lhs[1..].chars().next().is_some_and(|c| c.is_ascii_digit()) {
                let nums = parse_floats(rhs);
                if nums.len() >= 4 {
                    start2.push(nums[1]);
                    certified.push(nums[2]);
                }
            }
        }
        if let Some(tail) = s.strip_prefix("Data:")
            && tail.trim_start().starts_with('y')
        {
            in_data = true;
            continue;
        }
        if in_data {
            let nums = parse_floats(s);
            if nums.len() >= 2 {
                y.push(nums[0]);
                x.push(nums[1]);
            }
        }
    }

    assert_eq!(start2.len(), 7, "failed to parse Hahn1 Start 2");
    assert_eq!(
        certified.len(),
        7,
        "failed to parse Hahn1 certified parameters"
    );
    assert!(!x.is_empty(), "failed to parse Hahn1 observations");
    Hahn1Dataset {
        x,
        y,
        start2,
        certified,
    }
}

#[test]
fn hahn1_core_start2_reaches_certified_solution() {
    let ds = load_hahn1(&repo_path("examples/data/nist/Hahn1.dat"));
    let fcn = Hahn1Fcn {
        x: ds.x.clone(),
        y: ds.y.clone(),
    };

    let mut migrad = MnMigrad::new();
    for (i, &start) in ds.start2.iter().enumerate() {
        migrad = migrad.add(format!("b{}", i + 1), start, 0.1);
    }

    let min = migrad
        .with_strategy(1)
        .tolerance(0.1)
        .max_fcn(1_000_000)
        .minimize(&fcn);

    assert!(
        min.is_valid(),
        "Hahn1 core fit should be valid: fval={} edm={} nfcn={}",
        min.fval(),
        min.edm(),
        min.nfcn()
    );
    assert!(
        (min.fval() - 1.5324382854).abs() <= 1.5324382854e-2,
        "Hahn1 fval {} did not reach the certified RSS basin",
        min.fval()
    );

    let params = min.params();
    for (i, (&got, &cert)) in params.iter().zip(ds.certified.iter()).enumerate() {
        let rel = (got - cert).abs() / cert.abs().max(1e-300);
        assert!(
            rel <= 1e-2,
            "b{} fitted {:.10e} vs certified {:.10e}, relative error {:.3e}",
            i + 1,
            got,
            cert,
            rel
        );
    }
}
