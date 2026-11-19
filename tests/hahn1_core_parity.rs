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

// ROOT Minuit2 v6-36-08 result for this fit: strategy 1, tolerance 0.1, error 0.1, start 2.
// Source file: reports/verification/raw/ref/hahn1_core_start2_migrad_root.txt.
// The program that made it is tools/hahn1_probe/hahn1_probe.cpp.
const ROOT_FVAL: f64 = 1.532_556_828;
const ROOT_NFCN: usize = 456;
const ROOT_PARAMS: [f64; 7] = [
    1.092_152_452_1e+00,
    -1.236_442_692_8e-01,
    4.102_302_183_8e-03,
    -1.450_205_237_2e-06,
    -5.767_145_150_0e-03,
    2.412_354_167_9e-04,
    -1.242_833_394_9e-07,
];

#[test]
fn hahn1_core_start2_matches_root_and_stays_near_certified_solution() {
    let ds = load_hahn1(&repo_path("examples/data/nist/Hahn1.dat"));
    let fcn = Hahn1Fcn {
        x: ds.x.clone(),
        y: ds.y.clone(),
    };

    let mut migrad = MnMigrad::new();
    for (i, &start) in ds.start2.iter().enumerate() {
        migrad = migrad.add(format!("b{}", i + 1), start, 0.1);
    }

    let calls = std::cell::Cell::new(0_usize);
    let min = migrad
        .with_strategy(1)
        .tolerance(0.1)
        .max_fcn(1_000_000)
        .minimize(&|p: &[f64]| {
            calls.set(calls.get() + 1);
            fcn.value(p)
        });

    assert_eq!(min.nfcn(), calls.get());

    assert!(
        min.is_valid(),
        "Hahn1 core fit should be valid: fval={} edm={} nfcn={}",
        min.fval(),
        min.edm(),
        min.nfcn()
    );
    // The residual sum of squares (RSS) is 1.5325568, which is 1.2e-4 above the certified NIST
    // value 1.5324383, as in ROOT. The fit stops at an EDM of 1.4e-4, under the Migrad goal
    // 0.1 * 0.002 = 2e-4. The bound of 1 percent is 129 times the excess; a fit that stops early
    // ends far higher (the start RSS is 2.1e6).
    assert!(
        (min.fval() - 1.532_438_285_4).abs() <= 1.532_438_285_4e-2,
        "Hahn1 fval {} did not reach the certified RSS basin",
        min.fval()
    );

    let params = min.params();

    // ROOT v6-36-08 itself ends at 1.347e-2 (b1) and 1.678e-2 (b4) from the certified values.
    // The bound is ROOT's largest value rounded up.
    for (i, (&got, &cert)) in params.iter().zip(ds.certified.iter()).enumerate() {
        let rel = (got - cert).abs() / cert.abs().max(1e-300);
        assert!(
            rel <= 2e-2,
            "b{} fitted {:.10e} vs certified {:.10e}, relative error {:.3e}",
            i + 1,
            got,
            cert,
            rel
        );
    }

    // Measured port difference is 2.8e-8, and ROOT prints 10 digits.
    let fval_rel = (min.fval() - ROOT_FVAL).abs() / ROOT_FVAL;
    assert!(
        fval_rel <= 1e-6,
        "Hahn1 fval {} differs from ROOT {} by {:.3e}",
        min.fval(),
        ROOT_FVAL,
        fval_rel
    );

    // Measured largest difference is 2.7e-6, and the fit stops at an EDM of 1e-4.
    // The last digits of a parameter are not determined at that EDM.
    // The bound is 30 times the measured value and far below the 1e-2 distance to the certified values.
    for (i, (&got, &root)) in params.iter().zip(ROOT_PARAMS.iter()).enumerate() {
        let rel = (got - root).abs() / root.abs();
        assert!(
            rel <= 1e-4,
            "b{} fitted {:.10e} vs ROOT {:.10e}, relative difference {:.3e}",
            i + 1,
            got,
            root,
            rel
        );
    }

    // Both start Hesse after 294 calls and use fourteen gradient refinement calls.
    // Both continuations take 112 calls, reaching the same total.
    assert_eq!(min.nfcn(), ROOT_NFCN);
}
