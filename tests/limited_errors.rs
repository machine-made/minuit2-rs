//! Errors of parameters with limits against ROOT Minuit2 v6-36-08.
//!
//! The transform from internal to external values is not linear, so ROOT takes the external
//! error of a parameter with limits as the mean of two one-sided changes of its value.

use std::f64::consts::PI;
use std::path::PathBuf;
use std::process::Command;

use minuit2::MnUserParameters;
use serde_json::Value;

/// Limits of one parameter, as `ref_compare_runner` sets them for each workload.
#[derive(Clone, Copy)]
enum Limits {
    Lower(f64),
    Upper(f64),
    /// Both limits. The sine transform reaches every external value at two internal angles:
    /// `theta` on the rising side of the sine and `pi - theta` on the falling side.
    /// `falling` tells which side the fit ended on.
    Both {
        lower: f64,
        upper: f64,
        falling: bool,
    },
}

const STATE_CASES: [(&str, [Limits; 2]); 3] = [
    (
        "quadratic2_lower_limited_migrad",
        [Limits::Lower(0.0), Limits::Lower(-2.5)],
    ),
    (
        "quadratic2_upper_limited_migrad",
        [Limits::Upper(1.8), Limits::Upper(-1.5)],
    ),
    (
        "quadratic2_limited_migrad",
        [
            Limits::Both {
                lower: 0.0,
                upper: 2.0,
                falling: true,
            },
            Limits::Both {
                lower: -3.0,
                upper: -1.0,
                falling: false,
            },
        ],
    ),
];

fn numbers(value: &Value, what: &str) -> Vec<f64> {
    value
        .as_array()
        .unwrap_or_else(|| panic!("{what}: not an array"))
        .iter()
        .map(|v| {
            v.as_f64()
                .unwrap_or_else(|| panic!("{what}: entry is not a number"))
        })
        .collect()
}

fn rust_report(id: &str) -> Value {
    let output = Command::new(env!("CARGO_BIN_EXE_ref_compare_runner"))
        .args(["--workload", id])
        .output()
        .expect("run ref_compare_runner");
    assert!(output.status.success(), "runner failed for {id}");
    serde_json::from_slice(&output.stdout).expect("runner output is JSON")
}

fn root_report(id: &str) -> Value {
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "reports",
        "verification",
        "raw",
        "ref",
        &format!("{id}.json"),
    ]
    .iter()
    .collect();
    let text = std::fs::read_to_string(&path).expect("read the stored ROOT file");
    serde_json::from_str(&text).expect("stored ROOT file is JSON")
}

/// Print one error, and record it when it is further from ROOT than `tolerance` (relative).
fn compare(label: &str, got: f64, want: f64, tolerance: f64, mismatches: &mut Vec<String>) {
    let rel = (got - want).abs() / want.abs();
    println!("{label}: rust {got} root {want} relative difference {rel:e}");
    if rel.is_nan() || rel > tolerance {
        mismatches.push(format!("{label}: rust {got} root {want}"));
    }
}

fn assert_no_mismatch(mismatches: &[String], tolerance: f64) {
    assert!(
        mismatches.is_empty(),
        "errors differ from ROOT by more than {tolerance:e} (relative):\n{}",
        mismatches.join("\n")
    );
}

/// ROOT's stored fit state gives ROOT's stored errors through `int2ext_error`.
///
/// The stored file holds the external values and the external covariance. The internal error
/// follows from them: `sqrt(cov_ii)` is the external error scaled by the derivative of the
/// transform, so dividing by that derivative gives the internal error that ROOT passed in.
///
/// For two limits, the branch of the internal angle matters. When the internal error is above
/// one, ROOT replaces the change at `internal + err` by `upper - lower`. That is the change
/// on the rising side of the sine on one branch and on the falling side on the other, so the
/// two branches give different errors. The stored file keeps only the external value, so the
/// table above records the branch the fit ended on. Only `x` of `quadratic2_limited_migrad`
/// has an internal error above one, and it ended on the falling side.
///
/// The worst relative difference measured is 5e-16. The tolerance leaves a margin of about
/// 2000. A wrong formula moves an error by tens of percent.
#[test]
fn limited_error_formula_matches_root_state() {
    const TOLERANCE: f64 = 1e-12;
    let mut mismatches = Vec::new();
    for (id, limits) in STATE_CASES {
        let root = root_report(id);
        let values = numbers(&root["params"], id);
        let errors = numbers(&root["errors"], id);
        assert_eq!(errors.len(), limits.len(), "{id}: number of errors");

        let mut params = MnUserParameters::new();
        for (k, limit) in limits.iter().enumerate() {
            let name = format!("p{k}");
            match *limit {
                Limits::Lower(lower) => params.add_lower_limited(name, values[k], 0.1, lower),
                Limits::Upper(upper) => params.add_upper_limited(name, values[k], 0.1, upper),
                Limits::Both { lower, upper, .. } => {
                    params.add_limited(name, values[k], 0.1, lower, upper)
                }
            };
        }
        let trafo = params.trafo();

        for (k, want) in errors.iter().enumerate() {
            let mut internal = trafo.ext2int(k, values[k]);
            if matches!(limits[k], Limits::Both { falling: true, .. }) {
                internal = PI - internal;
            }
            let covariance = root["covariance"][k][k]
                .as_f64()
                .unwrap_or_else(|| panic!("{id}: covariance entry is not a number"));
            let internal_error = covariance.sqrt() / trafo.dint2ext(k, internal).abs();
            let got = trafo.int2ext_error(k, internal, internal_error);
            compare(
                &format!("{id} errors[{k}]"),
                got,
                *want,
                TOLERANCE,
                &mut mismatches,
            );
        }
    }
    assert_no_mismatch(&mismatches, TOLERANCE);
}

/// The Rust fit reports the errors that ROOT reports, for the workloads where the two fits
/// end at the same point.
///
/// The fit of `quadratic2_upper_limited_migrad` ends at a different point than ROOT's: it
/// takes 43 function calls where ROOT takes 61, and its covariance differs by about 3e-4. Its
/// errors are therefore not compared end to end here. The test above covers its
/// error formula.
///
/// The worst relative difference measured is 8e-11, against a tolerance of 1e-6.
#[test]
fn limited_errors_match_root_end_to_end() {
    const TOLERANCE: f64 = 1e-6;
    let mut mismatches = Vec::new();
    for id in [
        "quadratic2_lower_limited_migrad",
        "quadratic2_limited_migrad",
    ] {
        let rust = numbers(&rust_report(id)["errors"], id);
        let root = numbers(&root_report(id)["errors"], id);
        assert_eq!(rust.len(), root.len(), "{id}: number of errors differs");
        for (k, (r, o)) in rust.iter().zip(&root).enumerate() {
            compare(
                &format!("{id} errors[{k}]"),
                *r,
                *o,
                TOLERANCE,
                &mut mismatches,
            );
        }
    }
    assert_no_mismatch(&mismatches, TOLERANCE);
}

fn assert_close(got: f64, want: f64, what: &str) {
    assert!(
        (got - want).abs() <= 1e-15,
        "{what}: got {got}, want {want}"
    );
}

#[test]
fn both_limits_error_is_sine_of_internal_error() {
    let mut params = MnUserParameters::new();
    params.add_limited("x", 0.0, 0.1, -1.0, 1.0);
    let got = params.trafo().int2ext_error(0, 0.0, 0.1);
    assert_close(got, 0.1_f64.sin(), "limits -1 and 1");
    assert_close(got, 0.099_833_416_646_828_15, "limits -1 and 1, literal");
}

#[test]
fn both_limits_error_above_one_uses_the_whole_range() {
    let mut params = MnUserParameters::new();
    params.add_limited("x", 0.0, 0.1, -1.0, 1.0);
    let got = params.trafo().int2ext_error(0, 0.0, 2.0);
    // ROOT sets du1 to upper - lower = 2; du2 is sin(-2).
    assert_close(
        got,
        2.0_f64.midpoint(2.0_f64.sin()),
        "limits -1 and 1, error 2",
    );
}

#[test]
fn lower_limit_error_averages_the_two_one_sided_changes() {
    let mut params = MnUserParameters::new();
    params.add_lower_limited("x", 1.0, 0.1, 0.0);
    let got = params.trafo().int2ext_error(0, 1.0, 0.5);
    assert_close(got, 0.342_370_824_491_049_83, "lower limit 0");
}

#[test]
fn upper_limit_error_averages_the_two_one_sided_changes() {
    let mut params = MnUserParameters::new();
    params.add_upper_limited("x", 1.0, 0.1, 0.0);
    let got = params.trafo().int2ext_error(0, 1.0, 0.5);
    assert_close(got, 0.342_370_824_491_049_83, "upper limit 0");
}

#[test]
fn error_without_limits_is_unchanged() {
    let mut params = MnUserParameters::new();
    params.add("x", 1.0, 0.1);
    assert_eq!(
        params.trafo().int2ext_error(0, 1.0, 0.5).to_bits(),
        0.5_f64.to_bits()
    );
}
