//! Checks the upper limited quadratic objective call sequence against ROOT.
use std::cell::RefCell;

use minuit2::{FCN, MnMigrad};

#[derive(Default)]
struct RecordedQuadratic {
    calls: RefCell<Vec<Vec<f64>>>,
}

impl FCN for RecordedQuadratic {
    fn value(&self, p: &[f64]) -> f64 {
        self.calls.borrow_mut().push(p.to_vec());
        let (x, y) = (p[0], p[1]);
        let dx = x - 1.0;
        let dy = y + 2.0;
        dx * dx + 4.0 * dy * dy + 0.3 * x * y
    }
}

#[test]
fn upper_limited_calls_equal_root_calls() {
    let root: serde_json::Value = serde_json::from_str(include_str!(
        "../reports/verification/raw/trace/nfcn_accounting_quadratic2_upper_limited_migrad.root.json"
    ))
    .unwrap();
    let reference = &root["result"];
    let trace = reference["trace"].as_array().unwrap();
    let fcn = RecordedQuadratic::default();
    let minimum = MnMigrad::new()
        .add_upper_limited("x", 0.4, 0.1, 1.8)
        .add_upper_limited("y", -1.0, 0.1, -1.5)
        .tolerance(0.1)
        .minimize(&fcn);
    let calls = fcn.calls.borrow();

    for (actual, expected) in calls.iter().zip(trace) {
        let vector: Vec<f64> = expected["p"]
            .as_array()
            .unwrap()
            .iter()
            .map(|value| value.as_f64().unwrap())
            .collect();
        // The shared steps agree within relative 1e-9 despite transform rounding.
        assert!(
            actual.len() == vector.len()
                && actual
                    .iter()
                    .zip(&vector)
                    .all(|(a, b)| (a - b).abs() <= 1e-9 * a.abs().max(b.abs())),
            "first differing call {}: Rust {actual:?}, ROOT {vector:?}, phase {}",
            expected["call"],
            expected["phase"],
        );
    }
    assert_eq!(trace.len(), 61);
    assert_eq!(calls.len(), trace.len(), "call sequence length differs");
    assert_eq!(minimum.nfcn(), calls.len());
    assert_eq!(minimum.is_valid(), reference["valid"].as_bool().unwrap());
}
