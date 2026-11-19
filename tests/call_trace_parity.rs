use std::cell::RefCell;

use minuit2::{FCN, MnMigrad};

#[derive(Default)]
struct RecordedRosenbrock {
    calls: RefCell<Vec<Vec<f64>>>,
}

impl FCN for RecordedRosenbrock {
    fn value(&self, p: &[f64]) -> f64 {
        self.calls.borrow_mut().push(p.to_vec());
        let t1 = p[1] - p[0] * p[0];
        let t2 = 1.0 - p[0];
        100.0 * t1 * t1 + t2 * t2
    }
}

#[test]
#[expect(
    clippy::cast_possible_truncation,
    reason = "the reference call count is 140 and fits in usize"
)]
fn rosenbrock2_calls_equal_root_calls() {
    let root: serde_json::Value = serde_json::from_str(include_str!(
        "../reports/verification/raw/trace/nfcn_accounting_rosenbrock2_maxfcn140.root.json"
    ))
    .unwrap();
    let reference = &root["result"];
    let trace = reference["trace"].as_array().unwrap();
    let fcn = RecordedRosenbrock::default();
    let minimum = MnMigrad::new()
        .add("x", 0.0, 0.1)
        .add("y", 0.0, 0.1)
        .max_fcn(140)
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
        // The bound is relative to the larger coordinate. The measured largest difference
        // over the 140 calls is 2.5e-10 (call 104, second coordinate): both sides run the same
        // steps and differ in the last digits of the floating-point path. 1e-9 leaves a margin
        // of four. A wrong step moves a coordinate by far more, as the failure at call 4 showed
        // (relative difference 2.7e-7).
        assert!(
            actual.len() == vector.len()
                && actual
                    .iter()
                    .zip(&vector)
                    .all(|(a, b)| { (a - b).abs() <= 1e-9 * a.abs().max(b.abs()) }),
            "first differing call {}: Rust {actual:?}, ROOT {vector:?}, phase {}",
            expected["call"],
            expected["phase"],
        );
    }
    assert_eq!(trace.len(), 140);
    assert_eq!(calls.len(), trace.len(), "call sequence length differs");
    assert_eq!(minimum.nfcn(), calls.len());
    assert_eq!(
        minimum.nfcn(),
        reference["reported_nfcn"].as_u64().unwrap() as usize
    );
    assert_eq!(minimum.is_valid(), reference["valid"].as_bool().unwrap());
}
