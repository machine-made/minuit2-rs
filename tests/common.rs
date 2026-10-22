use std::fmt::Display;

/// Shared assertion for the FunctionMinimum `Display` contract.
pub fn assert_function_minimum_display<T: Display>(result: &T, required_params: &[&str]) {
    let output = format!("{result}");

    assert!(output.contains("FunctionMinimum:"));
    assert!(output.contains("valid:"));
    assert!(output.contains("parameters:"));

    let fval_line = output
        .lines()
        .find(|line| line.trim_start().starts_with("fval:"))
        .expect("fval line should exist");
    let fval_text = fval_line
        .split_whitespace()
        .next_back()
        .expect("fval value should be present");
    let fval = fval_text.parse::<f64>().expect("fval should parse as f64");
    assert!(fval.is_finite());

    for p in required_params {
        assert!(output.contains(p));
    }
}
