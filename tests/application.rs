use minuit2::application::{DEFAULT_TOLERANCE, default_max_fcn};

#[test]
fn default_max_fcn_matches_formula() {
    assert_eq!(default_max_fcn(0), 200);
    assert_eq!(default_max_fcn(1), 305);
    assert_eq!(default_max_fcn(2), 420);
    assert_eq!(default_max_fcn(5), 825);
}

#[test]
fn default_tolerance_is_root_compatible() {
    assert!((DEFAULT_TOLERANCE - 0.1).abs() < 1e-15);
}
