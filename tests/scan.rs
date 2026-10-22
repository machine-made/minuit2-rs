use minuit2::{MnMigrad, MnScan};

/// 1D scan of a quadratic: should produce parabolic profile.
#[test]
fn scan_quadratic_profile() {
    let result = MnMigrad::new()
        .add("x", 5.0, 1.0)
        .add("y", -3.0, 1.0)
        .minimize(&|p: &[f64]| p[0] * p[0] + p[1] * p[1]);

    assert!(result.is_valid());

    let scan = MnScan::new(&|p: &[f64]| p[0] * p[0] + p[1] * p[1], &result);
    let points = scan.scan(0, 20, -2.0, 2.0);

    assert!(!points.is_empty());
    assert!(points.len() >= 20);

    // Minimum should be near x=0
    let min_point = points
        .iter()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .unwrap();
    assert!(
        min_point.0.abs() < 0.3,
        "minimum x should be near 0, got {}",
        min_point.0
    );
}

/// Auto-range scan: default is ±2*error.
#[test]
fn scan_auto_range() {
    let result = MnMigrad::new()
        .add("x", 0.0, 1.0)
        .minimize(&|p: &[f64]| p[0] * p[0]);

    assert!(result.is_valid());

    let scan = MnScan::new(&|p: &[f64]| p[0] * p[0], &result);
    // low == high == 0.0 triggers auto-range
    let points = scan.scan(0, 10, 0.0, 0.0);

    assert!(!points.is_empty());
    // Check range covers roughly ±2*error
    let x_min = points.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
    let x_max = points.iter().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max);
    assert!(x_min < -0.5, "auto-range should go below 0, got {x_min}");
    assert!(x_max > 0.5, "auto-range should go above 0, got {x_max}");
}

/// Scan finds better minimum.
#[test]
fn scan_minimum_tracking() {
    use minuit2::scan::MnParameterScan;
    use minuit2::user_parameters::MnUserParameters;

    // Start far from minimum.
    let mut params = MnUserParameters::new();
    params.add("x", 10.0, 1.0);

    let fcn = |p: &[f64]| (p[0] - 3.0).powi(2);
    let initial_fval = fcn(&[10.0]);

    let mut scanner = MnParameterScan::new(&fcn, params, initial_fval);
    let _points = scanner.scan(0, 50, 0.0, 6.0);

    // Scanner should have found a better point near x=3.
    assert!(
        scanner.fval() < initial_fval,
        "scanner should find better fval"
    );
    assert!(
        (scanner
            .params()
            .value("x")
            .expect("parameter x should exist")
            - 3.0)
            .abs()
            < 0.5,
        "tracked minimum should move from 10 to near 3"
    );
}

/// nsteps is clamped to [2, 101], and bounded scan ranges are clamped to limits.
#[test]
fn scan_range_and_step_clamping() {
    let result = MnMigrad::new()
        .add_limited("x", 3.0, 1.0, 0.0, 2.0)
        .add("y", -1.0, 1.0)
        .minimize(&|p: &[f64]| (p[0] - 1.0).powi(2) + p[1] * p[1]);

    assert!(result.is_valid());

    let scan = MnScan::new(&|p: &[f64]| (p[0] - 1.0).powi(2) + p[1] * p[1], &result);

    // nsteps is clamped below range to 2.
    let points = scan.scan(0, 0, -10.0, 10.0);
    assert_eq!(points.len(), 3);

    let x_min = points.iter().map(|v| v.0).fold(f64::INFINITY, f64::min);
    let x_max = points.iter().map(|v| v.0).fold(f64::NEG_INFINITY, f64::max);
    assert!((x_min - 0.0).abs() < 1e-12);
    assert!((x_max - 2.0).abs() < 1e-12);
}

/// `scan()` and `scan_serial()` should be equivalent.
#[test]
fn scan_default_matches_serial() {
    let result = MnMigrad::new()
        .add("x", 5.0, 1.0)
        .add("y", -3.0, 1.0)
        .minimize(&|p: &[f64]| p[0] * p[0] + p[1] * p[1]);

    assert!(result.is_valid());

    let scan = MnScan::new(&|p: &[f64]| p[0] * p[0] + p[1] * p[1], &result);
    let points_default = scan.scan(0, 25, -2.5, 2.5);
    let points_serial = scan.scan_serial(0, 25, -2.5, 2.5);

    assert_eq!(points_default.len(), points_serial.len());
    for (a, b) in points_default.iter().zip(points_serial.iter()) {
        assert!((a.0 - b.0).abs() < 1e-12);
        assert!((a.1 - b.1).abs() < 1e-12);
    }
}

/// Parallel scan should match serial results.
#[cfg(feature = "parallel")]
#[test]
fn scan_parallel_matches_serial() {
    let result = MnMigrad::new()
        .add("x", 1.5, 0.5)
        .add("y", -0.5, 0.5)
        .minimize(&|p: &[f64]| {
            // Slightly non-trivial shape to exercise full scan path.
            (p[0] - 0.2).powi(2) + 2.0 * (p[1] + 0.4).powi(2) + 0.1 * p[0] * p[1]
        });

    assert!(result.is_valid());

    let fcn = |p: &[f64]| (p[0] - 0.2).powi(2) + 2.0 * (p[1] + 0.4).powi(2) + 0.1 * p[0] * p[1];
    let scan = MnScan::new(&fcn, &result);
    let serial = scan.scan_serial(0, 40, -2.0, 2.0);
    let parallel = scan.scan_parallel(0, 40, -2.0, 2.0);

    assert_eq!(serial.len(), parallel.len());
    for (a, b) in serial.iter().zip(parallel.iter()) {
        assert!((a.0 - b.0).abs() < 1e-12);
        assert!((a.1 - b.1).abs() < 1e-12);
    }
}
