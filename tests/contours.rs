use minuit2::MnMinos;
use minuit2::{MnContours, MnHesse, MnMigrad};

/// 2D quadratic: contour points should form approximate ellipse.
#[test]
fn contour_quadratic_ellipse() {
    let quadratic = |p: &[f64]| 2.0 * p[0] * p[0] + 8.0 * p[1] * p[1];

    let result = MnMigrad::new()
        .add("x", 5.0, 1.0)
        .add("y", -3.0, 1.0)
        .minimize(&quadratic);

    assert!(result.is_valid());

    let hesse_result = MnHesse::new().calculate(&quadratic, &result);

    let contours = MnContours::new(&quadratic, &hesse_result);
    let contour = contours.contour(0, 1, 8);
    let points = &contour.points;

    // Should have at least 4 points (the cardinal MINOS points)
    assert!(
        points.len() >= 4,
        "contour should have >= 4 points, got {}",
        points.len()
    );

    // All points should be at approximately the same function value: fmin + up
    let up = hesse_result.up();
    let fmin = hesse_result.fval();
    let target = fmin + up;

    assert_eq!(contour.xpar(), 0);
    assert_eq!(contour.ypar(), 1);
    assert_eq!(contour.nfcn(), 0);
    assert!(contour.x_min().is_finite(), "x minimum should be finite");
    assert!(contour.y_min().is_finite(), "y minimum should be finite");

    for (x, y) in points {
        let f = quadratic(&[*x, *y]);
        // Contour points should be approximately at F = fmin + up
        // Allow generous tolerance since contour fitting is approximate
        assert!(
            (f - target).abs() < 0.5 * up,
            "contour point ({x}, {y}) has f={f}, expected ~{target}"
        );
    }
}

/// Correlated 2D quadratic: contour should be tilted ellipse.
#[test]
fn contour_correlated() {
    // f(x,y) = x^2 + y^2 + x*y
    let correlated = |p: &[f64]| p[0] * p[0] + p[1] * p[1] + p[0] * p[1];

    let result = MnMigrad::new()
        .add("x", 5.0, 1.0)
        .add("y", -3.0, 1.0)
        .minimize(&correlated);

    assert!(result.is_valid());

    let hesse_result = MnHesse::new().calculate(&correlated, &result);

    let contours = MnContours::new(&correlated, &hesse_result);
    let points = contours.points(0, 1, 8);

    assert!(
        points.len() >= 4,
        "contour should have >= 4 points, got {}",
        points.len()
    );

    // Points should not all lie on axis-aligned ellipse (correlation tilts it)
    // Check that at least one point has both x and y non-zero
    let off_axis = points.iter().any(|(x, y)| x.abs() > 0.01 && y.abs() > 0.01);
    assert!(off_axis, "correlated contour should have off-axis points");
}

#[test]
fn contours_points_respect_minimum_cardinal_count() {
    let result = MnMigrad::new()
        .add("x", 5.0, 1.0)
        .add("y", -3.0, 1.0)
        .minimize(&|p: &[f64]| 2.0 * p[0] * p[0] + 8.0 * p[1] * p[1]);

    assert!(result.is_valid());

    let hesse_result =
        MnHesse::new().calculate(&|p: &[f64]| 2.0 * p[0] * p[0] + 8.0 * p[1] * p[1], &result);
    let contours = MnContours::new(
        &|p: &[f64]| 2.0 * p[0] * p[0] + 8.0 * p[1] * p[1],
        &hesse_result,
    );

    let points = contours.points(0, 1, 2);
    assert_eq!(points.len(), 4);

    let x_min = hesse_result.user_state().parameter(0).value();
    let y_min = hesse_result.user_state().parameter(1).value();

    let minos = MnMinos::new(
        &|p: &[f64]| 2.0 * p[0] * p[0] + 8.0 * p[1] * p[1],
        &hesse_result,
    );
    let x_minos = minos.minos_error(0);
    let y_minos = minos.minos_error(1);

    assert!((points[0].0 - (x_min + x_minos.upper_error())).abs() < 1e-8);
    assert!((points[2].0 - (x_min + x_minos.lower_error())).abs() < 1e-8);
    assert!((points[1].1 - (y_min + y_minos.upper_error())).abs() < 1e-8);
    assert!((points[3].1 - (y_min + y_minos.lower_error())).abs() < 1e-8);
}

#[test]
fn contour_contains_same_points_as_points_call() {
    let result = MnMigrad::new()
        .add("x", 5.0, 1.0)
        .add("y", -3.0, 1.0)
        .minimize(&|p: &[f64]| 2.0 * p[0] * p[0] + 8.0 * p[1] * p[1]);

    assert!(result.is_valid());

    let hesse_result =
        MnHesse::new().calculate(&|p: &[f64]| 2.0 * p[0] * p[0] + 8.0 * p[1] * p[1], &result);
    let contours = MnContours::new(
        &|p: &[f64]| 2.0 * p[0] * p[0] + 8.0 * p[1] * p[1],
        &hesse_result,
    );

    let points = contours.points(0, 1, 12);
    let contour = contours.contour(0, 1, 12);
    assert_eq!(points, contour.points);
}
