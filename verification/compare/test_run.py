"""Exercise matched-accuracy selection and the scientific verdict thresholds."""

import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "compare_run", Path(__file__).with_name("run.py")
)
harness = importlib.util.module_from_spec(spec)
spec.loader.exec_module(harness)


def row(
    group="EDM",
    *,
    rust_calls=100,
    python_calls=100,
    rust_error=1.0,
    python_error=1.0,
    distance=0.0,
    nonquadratic=True,
    **extra,
):
    return dict(
        group=group,
        case="example",
        tolerance=0.1,
        nonquadratic=nonquadratic,
        rust_true_distance=distance,
        iminuit_true_distance=distance,
        rust_value_calls=rust_calls,
        iminuit_value_calls=python_calls,
        rust_hesse_calls=rust_calls,
        iminuit_hesse_calls=python_calls,
        rust_cov_error=rust_error,
        iminuit_cov_error=python_error,
        **extra,
    )


def test_calls_threshold_includes_both_endpoints():
    assert harness.verdict([row(rust_calls=90)], "EDM").startswith("no difference;")
    assert harness.verdict([row(rust_calls=110)], "EDM").startswith("no difference;")
    assert harness.verdict([row(rust_calls=89)], "EDM").startswith("difference;")
    assert harness.verdict([row(rust_calls=111)], "EDM").startswith("difference;")


def test_covariance_threshold_includes_both_endpoints():
    for error in (0.5, 2.0):
        assert harness.verdict([row(rust_error=error)], "EDM").startswith(
            "no difference;"
        )
    for error in (0.499, 2.001):
        assert harness.verdict([row(rust_error=error)], "EDM").startswith("difference;")


def test_unmatched_rows_cannot_change_medians():
    rows = [row(), row(rust_calls=1000, distance=1.000001e-6)]
    assert harness.medians(rows, "EDM") == (1.0, 1.0, 1)
    rows[1]["rust_true_distance"] = 0.0
    assert harness.medians(rows, "EDM") == (1.0, 1.0, 1)
    rows[1]["iminuit_true_distance"] = 1e-6
    assert harness.medians(rows, "EDM") == (5.5, 1.0, 2)


def test_quadratics_cannot_change_hessian_verdict():
    rows = [row("Hessian"), row("Hessian", rust_calls=1000, nonquadratic=False)]
    assert harness.medians(rows, "Hessian") == (1.0, 1.0, 1)


def test_hessian_ratios_use_only_hesse_calls():
    record = row("Hessian", rust_calls=1000)
    record["rust_hesse_calls"] = 100
    assert harness.medians([record], "Hessian") == (1.0, 1.0, 1)


def test_any_bounds_success_difference_counts_even_without_matched_accuracy():
    record = row("bounds", distance=1.0, rust_success=False, iminuit_success=True)
    assert harness.verdict([record], "bounds").startswith("difference;")


def test_first_tolerance_is_selected_for_each_library_independently():
    loose = row(rust_calls=50, python_calls=60)
    loose["tolerance"] = 1.0
    loose["rust_true_distance"] = 1e-4
    tight = row(rust_calls=70, python_calls=80)
    assert harness.matched_counts([loose, tight]) == [
        ("example", {"rust": ("70", "0.1"), "iminuit": ("60", "1.0")})
    ]


def test_missing_and_zero_errors_have_explicit_conventions():
    assert harness.ratio(None, 1) is None
    assert harness.ratio(0, 0) == 1
    assert harness.ratio(1, 0) == float("inf")
    assert "not reached" in harness.verdict([], "gradient")


def test_nonfinite_covariance_remains_a_failed_row():
    from types import SimpleNamespace

    import mpmath as mp

    case = SimpleNamespace(
        free=[0],
        reference=[mp.mpf(0)],
        minimum=mp.mpf(0),
        covariance=mp.matrix([[1]]),
        value=lambda p, precise: p[0] ** 2,
    )
    snapshot = {
        "valid": True,
        "params": [0.0],
        "fval": 0.0,
        "edm": 0.0,
        "covariance": [[None]],
        "value_calls": 10,
        "gradient_calls": 0,
    }
    result = {
        "migrad": snapshot,
        "hesse": snapshot,
        "hesse_value_calls": 2,
        "hesse_gradient_calls": 0,
    }
    measured = harness.metrics(SimpleNamespace(mp=mp), case, "Hessian", result)
    assert measured["status"] == "failed"
    assert measured["message"] == "covariance contains a nonfinite entry"
    assert measured["cov_error"] is None
    assert measured["value_calls"] == 10


def test_fixed_covariance_full_and_reduced_forms_agree():
    from types import SimpleNamespace

    import mpmath as mp

    case = SimpleNamespace(
        free=[1], reference=[0.4, -2.0], covariance=mp.matrix([[0.25]])
    )
    assert harness.covariance_error(case, [[0.0, 0.0], [0.0, 0.5]]) == 1.0
    assert harness.covariance_error(case, [[0.5]]) == 1.0
