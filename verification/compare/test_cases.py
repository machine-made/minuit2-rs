"""Check objective values, analytic gradients, and reference curvature."""

import importlib.util
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

# Load only the adjacent module when Python runs in isolated mode.
_spec = importlib.util.spec_from_file_location(
    "compare_cases", Path(__file__).with_name("cases.py")
)
cases_module = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = cases_module
_spec.loader.exec_module(cases_module)
ALL = cases_module.cases() + cases_module.bounds_cases()


@pytest.mark.parametrize("case", ALL, ids=lambda case: case.name)
def test_reference_value(case):
    actual = case.value([float(x) for x in case.reference])
    assert abs(mp.mpf(actual) - case.minimum) <= mp.mpf("1e-12") * max(
        1, abs(case.minimum)
    )


@pytest.mark.parametrize("case", ALL, ids=lambda case: case.name)
def test_gradient(case):
    # Starts exercise nonzero derivatives; references exercise cancellation near the minimum.
    for p in (case.spec["start"], [float(x) for x in case.reference]):
        actual = case.gradient(p)
        for i, g in enumerate(actual):
            orders = tuple(int(k == i) for k in range(len(p)))
            reference = mp.diff(
                lambda *x: case.value(x, True), tuple(map(mp.mpf, p)), orders
            )
            assert abs(mp.mpf(g) - reference) <= mp.mpf("1e-10") * max(
                1, abs(reference)
            )


@pytest.mark.parametrize("case", ALL, ids=lambda case: case.name)
def test_covariance_convention(case):
    h = np.array([[float(case.hessian[i, j]) for j in case.free] for i in case.free])
    covariance = np.array(case.covariance.tolist(), dtype=float)
    assert np.all(np.linalg.eigvalsh(h) > 0)
    np.testing.assert_allclose(
        h @ covariance, 2 * case.spec["errordef"] * np.eye(len(case.free)), atol=1e-10
    )


@pytest.mark.parametrize("case", ALL, ids=lambda case: case.name)
def test_gradient_relative_away_from_minimum(case):
    for fraction in (0.0, 0.4):
        p = [
            (1 - fraction) * start + fraction * float(reference)
            for start, reference in zip(case.spec["start"], case.reference)
        ]
        actual = case.gradient(p)
        for i, g in enumerate(actual):
            orders = tuple(int(k == i) for k in range(len(p)))
            reference = mp.diff(
                lambda *x: case.value(x, True), tuple(map(mp.mpf, p)), orders
            )
            assert abs(mp.mpf(g) - reference) <= mp.mpf("1e-10") * abs(
                reference
            ) + mp.mpf("1e-30")
