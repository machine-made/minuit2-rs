"""Runtime tests for phase-2 iminuit-compatible scan/profile/contour methods.

Run after `maturin develop`:
    pytest python/tests/test_phase2.py -q
"""

import numpy as np
import pytest
from minuit2 import Minuit


def quad(x, y):
    # Minimum at (1, 2), fval = 0; x and y are independent.
    return (x - 1.0) ** 2 + (y - 2.0) ** 2


def fitted():
    m = Minuit(quad, x=0.0, y=0.0)
    m.migrad()
    return m


def test_profile_returns_two_arrays_others_fixed():
    m = fitted()
    xs, fs = m.profile("x", size=11, bound=2.0, subtract_min=True)
    assert isinstance(xs, np.ndarray) and isinstance(fs, np.ndarray)
    assert xs.shape == (11,) and fs.shape == (11,)
    # profile of (x-1)^2 with y fixed: parabola with interior minimum near x=1.
    # (The grid need not land exactly on x=1, so the min is small but not zero.)
    assert fs.min() < fs[0] and fs.min() < fs[-1]
    assert abs(xs[fs.argmin()] - 1.0) < 0.5


def test_mnprofile_reminimizes_and_reports_ok():
    m = fitted()
    xs, fs, ok = m.mnprofile("x", size=7, bound=2.0, subtract_min=True)
    assert xs.shape == (7,) and fs.shape == (7,) and ok.shape == (7,)
    assert ok.dtype == np.bool_
    assert ok.all()  # independent params => every re-minimization converges
    assert fs.min() < 1e-6


def test_mncontour_is_closed_nx2():
    m = fitted()
    pts = m.mncontour("x", "y", size=20)
    assert pts.ndim == 2 and pts.shape[1] == 2
    assert pts.shape[0] == 21  # closed: last point appended == first
    assert np.allclose(pts[0], pts[-1])


def test_mncontour_cl_not_implemented():
    m = fitted()
    with pytest.raises(NotImplementedError):
        m.mncontour("x", "y", cl=0.9)


def test_contour_grid_shapes():
    m = fitted()
    xg, yg, fval2d = m.contour("x", "y", size=8, bound=2.0, subtract_min=True)
    assert xg.shape == (8,) and yg.shape == (8,)
    assert fval2d.shape == (8, 8)
    # grid holds others fixed; the interior (near the minimum) is far lower than
    # the corners. The grid need not land exactly on (1,2), so min is not zero.
    assert fval2d.min() < fval2d[0, 0]
    assert fval2d.min() < fval2d[-1, -1]


def test_scan_finds_best_hypercube_point_from_bad_start():
    m = Minuit(quad, x=-4.0, y=7.0)
    m.errors["x"] = 6.0
    m.errors["y"] = 6.0
    m.scan(ncall=100)
    try:
        assert abs(m.values["x"] - 0.6666666667) < 1e-9
        assert abs(m.values["y"] - 2.3333333333) < 1e-9
        assert m.fval < 0.23
    except AssertionError:
        print("minuit2 scan values", m.values.to_dict(), "fval", m.fval)
        raise


def test_scan_uses_limits_and_includes_boundaries():
    m = Minuit(quad, x=-4.0, y=7.0)
    m.limits["x"] = (0.0, 2.0)
    m.errors["y"] = 6.0
    m.scan(ncall=100)
    try:
        assert abs(m.values["x"] - 0.8888888889) < 1e-9
        assert abs(m.values["y"] - 2.3333333333) < 1e-9
        assert m.fval < 0.13
        assert m.covariance is None
    except AssertionError:
        print("minuit2 scan values", m.values.to_dict(), "fval", m.fval)
        raise
