"""Runtime smoke test for the phase-1 iminuit-compatible API.

Run after `maturin develop --features python`:
    pytest python/tests/test_smoke.py -q
"""

import minuit2
from minuit2 import Minuit


def quad(x, y):
    # Minimum at (1, 2), fval = 0.
    return (x - 1.0) ** 2 + (y - 2.0) ** 2


def test_version_matches_cargo():
    assert minuit2.__version__  # non-empty, sourced from CARGO_PKG_VERSION
    assert minuit2.__version__[0].isdigit()


def test_migrad_chaining_and_results():
    m = Minuit(quad, x=0.0, y=0.0)
    out = m.migrad()
    assert out is m  # chaining: migrad returns self
    assert m.valid
    assert abs(m.values["x"] - 1.0) < 1e-3
    assert abs(m.values["y"] - 2.0) < 1e-3
    assert m.fval < 1e-6
    assert m.npar == 2
    assert m.nfit == 2
    assert tuple(m.parameters) == ("x", "y")
    assert m.nfcn > 0


def test_fmin_and_params_objects():
    m = Minuit(quad, x=0.0, y=0.0).migrad()
    fm = m.fmin
    assert fm is not None
    assert fm.is_valid
    assert fm.errordef == 1.0
    assert fm.nfcn > 0
    params = m.params
    assert len(params) == 2
    assert params[0].name == "x"
    assert params[0].is_fixed is False
    assert abs(params[0].value - 1.0) < 1e-3


def test_errordef_strategy_tol_roundtrip():
    m = Minuit(quad, x=0.0, y=0.0)
    m.errordef = 0.5
    assert m.errordef == 0.5
    m.strategy = 2
    assert m.strategy == 2
    m.tol = 1e-4
    assert m.tol == 1e-4
    try:
        m.strategy = 5
        raised = False
    except ValueError:
        raised = True
    assert raised  # strategy must be 0/1/2


def test_hesse_then_minos_populates_merrors():
    m = Minuit(quad, x=0.0, y=0.0).migrad().hesse()
    assert m.covariance is not None
    out = m.minos()
    assert out is m  # minos returns self
    assert set(m.merrors.keys()) == {"x", "y"}
    ex = m.merrors["x"]
    assert ex.is_valid
    assert ex.lower < 0.0 < ex.upper  # iminuit sign convention
    assert ex.upper == abs(ex.lower) or ex.upper > 0.0
