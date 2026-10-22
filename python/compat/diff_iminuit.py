"""Differential drop-in check: run identical user code against iminuit.Minuit
and minuit2.Minuit, then report where minuit2 diverges.

A check is a callable `fn(Minuit) -> observable`. For each check we run it with
the real iminuit class (reference) and the minuit2 class, then classify:

  PASS      both succeed and observables match numerically
  MISMATCH  both succeed but numbers differ
  GAP       iminuit succeeds, minuit2 raises (missing/incompatible API)
  BOTH-FAIL both raise (not a drop-in concern)

Run:  python python/compat/diff_iminuit.py
"""

import numpy as np
from iminuit import Minuit as IMinuit
from minuit2 import Minuit as RMinuit


def cost(a, b):
    # least-squares-like; minimum at (1, 2), fval = 0
    return (a - 1.0) ** 2 + (b - 2.0) ** 2


def _fit(M):
    m = M(cost, a=0.0, b=0.0)
    m.migrad()
    return m


# --- checks -----------------------------------------------------------------
def c_construct_kwargs(M):
    m = M(cost, a=0.0, b=0.0)
    m.migrad()
    return [m.values["a"], m.values["b"], m.fval]


def c_construct_positional(M):
    m = M(cost, 0.0, 0.0)  # positional starting values, names from signature
    m.migrad()
    return [m.values["a"], m.values["b"]]


def c_values_by_index(M):
    m = _fit(M)
    return [m.values[0], m.values[1]]


def c_values_setitem_persists(M):
    m = M(cost, a=0.0, b=0.0)
    m.values["a"] = 3.0  # item assignment must persist on the object
    return [m.values["a"]]


def c_limits_setitem_persists(M):
    m = M(cost, a=0.0, b=0.0)
    m.limits["a"] = (-5.0, 5.0)
    lo, hi = m.limits["a"]
    return [float(lo), float(hi)]


def c_errordef_constant(M):
    # class constant + assignment
    _ = M.LEAST_SQUARES
    m = M(cost, a=0.0, b=0.0)
    m.errordef = M.LEAST_SQUARES
    return [m.errordef]


def c_errors_and_valid(M):
    m = _fit(M)
    return [m.errors["a"], m.errors["b"], float(m.valid)]


def c_fixed_view_assign(M):
    m = M(cost, a=0.0, b=5.0)
    m.fixed["b"] = True
    m.migrad()
    return [m.values["a"], m.values["b"]]  # b should stay 5.0


def c_limits_view_assign(M):
    m = M(cost, a=0.5, b=0.0)  # off-boundary start (realistic usage)
    m.limits["a"] = (0.0, None)  # one-sided lower limit
    m.migrad()
    return [m.values["a"]]


def c_param_started_at_limit(M):
    # Starting a parameter exactly on its bound: the negative-curvature seed
    # escape (see src/migrad/seed.rs) moves it off the singular point so it
    # reaches the interior minimum, matching ROOT/iminuit.
    m = M(cost, a=0.0, b=0.0)
    m.limits["a"] = (0.0, None)
    m.migrad()
    return [m.values["a"]]


def c_fmin_fields(M):
    m = _fit(M)
    fm = m.fmin
    return [float(fm.is_valid), fm.edm]


def c_params_objects(M):
    m = _fit(M)
    p = m.params[0]
    return [p.name == "a", p.value]


def c_counts(M):
    m = _fit(M)
    return [m.npar, m.nfit, list(m.parameters) == ["a", "b"]]


def c_hesse_covariance(M):
    m = _fit(M).hesse()
    cov = m.covariance
    return [cov[0][0], cov[1][1]]


def c_minos_merrors(M):
    m = _fit(M)
    m.minos()
    e = m.merrors["a"]
    return [e.lower, e.upper, float(e.is_valid)]


def c_profile(M):
    m = _fit(M)
    x, y = m.profile("a", size=20, bound=2.0, subtract_min=True)
    return [len(x), len(y), float(np.min(y))]


def c_mncontour(M):
    m = _fit(M)
    pts = m.mncontour("a", "b", size=20)
    return [pts.shape[1], pts.shape[0] >= 20]


def c_contour_grid(M):
    m = _fit(M)
    xg, yg, fg = m.contour("a", "b", size=10)
    return [len(xg), len(yg), fg.shape == (10, 10)]


def c_reset(M):
    m = _fit(M)
    m.reset()
    return [m.values["a"]]  # back to initial 0.0


def c_fixto(M):
    m = M(cost, a=0.0, b=0.0)
    m.fixto("b", 3.0)
    m.migrad()
    return [m.values["b"]]  # fixed to 3.0


def c_values_iteration(M):
    m = _fit(M)
    return list(m.values)  # iterating a ValueView yields the values in order


def c_values_to_dict(M):
    m = _fit(M)
    d = m.values.to_dict()
    return [d["a"], d["b"]]


def c_errors_by_index(M):
    m = _fit(M)
    return [m.errors[0], m.errors[1]]


def c_negative_index(M):
    m = _fit(M)
    return [m.values[-1]]  # -1 -> last parameter (b)


def c_default_errordef(M):
    m = M(cost, a=0.0, b=0.0)
    return [m.errordef]  # iminuit defaults to 1.0 (least-squares)


def c_multi_param_minos(M):
    m = _fit(M)
    m.minos()
    ea, eb = m.merrors["a"], m.merrors["b"]
    return [ea.lower, ea.upper, eb.lower, eb.upper]


def c_fix_release_roundtrip(M):
    m = M(cost, a=0.0, b=5.0)
    m.fixed["b"] = True
    m.migrad()  # b held at 5.0
    m.fixed["b"] = False
    m.migrad()  # b now free -> 2.0
    return [m.values["a"], m.values["b"]]


def c_scan_quadratic_bad_start(M):
    m = M(cost, a=-4.0, b=7.0)
    m.errors["a"] = 6.0
    m.errors["b"] = 6.0
    m.scan(ncall=100)
    return [m.values["a"], m.values["b"], m.fval]


def c_scan_bounded_parameter(M):
    m = M(cost, a=-4.0, b=7.0)
    m.limits["a"] = (0.0, 2.0)
    m.errors["b"] = 6.0
    m.scan(ncall=100)
    return [m.values["a"], m.values["b"], m.fval]


CHECKS = [
    ("construct_kwargs", c_construct_kwargs),
    ("construct_positional", c_construct_positional),
    ("values_by_index", c_values_by_index),
    ("values_setitem_persists", c_values_setitem_persists),
    ("limits_setitem_persists", c_limits_setitem_persists),
    ("errordef_constant", c_errordef_constant),
    ("errors_and_valid", c_errors_and_valid),
    ("fixed_view_assign", c_fixed_view_assign),
    ("limits_view_assign", c_limits_view_assign),
    ("start_at_limit", c_param_started_at_limit),
    ("fmin_fields", c_fmin_fields),
    ("params_objects", c_params_objects),
    ("counts", c_counts),
    ("hesse_covariance", c_hesse_covariance),
    ("minos_merrors", c_minos_merrors),
    ("profile", c_profile),
    ("mncontour", c_mncontour),
    ("contour_grid", c_contour_grid),
    ("reset", c_reset),
    ("fixto", c_fixto),
    ("values_iteration", c_values_iteration),
    ("values_to_dict", c_values_to_dict),
    ("errors_by_index", c_errors_by_index),
    ("negative_index", c_negative_index),
    ("default_errordef", c_default_errordef),
    ("multi_param_minos", c_multi_param_minos),
    ("fix_release_roundtrip", c_fix_release_roundtrip),
    ("scan_quadratic_bad_start", c_scan_quadratic_bad_start),
    ("scan_bounded_parameter", c_scan_bounded_parameter),
]


def _run(fn, M):
    try:
        return ("ok", fn(M))
    except Exception as e:  # noqa: BLE001
        return ("err", f"{type(e).__name__}: {e}")


def _to_floats(v):
    out = []
    for x in v:
        out.append(float(x) if not isinstance(x, bool) else float(x))
    return out


def main():
    rows = []
    for name, fn in CHECKS:
        si, ri = _run(fn, IMinuit)
        sr, rr = _run(fn, RMinuit)
        if si == "ok" and sr == "ok":
            try:
                a = _to_floats(ri)
                b = _to_floats(rr)
                if len(a) == len(b) and np.allclose(a, b, rtol=1e-2, atol=1e-2, equal_nan=True):
                    status = "PASS"
                    detail = ""
                else:
                    status = "MISMATCH"
                    detail = f"iminuit={a} minuit2={b}"
            except Exception as e:  # noqa: BLE001
                status = "MISMATCH"
                detail = f"compare-error {e}: iminuit={ri} minuit2={rr}"
        elif si == "ok" and sr == "err":
            status = "GAP"
            detail = f"minuit2 -> {rr}"
        elif si == "err" and sr == "err":
            status = "BOTH-FAIL"
            detail = f"iminuit -> {ri}"
        else:
            status = "MINUIT2-ONLY"
            detail = f"iminuit -> {ri}; minuit2 ok"
        rows.append((name, status, detail))

    width = max(len(n) for n, _, _ in rows)
    counts = {}
    print(f"\n{'CHECK':<{width}}  STATUS     DETAIL")
    print("-" * (width + 60))
    for name, status, detail in rows:
        counts[status] = counts.get(status, 0) + 1
        print(f"{name:<{width}}  {status:<9}  {detail[:90]}")
    print("-" * (width + 60))
    summary = "  ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    print(f"SUMMARY: {summary}")


if __name__ == "__main__":
    main()
