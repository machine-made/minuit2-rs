"""Run a batch with counters inside the Python objective and gradient."""

import importlib.util
import json
import math
import sys
from pathlib import Path

import iminuit
from iminuit import Minuit
from iminuit._core import MnStrategy

module_spec = importlib.util.spec_from_file_location(
    "compare_cases", Path(__file__).with_name("cases.py")
)
cases = importlib.util.module_from_spec(module_spec)
sys.modules[module_spec.name] = cases
module_spec.loader.exec_module(cases)


def run(spec):
    values = 0
    gradients = 0

    def objective(*p):
        nonlocal values
        values += 1
        return cases.evaluate(spec, p)

    def gradient(*p):
        nonlocal gradients
        gradients += 1
        return cases.gradient(spec, p)

    def snapshot(m):
        return {
            "valid": bool(m.valid),
            "fval": float(m.fval),
            "params": list(m.values),
            "edm": float(m.fmin.edm),
            "value_calls": values,
            "gradient_calls": gradients,
            "covariance": None
            if m.covariance is None
            else [
                [float(m.covariance[i, j]) for j in range(len(spec["start"]))]
                for i in range(len(spec["start"]))
            ],
        }

    try:
        m = Minuit(
            objective,
            *spec["start"],
            grad=gradient if spec.get("gradient", False) else None,
        )
        m.errordef = spec["errordef"]
        m.errors = spec["steps"]
        m.limits = spec["limits"]
        for i in spec["fixed"]:
            m.fixed[i] = True
        # The pinned public setter rejects level 3; the bundled constructor accepts it.
        if spec["strategy"] == 3:
            m._strategy = MnStrategy(3)
            assert m.strategy.strategy == 3
        else:
            m.strategy = spec["strategy"]
        m.tol = spec["tolerance"]
        # One Migrad invocation matches the Rust API and avoids wrapper retries.
        m.migrad(
            ncall=200 + 100 * m.nfit + 5 * m.nfit * m.nfit, iterate=1, use_simplex=False
        )
        before = snapshot(m)
        old_values, old_gradients = values, gradients
        try:
            m.hesse(ncall=200 + 100 * m.nfit + 5 * m.nfit * m.nfit)
            return {
                "migrad": before,
                "hesse": snapshot(m),
                "hesse_value_calls": values - old_values,
                "hesse_gradient_calls": gradients - old_gradients,
            }
        except Exception as error:  # noqa: BLE001
            # Every library exception must remain a measured failed row.
            return {
                "migrad": before,
                "error": f"Hesse failed: {type(error).__name__}: {error}",
                "hesse_value_calls": values - old_values,
                "hesse_gradient_calls": gradients - old_gradients,
            }
    except Exception as error:  # noqa: BLE001
        # Every library exception must remain a measured failed row.
        return {
            "error": f"{type(error).__name__}: {error}",
            "value_calls": values,
            "gradient_calls": gradients,
        }


def finite_json(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: finite_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [finite_json(item) for item in value]
    return value


if __name__ == "__main__":
    specs = json.load(sys.stdin)
    print(
        json.dumps(
            finite_json(
                {
                    "version": iminuit.__version__,
                    "results": [run(spec) for spec in specs],
                }
            ),
            allow_nan=False,
        )
    )
