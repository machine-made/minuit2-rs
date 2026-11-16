#!/usr/bin/env python3
"""
Run eight start-on-limit cases of quadratic2 through iminuit.

Each case minimises quadratic2(x, y) = (x-1)^2 + 4(y+2)^2 + 0.3xy. The closed-form
minimum is x = 1.30735, y = -2.04903, f = -0.69956003. The Rust results are
checked in tests/limit_start_two_params.rs: fval within 1e-4 of that value,
y below -1.5 and x within 1e-2 of 1.3074.

Run:
    uv run --with iminuit==2.33.0 python scripts/limit_start_iminuit.py
In the repository root, add --no-project after `uv run`: pyproject.toml would
otherwise make uv build the Rust extension first.

Output: one line per case with fval, values, nfcn and valid.
"""

from iminuit import Minuit

NONE = None

# name, x start, x limits, y start, y limits, step of x, step of y
CASES = [
    ("1 both limited, y on upper end", 0.4, (0.0, 2.0), -1.0, (-3.0, -1.0), 0.1, 0.1),
    ("2 case 1 with step 0.5 for y", 0.4, (0.0, 2.0), -1.0, (-3.0, -1.0), 0.1, 0.5),
    ("3 x free, y limited on upper end", 0.4, (NONE, NONE), -1.0, (-3.0, -1.0), 0.1, 0.1),
    ("4 x lower limit, y limited on upper end", 0.4, (0.0, NONE), -1.0, (-3.0, -1.0), 0.1, 0.1),
    ("5 both limited, y on lower end", 0.4, (0.0, 2.0), -3.0, (-3.0, -1.0), 0.1, 0.1),
    ("6 x free, y upper limit only", 0.4, (NONE, NONE), -1.0, (NONE, -1.0), 0.1, 0.1),
    ("7 both on one-sided limits", 0.0, (0.0, NONE), -1.0, (NONE, -1.0), 0.1, 0.1),
    ("8 both on double limits", 0.0, (0.0, 2.0), -1.0, (-3.0, -1.0), 0.1, 0.1),
]


def quadratic2(x, y):
    return (x - 1.0) ** 2 + 4.0 * (y + 2.0) ** 2 + 0.3 * x * y


def run_case(x0, xlim, y0, ylim, xstep, ystep):
    m = Minuit(quadratic2, x=x0, y=y0)
    m.errordef = 1.0
    m.limits = [xlim, ylim]
    m.errors = (xstep, ystep)
    m.tol = 0.1
    m.migrad()
    return m


def main():
    import iminuit

    print(f"iminuit {iminuit.__version__}")
    print("closed form: fval -0.69956003, x 1.30735, y -2.04903")
    for name, x0, xlim, y0, ylim, xstep, ystep in CASES:
        m = run_case(x0, xlim, y0, ylim, xstep, ystep)
        values = ", ".join(f"{v:.9g}" for v in m.values)
        print(
            f"case {name}: fval {m.fval:.9g}, values ({values}), "
            f"nfcn {m.nfcn}, valid {m.valid}"
        )


if __name__ == "__main__":
    main()
