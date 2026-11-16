#!/usr/bin/env python3
"""
Run seven start-on-limit cases of quadratic2 through iminuit and check them.

Each case minimises quadratic2(x, y) = (x-1)^2 + 4(y+2)^2 + 0.3xy. The closed-form
minimum is x = 1.30735, y = -2.04903, f = -0.69956003. The script exits with
status 1 unless every case is valid and its fval is within 1e-4 of that value.
The same value and tolerance are the check of the Rust results in
tests/limit_start_two_params.rs.

Run:
    uv run --with iminuit==2.33.0 python scripts/limit_start_iminuit.py
In the repository root, add --no-project after `uv run`: pyproject.toml would
otherwise make uv build the Rust extension first.

Output: one line per case with fval, values, nfcn and valid, then one line with
the number of cases that failed the check.
"""

from collections.abc import Sequence
from typing import Final, TypeAlias, cast

# iminuit 2.33.0 has no type stubs, so basedpyright strict reports them.
import iminuit  # pyright: ignore[reportMissingTypeStubs]
from iminuit import Minuit  # pyright: ignore[reportMissingTypeStubs]

NONE: Final[None] = None

# lower and upper limit of one parameter; None leaves that side open
Limits: TypeAlias = tuple[float | None, float | None]
# name, x start, x limits, y start, y limits, step of x, step of y
Case: TypeAlias = tuple[str, float, Limits, float, Limits, float, float]

CASES: Final[list[Case]] = [
    ("1 both limited, y on upper end", 0.4, (0.0, 2.0), -1.0, (-3.0, -1.0), 0.1, 0.1),
    ("2 case 1 with step 0.5 for y", 0.4, (0.0, 2.0), -1.0, (-3.0, -1.0), 0.1, 0.5),
    ("3 x free, y limited on upper end", 0.4, (NONE, NONE), -1.0, (-3.0, -1.0), 0.1, 0.1),
    ("4 x lower limit, y limited on upper end", 0.4, (0.0, NONE), -1.0, (-3.0, -1.0), 0.1, 0.1),
    ("5 x free, y upper limit only", 0.4, (NONE, NONE), -1.0, (NONE, -1.0), 0.1, 0.1),
    ("6 both on one-sided limits", 0.0, (0.0, NONE), -1.0, (NONE, -1.0), 0.1, 0.1),
    ("7 both on double limits", 0.0, (0.0, 2.0), -1.0, (-3.0, -1.0), 0.1, 0.1),
]

CLOSED_FORM_FVAL: Final[float] = -0.69956003
# 1e-4: the bound on fval of tests/limit_start_two_params.rs. iminuit 2.33.0 ends
# within 2.9e-5 of the closed form in these cases, so 1e-4 is 3 times that.
FVAL_TOLERANCE: Final[float] = 1e-4


def quadratic2(x: float, y: float) -> float:
    return (x - 1.0) ** 2 + 4.0 * (y + 2.0) ** 2 + 0.3 * x * y


def run_case(
    x0: float,
    xlim: Limits,
    y0: float,
    ylim: Limits,
    xstep: float,
    ystep: float,
) -> Minuit:
    # iminuit reads the names x and y from the signature of quadratic2; its type
    # hint for the function is (*args: float) -> float.
    m = Minuit(quadratic2, x=x0, y=y0)  # pyright: ignore[reportArgumentType]
    m.errordef = 1.0
    m.limits = [xlim, ylim]
    m.errors = (xstep, ystep)
    m.tol = 0.1
    m.migrad()
    return m


def main() -> None:
    print(f"iminuit {iminuit.__version__}")
    print("closed form: fval -0.69956003, x 1.30735, y -2.04903")
    failed = 0
    for name, x0, xlim, y0, ylim, xstep, ystep in CASES:
        m = run_case(x0, xlim, y0, ylim, xstep, ystep)
        fval = m.fval
        if fval is None:
            raise SystemExit(f"case {name}: iminuit returned no fval")
        values = ", ".join(f"{v:.9g}" for v in cast("Sequence[float]", m.values))
        print(
            f"case {name}: fval {fval:.9g}, values ({values}), "
            f"nfcn {m.nfcn}, valid {m.valid}"
        )
        if not m.valid or abs(fval - CLOSED_FORM_FVAL) >= FVAL_TOLERANCE:
            failed += 1
    print(f"cases failing the check: {failed}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
