"""Deterministic objectives and independent references at 40 decimal digits."""

import math
from dataclasses import dataclass
from functools import lru_cache

import mpmath as mp
import numpy as np

SEED = 20261113
mp.mp.dps = 40


@dataclass
class Case:
    spec: dict
    reference: list
    minimum: object
    hessian: object
    covariance: object
    nonquadratic: bool

    @property
    def name(self):
        return self.spec["name"]

    @property
    def free(self):
        return [i for i in range(len(self.reference)) if i not in self.spec["fixed"]]

    def value(self, p, precise=False):
        return evaluate(self.spec, p, precise)

    def gradient(self, p, precise=False):
        return gradient(self.spec, p, precise)


def evaluate(spec, p, precise=False):
    conv = mp.mpf if precise else float
    exp = mp.exp if precise else math.exp
    log1p = mp.log1p if precise else math.log1p
    p = list(map(conv, p))
    kind = spec["function"]
    if kind == "dense":
        d = [v - conv(c) for v, c in zip(p, spec["center"])]
        return sum(
            conv(0.5) * d[i] * conv(h) * d[j]
            for i, row in enumerate(spec["hessian"])
            for j, h in enumerate(row)
        )
    if kind == "quadratic2":
        x, y = p
        return (x - 1) ** 2 + 4 * (y + 2) ** 2 + conv(0.3) * x * y
    if kind == "rosenbrock":
        x, y = p
        return 100 * (y - x * x) ** 2 + (1 - x) ** 2
    if kind == "bound":
        return (p[0] - conv(spec["target"])) ** 2
    value = conv(0)
    for row in spec["data"]:
        x, y = map(conv, row[:2])
        z = (x - p[1]) / p[2]
        peak = p[0] * exp(-z * z / 2)
        if kind == "gaussian":
            value += ((peak - y) / conv(row[2])) ** 2
        else:
            mu = peak + p[3] * exp(-p[4] * (x + 5))
            if mu <= 0:
                return mp.inf if precise else math.inf
            d = (mu - y) / y
            value += y * (d - log1p(d))
    return value


def gradient(spec, p, precise=False):
    conv = mp.mpf if precise else float
    exp = mp.exp if precise else math.exp
    p = list(map(conv, p))
    kind = spec["function"]
    if kind == "dense":
        d = [v - conv(c) for v, c in zip(p, spec["center"])]
        return [sum(conv(h) * v for h, v in zip(row, d)) for row in spec["hessian"]]
    if kind == "quadratic2":
        x, y = p
        return [2 * (x - 1) + conv(0.3) * y, 8 * (y + 2) + conv(0.3) * x]
    if kind == "rosenbrock":
        x, y = p
        return [-400 * x * (y - x * x) + 2 * (x - 1), 200 * (y - x * x)]
    if kind == "bound":
        return [2 * (p[0] - conv(spec["target"]))]
    out = [conv(0) for _ in p]
    for row in spec["data"]:
        x, y = map(conv, row[:2])
        z = (x - p[1]) / p[2]
        shape = exp(-z * z / 2)
        peak = p[0] * shape
        dm = [shape, peak * z / p[2], peak * z * z / p[2]]
        if kind == "gaussian":
            weight = 2 * (peak - y) / conv(row[2]) ** 2
        else:
            background = exp(-p[4] * (x + 5))
            mu = peak + p[3] * background
            dm += [background, -p[3] * (x + 5) * background]
            weight = 1 - y / mu
        out = [g + weight * d for g, d in zip(out, dm)]
    return out


def make_case(
    name,
    kind,
    start,
    reference,
    *,
    steps=None,
    limits=None,
    fixed=(),
    errordef=1.0,
    **data,
):
    spec = dict(
        name=name,
        function=kind,
        start=start,
        steps=steps or [0.1] * len(start),
        limits=limits or [[None, None] for _ in start],
        fixed=list(fixed),
        errordef=errordef,
        **data,
    )
    reference = list(map(mp.mpf, reference))
    free = [i for i in range(len(start)) if i not in fixed]
    if kind in ("gaussian", "poisson"):
        # Solve stationarity independently of either minimizer.
        reference = list(
            mp.findroot(
                lambda *p: gradient(spec, p, True),
                reference,
                tol=mp.mpf("1e-36"),
                maxsteps=100,
            )
        )
    elif kind == "quadratic2":
        c = mp.mpf(0.3)
        if fixed:
            reference[1] = (-16 - c * reference[0]) / 8
        else:
            reference = list(
                mp.lu_solve(mp.matrix([[2, c], [c, 8]]), mp.matrix([2, -16]))
            )
    minimum = evaluate(spec, reference, True)
    n = len(reference)
    hessian = mp.matrix(n)
    for i in range(n):
        for j in range(i, n):
            orders = tuple(int(k == i) + int(k == j) for k in range(n))
            hessian[i, j] = hessian[j, i] = mp.diff(
                lambda *p: evaluate(spec, p, True), tuple(reference), orders
            )
    reduced = mp.matrix([[hessian[i, j] for j in free] for i in free])
    covariance = 2 * mp.mpf(errordef) * reduced**-1
    return Case(
        spec,
        reference,
        minimum,
        hessian,
        covariance,
        kind in ("rosenbrock", "gaussian", "poisson"),
    )


@lru_cache(maxsize=1)
def cases():
    rng = np.random.default_rng(SEED)
    result = []
    for n in (2, 5, 10, 20):
        a = rng.normal(size=(n, n)) / math.sqrt(n)
        hessian = (a.T @ a + np.eye(n)).tolist()
        center = rng.uniform(-0.5, 0.5, n).tolist()
        result.append(
            make_case(
                f"dense{n}", "dense", [1.5] * n, center, center=center, hessian=hessian
            )
        )
    result += [
        make_case("rosenbrock", "rosenbrock", [0.0, 0.0], [1.0, 1.0]),
        make_case("quadratic2", "quadratic2", [0.4, -1.0], [1.0, -2.0]),
        make_case(
            "quadratic2_fixed", "quadratic2", [0.4, -1.0], [0.4, -2.0], fixed=(0,)
        ),
    ]
    xs = np.linspace(-4.0, 4.0, 200)
    ys = 8 * np.exp(-0.5 * ((xs - 0.25) / 0.9) ** 2) + rng.normal(0, 0.03, 200)
    result.append(
        make_case(
            "gaussian200",
            "gaussian",
            [7.0, 0.1, 1.1],
            [8.0, 0.25, 0.9],
            limits=[[0.01, None], [-2.0, 2.0], [0.1, 3.0]],
            data=[[float(x), float(y), 0.2] for x, y in zip(xs, ys)],
        )
    )
    xs = np.linspace(-4.95, 4.95, 100)
    shape = 3 * np.exp(-0.5 * ((xs - 0.3) / 0.8) ** 2) + np.exp(-0.2 * (xs + 5))
    for total in (10000, 1000000):
        scale = total / float(sum(shape))
        truth = [3 * scale, 0.3, 0.8, scale, 0.2]
        result.append(
            make_case(
                f"poisson{total}",
                "poisson",
                [2.4 * scale, 0.45, 0.88, 1.1 * scale, 0.18],
                truth,
                steps=[0.1 * scale, 0.1, 0.1, 0.1 * scale, 0.02],
                limits=[
                    [0.001, None],
                    [-2.0, 2.0],
                    [0.1, 3.0],
                    [0.001, None],
                    [0.01, 1.0],
                ],
                errordef=0.5,
                data=[[float(x), float(y * scale)] for x, y in zip(xs, shape)],
            )
        )
    return tuple(result)


@lru_cache(maxsize=1)
def bounds_cases():
    result = []
    # One-sided limits use five units as their distance scale.
    for variant, limit, edge, direction, target in (
        ("lower", [0.0, None], 0.0, 1, 1.0),
        ("upper", [None, 5.0], 5.0, -1, 1.0),
        ("both", [0.0, 5.0], 0.0, 1, 1.0),
        ("at_limit", [0.0, None], 0.0, 1, -1.0),
    ):
        for distance in (0.0, 1e-9, 1e-6, 1e-3, 1e-1):
            result.append(
                make_case(
                    f"{variant}_{distance!r}",
                    "bound",
                    [edge + direction * 5 * distance],
                    [max(0.0, target)],
                    limits=[limit],
                    target=target,
                    distance=distance,
                )
            )
    for y in (-1.0, -1.2):
        result.append(
            make_case(
                f"quadratic2_y{y!r}",
                "quadratic2",
                [0.4, y],
                [1.0, -2.0],
                limits=[[0.0, 2.0], [-3.0, -1.0]],
            )
        )
    return tuple(result)
