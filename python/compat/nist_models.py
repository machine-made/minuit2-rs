"""Shared NIST StRD models/data for hard-dataset baseline scripts.

Parses committed examples/data/nist/*.dat files for Start 2 values, certified
parameters, certified residual SS, and observations. Model formulas are
transcribed from the NIST Model sections (also mirrored in
``tests/nist_strd_certified.rs`` tolerance scheme).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable
import math

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "examples" / "data" / "nist"


@dataclass(frozen=True)
class NistDataset:
    name: str
    x: tuple[float, ...]
    y: tuple[float, ...]
    start2: tuple[float, ...]
    certified: tuple[float, ...]
    certified_rss: float
    rel_tol: float
    model: Callable[[tuple[float, ...], float], float]


def _parse_float_token(token: str) -> float | None:
    t = token.strip().rstrip(",")
    if not t:
        return None
    try:
        return float(t)
    except ValueError:
        if t.startswith("."):
            return float("0" + t)
        if t.startswith("-."):
            return float("-0" + t[1:])
    return None


def _parse_floats(text: str) -> list[float]:
    return [v for tok in text.split() if (v := _parse_float_token(tok)) is not None]


def _parse_dat(name: str, nparams: int) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...], float]:
    path = DATA_DIR / f"{name}.dat"
    if not path.exists():
        raise FileNotFoundError(f"missing NIST data file: {path}")
    start2: list[float] = []
    certified: list[float] = []
    x: list[float] = []
    y: list[float] = []
    rss = math.nan
    in_data = False
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = raw.strip()
        if "=" in s:
            lhs, rhs = s.split("=", 1)
            lhs = lhs.strip()
            if lhs.startswith("b") and len(lhs) >= 2 and lhs[1].isdigit():
                nums = _parse_floats(rhs)
                if len(nums) >= 4:
                    start2.append(nums[1])
                    certified.append(nums[2])
        if s.startswith("Residual Sum of Squares:"):
            nums = _parse_floats(s)
            if nums:
                rss = nums[-1]
        if s.startswith("Data:") and s.removeprefix("Data:").lstrip().startswith("y"):
            in_data = True
            continue
        if in_data:
            nums = _parse_floats(s)
            if len(nums) >= 2:
                y.append(nums[0])
                x.append(nums[1])
    if len(start2) != nparams or len(certified) != nparams or not x or not math.isfinite(rss):
        raise ValueError(f"failed to parse {name}: start2={len(start2)} certified={len(certified)} points={len(x)} rss={rss}")
    return tuple(x), tuple(y), tuple(start2), tuple(certified), rss


def model_lanczos3(p: tuple[float, ...], x: float) -> float:
    return p[0] * math.exp(-p[1] * x) + p[2] * math.exp(-p[3] * x) + p[4] * math.exp(-p[5] * x)


def model_boxbod(p: tuple[float, ...], x: float) -> float:
    return p[0] * (1.0 - math.exp(-p[1] * x))


def model_mgh09(p: tuple[float, ...], x: float) -> float:
    den = x * x + x * p[2] + p[3]
    if abs(den) < 1e-300:
        return math.nan
    return p[0] * (x * x + x * p[1]) / den


def model_hahn1(p: tuple[float, ...], x: float) -> float:
    x2 = x * x
    x3 = x2 * x
    den = 1.0 + p[4] * x + p[5] * x2 + p[6] * x3
    if abs(den) < 1e-300:
        return math.nan
    return (p[0] + p[1] * x + p[2] * x2 + p[3] * x3) / den


SPECS = {
    "Lanczos3": (6, 1e-2, model_lanczos3),
    "BoxBOD": (2, 1e-2, model_boxbod),
    "MGH09": (4, 1e-2, model_mgh09),
    "Hahn1": (7, 1e-2, model_hahn1),
}


def load_dataset(name: str) -> NistDataset:
    nparams, rel_tol, model = SPECS[name]
    x, y, start2, certified, rss = _parse_dat(name, nparams)
    return NistDataset(name, x, y, start2, certified, rss, rel_tol, model)


def dataset_names() -> tuple[str, ...]:
    return tuple(SPECS)
