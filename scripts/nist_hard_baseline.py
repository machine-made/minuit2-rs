#!/usr/bin/env python3
"""Baseline NIST hard StRD datasets: iminuit vs minuit2-rs from Start 2."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from python.compat.nist_models import dataset_names, load_dataset  # noqa: E402
from iminuit import Minuit as IMinuit  # noqa: E402
from minuit2 import Minuit as RMinuit  # noqa: E402

REPORT = ROOT / "reports" / "parity" / "nist_hard_baseline.md"


def _rss_fn(ds):
    def fcn(*pars: float) -> float:
        p = tuple(float(v) for v in pars)
        rss = 0.0
        for x, y in zip(ds.x, ds.y):
            pred = ds.model(p, x)
            if not math.isfinite(pred):
                return 1e30
            r = y - pred
            rss += r * r
        return rss

    return fcn


def _params_ok(vals: list[float], certified: tuple[float, ...], rel_tol: float) -> bool:
    if len(vals) != len(certified):
        return False
    for got, cert in zip(vals, certified):
        if not math.isfinite(got):
            return False
        if abs(got - cert) / max(abs(cert), 1e-300) > rel_tol:
            return False
    return True


def run_one(cls: Any, ds, strategy: int) -> dict[str, Any]:
    names = [f"b{i + 1}" for i in range(len(ds.start2))]
    m = cls(_rss_fn(ds), *ds.start2, name=tuple(names))
    m.errordef = 1.0
    m.strategy = strategy
    try:
        m.migrad(ncall=1_000_000)
    except TypeError:
        m.migrad()
    vals = [float(m.values[name]) for name in names]
    fval = float(m.fval)
    valid = bool(m.valid)
    nfcn = int(getattr(m, "nfcn", -1))
    params_ok = _params_ok(vals, ds.certified, ds.rel_tol)
    rss_ok = math.isfinite(fval) and abs(fval - ds.certified_rss) <= max(1e-8, abs(ds.certified_rss) * 1e-2)
    return {"valid": valid, "fval": fval, "rss_ok": rss_ok, "params_ok": params_ok, "nfcn": nfcn, "values": vals}


def label(r: dict[str, Any]) -> str:
    status = "OK" if r["valid"] and r["rss_ok"] and r["params_ok"] else "FAIL"
    return f"{status} (valid={r['valid']}, fval={r['fval']:.6g}, params={r['params_ok']}, nfcn={r['nfcn']})"


def run(datasets: list[str]) -> list[dict[str, Any]]:
    rows = []
    for name in datasets:
        ds = load_dataset(name)
        row = {"dataset": name, "certified_rss": ds.certified_rss, "tol": ds.rel_tol, "results": {}}
        for impl, cls in (("iminuit", IMinuit), ("minuit2-rs", RMinuit)):
            for strategy in (1, 2):
                row["results"][(impl, strategy)] = run_one(cls, ds, strategy)
        rows.append(row)
    return rows


def print_matrix(rows: list[dict[str, Any]]) -> None:
    print("dataset | iminuit s1 | iminuit s2 | minuit2-rs s1 | minuit2-rs s2")
    print("--- | --- | --- | --- | ---")
    for row in rows:
        r = row["results"]
        print(f"{row['dataset']} | {label(r[('iminuit', 1)])} | {label(r[('iminuit', 2)])} | {label(r[('minuit2-rs', 1)])} | {label(r[('minuit2-rs', 2)])}")


def succeeded(r: dict[str, Any]) -> bool:
    return bool(r["valid"] and r["rss_ok"] and r["params_ok"])


def write_report(rows: list[dict[str, Any]]) -> None:
    gaps = []
    for row in rows:
        r = row["results"]
        if (succeeded(r[("iminuit", 1)]) or succeeded(r[("iminuit", 2)])) and not (succeeded(r[("minuit2-rs", 1)]) or succeeded(r[("minuit2-rs", 2)])):
            gaps.append(row["dataset"])
    lines = ["# NIST hard-dataset baseline (iminuit vs minuit2-rs)", ""]
    lines.append("Priority genuine-gap targets: " + (", ".join(gaps) if gaps else "none"))
    lines += ["", "ROOT runner: skipped; upstream iminuit is the C++ Minuit2 comparator for this bead.", "", "| Dataset | iminuit s1 | iminuit s2 | minuit2-rs s1 | minuit2-rs s2 |", "| --- | --- | --- | --- | --- |"]
    for row in rows:
        r = row["results"]
        lines.append(f"| {row['dataset']} | {label(r[('iminuit', 1)])} | {label(r[('iminuit', 2)])} | {label(r[('minuit2-rs', 1)])} | {label(r[('minuit2-rs', 2)])} |")
    lines.append("")
    for row in rows:
        r = row["results"]
        i_ok = succeeded(r[("iminuit", 1)]) or succeeded(r[("iminuit", 2)])
        rust_ok = succeeded(r[("minuit2-rs", 1)]) or succeeded(r[("minuit2-rs", 2)])
        verdict = "genuine gap" if i_ok and not rust_ok else "parity-failure" if not i_ok and not rust_ok else "parity-success"
        lines.append(f"- **{row['dataset']}**: {verdict}. Certified residual SS={row['certified_rss']:.6g}; parameter tolerance={row['tol']:.0e}. iminuit {'reaches' if i_ok else 'does not reach'} the certified solution from Start 2; minuit2-rs {'reaches' if rust_ok else 'does not reach'} it under the same strategies.")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=dataset_names())
    ns = ap.parse_args()
    names = [ns.dataset] if ns.dataset else list(dataset_names())
    rows = run(names)
    print_matrix(rows)
    write_report(rows)
    print(f"\nWrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
