#!/usr/bin/env python3
"""
Executed-surface gate checks for CI.

Modes:
- strict: fail if any P0/P1 unmapped executed gaps remain
- non-regression: compare current high-priority gaps against a baseline and fail on regressions
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GAPS = REPO_ROOT / "reports" / "verification" / "executed_surface_gaps.csv"
DEFAULT_BASELINE = REPO_ROOT / "verification" / "traceability" / "executed_surface_gaps_baseline.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check executed-surface mapping gate")
    parser.add_argument("--gaps", default=str(DEFAULT_GAPS))
    parser.add_argument("--baseline", default=str(DEFAULT_BASELINE))
    parser.add_argument("--mode", choices=["strict", "non-regression"], default="non-regression")
    return parser.parse_args()


def read_gaps(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"gaps csv not found: {path}")
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    required = {"upstream_file", "upstream_symbol", "function_mangled", "gap_priority"}
    missing = required - set(rows[0].keys() if rows else required)
    if rows and missing:
        raise ValueError(f"gaps csv missing columns: {sorted(missing)}")
    return rows


def gap_id(row: dict[str, str]) -> str:
    return (
        f"{row.get('upstream_file', '').strip()}::"
        f"{row.get('upstream_symbol', '').strip()}::"
        f"{row.get('function_mangled', '').strip()}"
    )


def high_priority_ids(rows: list[dict[str, str]]) -> set[str]:
    return {
        gap_id(r)
        for r in rows
        if r.get("gap_priority", "").strip() in {"P0", "P1"}
    }


def priority_counts(rows: list[dict[str, str]]) -> Counter[str]:
    out = Counter()
    for row in rows:
        out[row.get("gap_priority", "").strip()] += 1
    return out


def run_strict(rows: list[dict[str, str]]) -> int:
    counts = priority_counts(rows)
    high = counts["P0"] + counts["P1"]
    print(
        "Executed-surface summary: "
        f"P0={counts['P0']} P1={counts['P1']} P2={counts['P2']}"
    )
    if high > 0:
        print("FAIL: strict executed-surface gate requires P0 == 0 and P1 == 0")
        return 1
    print("PASS: strict executed-surface gate")
    return 0


def run_non_regression(current_rows: list[dict[str, str]], baseline_rows: list[dict[str, str]]) -> int:
    cur_counts = priority_counts(current_rows)
    base_counts = priority_counts(baseline_rows)
    cur_ids = high_priority_ids(current_rows)
    base_ids = high_priority_ids(baseline_rows)

    print(
        "Current summary: "
        f"P0={cur_counts['P0']} P1={cur_counts['P1']} P2={cur_counts['P2']}"
    )
    print(
        "Baseline summary: "
        f"P0={base_counts['P0']} P1={base_counts['P1']} P2={base_counts['P2']}"
    )

    failed = False

    new_high = sorted(cur_ids - base_ids)
    if new_high:
        failed = True
        print(f"FAIL: {len(new_high)} new P0/P1 executed-surface gaps introduced")
        print(f"Sample: {', '.join(new_high[:5])}")

    if cur_counts["P0"] > base_counts["P0"]:
        failed = True
        print(f"FAIL: P0 gap count increased ({base_counts['P0']} -> {cur_counts['P0']})")

    if cur_counts["P1"] > base_counts["P1"]:
        failed = True
        print(f"FAIL: P1 gap count increased ({base_counts['P1']} -> {cur_counts['P1']})")

    if failed:
        return 1

    print("PASS: non-regression executed-surface gate")
    return 0


def main() -> int:
    args = parse_args()
    current_rows = read_gaps(Path(args.gaps))

    if args.mode == "strict":
        return run_strict(current_rows)

    baseline_rows = read_gaps(Path(args.baseline))
    return run_non_regression(current_rows, baseline_rows)


if __name__ == "__main__":
    raise SystemExit(main())
