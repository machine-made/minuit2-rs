#!/usr/bin/env python3
"""
Traceability gate checks for CI.

Modes:
- strict: fail if any unresolved rows remain
- non-regression: compare current matrix vs baseline matrix and fail on regressions
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MATRIX = REPO_ROOT / "reports" / "verification" / "traceability_matrix.csv"
DEFAULT_BASELINE = REPO_ROOT / "verification" / "traceability" / "traceability_baseline.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check legacy->Rust traceability gate")
    parser.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    parser.add_argument("--baseline", default=str(DEFAULT_BASELINE))
    parser.add_argument("--mode", choices=["strict", "non-regression"], default="non-regression")
    return parser.parse_args()


def read_matrix(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"matrix not found: {path}")
    out: dict[str, dict[str, str]] = {}
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        required = {"legacy_id", "effective_status"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"matrix missing columns: {sorted(missing)}")
        for row in reader:
            legacy_id = row["legacy_id"]
            if legacy_id in out:
                raise ValueError(f"duplicate legacy_id in matrix: {legacy_id}")
            out[legacy_id] = row
    return out


def unresolved_ids(rows: dict[str, dict[str, str]]) -> set[str]:
    return {
        legacy_id
        for legacy_id, row in rows.items()
        if row.get("effective_status", "").strip() == "unresolved"
    }


def implemented_or_waived_ids(rows: dict[str, dict[str, str]]) -> set[str]:
    return {
        legacy_id
        for legacy_id, row in rows.items()
        if row.get("effective_status", "").strip() in {"implemented", "waived"}
    }


def summarize(rows: dict[str, dict[str, str]]) -> Counter[str]:
    return Counter(row.get("effective_status", "unknown") for row in rows.values())


def run_strict(matrix: dict[str, dict[str, str]]) -> int:
    counts = summarize(matrix)
    unresolved = counts["unresolved"]
    print(
        "Traceability summary: "
        f"implemented={counts['implemented']} "
        f"waived={counts['waived']} "
        f"unresolved={unresolved}"
    )
    if unresolved > 0:
        print("FAIL: strict traceability gate requires unresolved == 0")
        return 1
    print("PASS: strict traceability gate")
    return 0


def run_non_regression(
    matrix: dict[str, dict[str, str]],
    baseline: dict[str, dict[str, str]],
) -> int:
    current_ids = set(matrix.keys())
    baseline_ids = set(baseline.keys())

    missing_from_current = sorted(baseline_ids - current_ids)
    if missing_from_current:
        print(f"FAIL: {len(missing_from_current)} baseline legacy IDs missing in current matrix")
        print(f"Sample: {', '.join(missing_from_current[:5])}")
        return 1

    baseline_unresolved = unresolved_ids(baseline)
    current_unresolved = unresolved_ids(matrix)

    new_unresolved = sorted(current_unresolved - baseline_unresolved)
    regressed = sorted(
        legacy_id
        for legacy_id in (implemented_or_waived_ids(baseline) & current_ids)
        if matrix[legacy_id].get("effective_status", "").strip() == "unresolved"
    )

    counts_cur = summarize(matrix)
    counts_base = summarize(baseline)
    print(
        "Current summary: "
        f"implemented={counts_cur['implemented']} "
        f"waived={counts_cur['waived']} "
        f"unresolved={counts_cur['unresolved']}"
    )
    print(
        "Baseline summary: "
        f"implemented={counts_base['implemented']} "
        f"waived={counts_base['waived']} "
        f"unresolved={counts_base['unresolved']}"
    )

    failed = False
    if new_unresolved:
        failed = True
        print(f"FAIL: {len(new_unresolved)} new unresolved legacy IDs introduced")
        print(f"Sample: {', '.join(new_unresolved[:5])}")

    if regressed:
        failed = True
        print(f"FAIL: {len(regressed)} IDs regressed from implemented/waived to unresolved")
        print(f"Sample: {', '.join(regressed[:5])}")

    if counts_cur["unresolved"] > counts_base["unresolved"]:
        failed = True
        print(
            "FAIL: unresolved count increased "
            f"({counts_base['unresolved']} -> {counts_cur['unresolved']})"
        )

    if failed:
        return 1

    improved = counts_base["unresolved"] - counts_cur["unresolved"]
    print(f"PASS: non-regression traceability gate (unresolved improvement: {improved})")
    return 0


def main() -> int:
    args = parse_args()
    matrix = read_matrix(Path(args.matrix))

    if args.mode == "strict":
        return run_strict(matrix)

    baseline = read_matrix(Path(args.baseline))
    return run_non_regression(matrix, baseline)


if __name__ == "__main__":
    raise SystemExit(main())
