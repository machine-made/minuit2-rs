#!/usr/bin/env python3
"""
Generate a claim-oriented legacy->Rust traceability matrix from parity symbols.

Inputs:
  - reports/parity/functions.csv
  - verification/traceability/waivers.csv (optional)
  - verification/traceability/waiver_rules.csv (optional)

Outputs:
  - reports/verification/traceability_matrix.csv
  - reports/verification/traceability_summary.md
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PARITY_CSV = REPO_ROOT / "reports" / "parity" / "functions.csv"
DEFAULT_WAIVERS_CSV = REPO_ROOT / "verification" / "traceability" / "waivers.csv"
DEFAULT_WAIVER_RULES_CSV = REPO_ROOT / "verification" / "traceability" / "waiver_rules.csv"
DEFAULT_OUT_CSV = REPO_ROOT / "reports" / "verification" / "traceability_matrix.csv"
DEFAULT_OUT_MD = REPO_ROOT / "reports" / "verification" / "traceability_summary.md"

WAIVER_TYPES_RESOLVE = {
    "intentional",
    "out-of-scope",
    "tooling",
    "architectural",
    "upstream-removed",
    "api-shape-drift",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate legacy->Rust traceability matrix")
    parser.add_argument("--parity-csv", default=str(DEFAULT_PARITY_CSV))
    parser.add_argument("--waivers-csv", default=str(DEFAULT_WAIVERS_CSV))
    parser.add_argument("--waiver-rules-csv", default=str(DEFAULT_WAIVER_RULES_CSV))
    parser.add_argument("--out-csv", default=str(DEFAULT_OUT_CSV))
    parser.add_argument("--out-md", default=str(DEFAULT_OUT_MD))
    return parser.parse_args()


def read_waivers(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}

    out: dict[str, dict[str, str]] = {}
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        required = {"legacy_id", "waiver_type", "reason"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"waiver file missing columns: {sorted(missing)}")
        for row in reader:
            legacy_id = row["legacy_id"].strip()
            if not legacy_id:
                continue
            if legacy_id in out:
                raise ValueError(f"duplicate waiver for legacy_id: {legacy_id}")
            out[legacy_id] = {
                "waiver_type": row["waiver_type"].strip() or "unspecified",
                "reason": row["reason"].strip() or "unspecified",
                "source": "explicit",
            }
    return out


def read_waiver_rules(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []

    out: list[dict[str, str]] = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        required = {
            "raw_status",
            "rationale_contains",
            "upstream_file_regex",
            "upstream_symbol_regex",
            "waiver_type",
            "reason",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"waiver rules file missing columns: {sorted(missing)}")
        for row in reader:
            out.append(
                {
                    "raw_status": row["raw_status"].strip(),
                    "rationale_contains": row["rationale_contains"].strip(),
                    "upstream_file_regex": row["upstream_file_regex"].strip(),
                    "upstream_symbol_regex": row["upstream_symbol_regex"].strip(),
                    "waiver_type": row["waiver_type"].strip(),
                    "reason": row["reason"].strip(),
                }
            )
    return out


def mk_legacy_id(upstream_file: str, upstream_symbol: str, upstream_line: str) -> str:
    line = upstream_line.strip() if upstream_line is not None else ""
    line_part = line if line else "na"
    return f"{upstream_file}::{upstream_symbol}@{line_part}"


def classify_ambiguous_as_implemented(row: dict[str, str]) -> bool:
    return (
        row.get("status") == "needs-review"
        and row.get("rationale") == "multiple Rust symbol candidates"
        and bool(row.get("rust_file", "").strip())
        and bool(row.get("rust_symbol", "").strip())
        and row.get("upstream_symbol") != "<no_symbol_extracted>"
    )


def map_effective_status(
    raw_status: str,
    waiver_type: str | None,
    ambiguous_implemented: bool,
) -> str:
    if ambiguous_implemented:
        return "implemented"
    if raw_status == "implemented":
        return "implemented"
    if waiver_type and waiver_type in WAIVER_TYPES_RESOLVE:
        return "waived"
    if raw_status == "intentionally-skipped":
        return "waived"
    return "unresolved"


def apply_rule_waiver(row: dict[str, str], rules: list[dict[str, str]]) -> dict[str, str] | None:
    import re

    for rule in rules:
        if rule["raw_status"] and rule["raw_status"] != row["status"]:
            continue
        if rule["rationale_contains"] and rule["rationale_contains"] not in row.get("rationale", ""):
            continue
        if rule["upstream_file_regex"] and not re.search(rule["upstream_file_regex"], row.get("upstream_file", "")):
            continue
        if rule["upstream_symbol_regex"] and not re.search(rule["upstream_symbol_regex"], row.get("upstream_symbol", "")):
            continue
        return {
            "waiver_type": rule["waiver_type"],
            "reason": rule["reason"],
            "source": "rule",
        }
    return None


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("no rows to write")
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def render_summary(rows: list[dict[str, str]], waivers_path: Path, waiver_rules_path: Path) -> str:
    by_effective = Counter(row["effective_status"] for row in rows)
    by_raw = Counter(row["raw_status"] for row in rows)

    unresolved_by_file: dict[str, int] = defaultdict(int)
    for row in rows:
        if row["effective_status"] == "unresolved":
            unresolved_by_file[row["upstream_file"]] += 1

    top_unresolved = sorted(unresolved_by_file.items(), key=lambda x: (-x[1], x[0]))[:20]

    lines: list[str] = []
    lines.append("# Traceability Summary")
    lines.append("")
    lines.append("Source parity file: `reports/parity/functions.csv`")
    lines.append(f"Waivers file: `{waivers_path.relative_to(REPO_ROOT)}` (optional)")
    lines.append(f"Waiver rules file: `{waiver_rules_path.relative_to(REPO_ROOT)}` (optional)")
    lines.append("")
    lines.append("## Effective Status Counts")
    lines.append("")
    lines.append(f"- implemented: **{by_effective['implemented']}**")
    lines.append(f"- waived: **{by_effective['waived']}**")
    lines.append(f"- unresolved: **{by_effective['unresolved']}**")
    lines.append("")
    lines.append("## Raw Status Counts")
    lines.append("")
    for key in ["implemented", "missing", "needs-review", "intentionally-skipped"]:
        lines.append(f"- {key}: **{by_raw[key]}**")
    lines.append("")
    lines.append("## Top Unresolved Files")
    lines.append("")
    if top_unresolved:
        for path, count in top_unresolved:
            lines.append(f"- `{path}`: {count}")
    else:
        lines.append("- none")
    lines.append("")
    lines.append("## Gate Hint")
    lines.append("")
    if by_effective["unresolved"] == 0:
        lines.append("- Strict gate (`unresolved == 0`) is currently satisfied.")
    else:
        lines.append("- Strict gate is **not** satisfied; use non-regression gate in CI.")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    parity_csv = Path(args.parity_csv)
    waivers_csv = Path(args.waivers_csv)
    waiver_rules_csv = Path(args.waiver_rules_csv)
    out_csv = Path(args.out_csv)
    out_md = Path(args.out_md)

    if not parity_csv.exists():
        raise FileNotFoundError(f"parity csv not found: {parity_csv}")

    waivers = read_waivers(waivers_csv)
    waiver_rules = read_waiver_rules(waiver_rules_csv)

    rows_out: list[dict[str, str]] = []
    seen_legacy_ids: set[str] = set()

    with parity_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        required = {
            "upstream_repo",
            "upstream_subdir",
            "upstream_ref",
            "upstream_commit",
            "upstream_file",
            "upstream_symbol",
            "upstream_line",
            "rust_file",
            "rust_symbol",
            "rust_line",
            "status",
            "rationale",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"parity csv missing columns: {sorted(missing)}")

        for row in reader:
            legacy_id = mk_legacy_id(
                upstream_file=row["upstream_file"],
                upstream_symbol=row["upstream_symbol"],
                upstream_line=row["upstream_line"],
            )
            if legacy_id in seen_legacy_ids:
                raise ValueError(f"duplicate legacy id in parity csv: {legacy_id}")
            seen_legacy_ids.add(legacy_id)

            explicit_waiver = waivers.get(legacy_id)
            implicit_waiver = apply_rule_waiver(row, waiver_rules)
            waiver = explicit_waiver or implicit_waiver
            waiver_type = (waiver or {}).get("waiver_type")
            raw_status = row["status"]
            ambiguous_implemented = classify_ambiguous_as_implemented(row)
            effective_status = map_effective_status(raw_status, waiver_type, ambiguous_implemented)

            rows_out.append(
                {
                    "legacy_id": legacy_id,
                    "upstream_repo": row["upstream_repo"],
                    "upstream_subdir": row["upstream_subdir"],
                    "upstream_ref": row["upstream_ref"],
                    "upstream_commit": row["upstream_commit"],
                    "upstream_file": row["upstream_file"],
                    "upstream_symbol": row["upstream_symbol"],
                    "upstream_line": row["upstream_line"],
                    "rust_file": row["rust_file"],
                    "rust_symbol": row["rust_symbol"],
                    "rust_line": row["rust_line"],
                    "raw_status": raw_status,
                    "effective_status": effective_status,
                    "waived": "true" if effective_status == "waived" else "false",
                    "waiver_type": waiver_type or ("intentional" if raw_status == "intentionally-skipped" else ""),
                    "waiver_reason": (waiver or {}).get(
                        "reason",
                        "constructor/destructor/operator handled idiomatically in Rust"
                        if raw_status == "intentionally-skipped"
                        else "",
                    ),
                    "waiver_source": (waiver or {}).get("source", "auto-intentional" if raw_status == "intentionally-skipped" else ""),
                    "ambiguous_implemented": "true" if ambiguous_implemented else "false",
                    "rationale": row["rationale"],
                }
            )

    # Validate that all waivers point to known legacy IDs.
    unknown_waivers = sorted(set(waivers.keys()) - seen_legacy_ids)
    if unknown_waivers:
        sample = ", ".join(unknown_waivers[:5])
        raise ValueError(f"waiver legacy_id not found in parity matrix: {sample}")

    write_csv(out_csv, rows_out)
    summary = render_summary(rows_out, waivers_csv, waiver_rules_csv)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text(summary + "\n")

    effective = Counter(row["effective_status"] for row in rows_out)
    print(f"Wrote {out_csv.relative_to(REPO_ROOT)}")
    print(f"Wrote {out_md.relative_to(REPO_ROOT)}")
    print(
        "Effective status counts: "
        f"implemented={effective['implemented']} "
        f"waived={effective['waived']} "
        f"unresolved={effective['unresolved']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
