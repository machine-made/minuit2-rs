#!/usr/bin/env python3
"""
Generate markdown coverage reports from cargo-llvm-cov raw output logs.

Inputs:
  - reports/coverage/core_coverage_raw.txt
  - reports/coverage/all_features_coverage_raw.txt

Outputs:
  - reports/coverage/core_coverage.md
  - reports/coverage/all_features_coverage.md
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


ROW_RE = re.compile(
    r"^(?P<file>\S+)\s+"
    r"(?P<regions>\d+)\s+(?P<missed_regions>\d+)\s+(?P<regions_cov>[0-9.]+)%\s+"
    r"(?P<functions>\d+)\s+(?P<missed_functions>\d+)\s+(?P<functions_cov>[0-9.]+)%\s+"
    r"(?P<lines>\d+)\s+(?P<missed_lines>\d+)\s+(?P<lines_cov>[0-9.]+)%\s+"
    r"(?P<branches>\d+)\s+(?P<missed_branches>\d+)\s+(?P<branches_cov>[0-9.\-]+%|-)\s*$"
)


@dataclass
class CoverageRow:
    file: str
    regions_cov: float
    functions_cov: float
    lines_cov: float


def parse_rows(raw: str) -> list[CoverageRow]:
    rows: list[CoverageRow] = []
    for line in raw.splitlines():
        m = ROW_RE.match(line.strip())
        if not m:
            continue
        rows.append(
            CoverageRow(
                file=m.group("file"),
                regions_cov=float(m.group("regions_cov")),
                functions_cov=float(m.group("functions_cov")),
                lines_cov=float(m.group("lines_cov")),
            )
        )
    return rows


def render_report(title: str, command: str, rows: list[CoverageRow]) -> str:
    if not rows:
        raise ValueError("No coverage rows parsed from raw log")

    total = next((r for r in rows if r.file == "TOTAL"), None)
    if total is None:
        raise ValueError("TOTAL row missing from coverage log")

    module_rows = [r for r in rows if r.file != "TOTAL"]
    low_rows = sorted(module_rows, key=lambda r: r.lines_cov)[:12]

    lines = []
    lines.append(f"# {title}")
    lines.append("")
    lines.append(f"- Command: `{command}`")
    lines.append("")
    lines.append("## Totals")
    lines.append("")
    lines.append(f"- Regions: **{total.regions_cov:.2f}%**")
    lines.append(f"- Functions: **{total.functions_cov:.2f}%**")
    lines.append(f"- Lines: **{total.lines_cov:.2f}%**")
    lines.append("")
    lines.append("## Lowest Line Coverage Modules")
    lines.append("")
    lines.append("| Module | Line Coverage | Function Coverage | Region Coverage |")
    lines.append("|---|---:|---:|---:|")
    for row in low_rows:
        lines.append(
            f"| `{row.file}` | {row.lines_cov:.2f}% | {row.functions_cov:.2f}% | {row.regions_cov:.2f}% |"
        )
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    coverage_dir = root / "reports" / "coverage"
    coverage_dir.mkdir(parents=True, exist_ok=True)

    core_raw_path = coverage_dir / "core_coverage_raw.txt"
    all_raw_path = coverage_dir / "all_features_coverage_raw.txt"

    core_rows = parse_rows(core_raw_path.read_text())
    all_rows = parse_rows(all_raw_path.read_text())

    core_report = render_report(
        title="Core Coverage Report (`--no-default-features`)",
        command="cargo llvm-cov --no-default-features --summary-only",
        rows=core_rows,
    )
    all_report = render_report(
        title="All Features Coverage Report (`--all-features`)",
        command="PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo llvm-cov --all-features --summary-only",
        rows=all_rows,
    )

    (coverage_dir / "core_coverage.md").write_text(core_report + "\n")
    (coverage_dir / "all_features_coverage.md").write_text(all_report + "\n")


if __name__ == "__main__":
    main()
