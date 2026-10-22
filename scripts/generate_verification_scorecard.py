#!/usr/bin/env python3
"""
Generate a claim-oriented verification scorecard for minuit2-rs.

Inputs:
  - verification/workloads/root_minuit2_v6_36_08.json
  - reports/verification/diff_results.csv
  - reports/verification/traceability_matrix.csv
  - reports/coverage/core_coverage_raw.txt
  - reports/coverage/all_features_coverage_raw.txt
  - reports/benchmarks/default_raw.txt
  - reports/benchmarks/parallel_raw.txt
  - reports/verification/root_test_port_status.md

Outputs:
  - reports/verification/manifest.json
  - reports/verification/scorecard.md
"""

from __future__ import annotations

import csv
import json
import re
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WORKLOADS_PATH = REPO_ROOT / "verification" / "workloads" / "root_minuit2_v6_36_08.json"
DIFF_RESULTS_PATH = REPO_ROOT / "reports" / "verification" / "diff_results.csv"
TRACEABILITY_PATH = REPO_ROOT / "reports" / "verification" / "traceability_matrix.csv"
ROOT_TEST_STATUS_PATH = REPO_ROOT / "reports" / "verification" / "root_test_port_status.md"
CORE_COVERAGE_RAW_PATH = REPO_ROOT / "reports" / "coverage" / "core_coverage_raw.txt"
ALL_COVERAGE_RAW_PATH = REPO_ROOT / "reports" / "coverage" / "all_features_coverage_raw.txt"
BENCH_DEFAULT_RAW_PATH = REPO_ROOT / "reports" / "benchmarks" / "default_raw.txt"
BENCH_PARALLEL_RAW_PATH = REPO_ROOT / "reports" / "benchmarks" / "parallel_raw.txt"
REF_COVERAGE_MANIFEST_PATH = REPO_ROOT / "reports" / "verification" / "reference_coverage" / "manifest.json"
EXEC_SURFACE_MANIFEST_PATH = REPO_ROOT / "reports" / "verification" / "executed_surface_manifest.json"
OUT_MANIFEST_PATH = REPO_ROOT / "reports" / "verification" / "manifest.json"
OUT_SCORECARD_PATH = REPO_ROOT / "reports" / "verification" / "scorecard.md"


TOTAL_RE = re.compile(
    r"^TOTAL\s+\d+\s+\d+\s+(?P<regions>[0-9.]+)%\s+"
    r"\d+\s+\d+\s+(?P<functions>[0-9.]+)%\s+"
    r"\d+\s+\d+\s+(?P<lines>[0-9.]+)%\s+"
)

BENCH_BLOCK_RE = re.compile(
    r"(?m)^(?P<name>[^\n].+?)\n"
    r"\s*time:\s+\["
    r"(?P<low>[0-9.]+)\s*(?P<low_unit>ns|µs|ms|s)\s+"
    r"(?P<mid>[0-9.]+)\s*(?P<mid_unit>ns|µs|ms|s)\s+"
    r"(?P<high>[0-9.]+)\s*(?P<high_unit>ns|µs|ms|s)\]\s*$"
)

UNIT_TO_US = {
    "ns": 0.001,
    "µs": 1.0,
    "ms": 1000.0,
    "s": 1_000_000.0,
}


def run_cmd(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, cwd=REPO_ROOT, text=True).strip()


def to_us(value: float, unit: str) -> float:
    return value * UNIT_TO_US[unit]


def parse_diff_results(path: Path) -> tuple[dict[str, int], int]:
    counts = Counter()
    total = 0
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            status = row.get("status", "").strip()
            if status:
                counts[status] += 1
            total += 1
    return {"pass": counts["pass"], "warn": counts["warn"], "fail": counts["fail"]}, total


def parse_traceability(path: Path) -> tuple[dict[str, int], int, list[tuple[str, int]]]:
    counts = Counter()
    unresolved_by_file: dict[str, int] = defaultdict(int)
    total = 0
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            effective = row.get("effective_status", "").strip()
            if effective:
                counts[effective] += 1
            if effective == "unresolved":
                unresolved_by_file[row.get("upstream_file", "")] += 1
            total += 1
    top_unresolved = sorted(unresolved_by_file.items(), key=lambda x: (-x[1], x[0]))[:10]
    return (
        {
            "implemented": counts["implemented"],
            "waived": counts["waived"],
            "unresolved": counts["unresolved"],
        },
        total,
        top_unresolved,
    )


def parse_total_coverage(path: Path) -> dict[str, float]:
    for line in path.read_text().splitlines():
        m = TOTAL_RE.match(line.strip())
        if m:
            return {
                "regions_pct": float(m.group("regions")),
                "functions_pct": float(m.group("functions")),
                "lines_pct": float(m.group("lines")),
            }
    raise ValueError(f"TOTAL coverage row not found in {path}")


def parse_benchmarks(path: Path) -> dict[str, tuple[float, str]]:
    raw = path.read_text()
    out: dict[str, tuple[float, str]] = {}
    for m in BENCH_BLOCK_RE.finditer(raw):
        name = m.group("name").strip()
        median = float(m.group("mid"))
        unit = m.group("mid_unit")
        out[name] = (to_us(median, unit), f"{median:.4g} {unit}")
    return out


def parse_p0_gap_count(path: Path) -> int:
    lines = path.read_text().splitlines()
    in_section = False
    count = 0
    for line in lines:
        if line.startswith("## ") and "Known P0 Gap" not in line:
            in_section = False
        if line.startswith("## Known P0 Gap"):
            in_section = True
            continue
        if not in_section:
            continue
        stripped = line.strip()
        if not stripped.startswith("|"):
            continue
        if "---" in stripped or "ROOT test intent" in stripped:
            continue
        count += 1
    return count


def parse_reference_coverage_manifest(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    data = json.loads(path.read_text())
    counts = data.get("counts") or {}
    return {
        "functions_in_scope": int(counts.get("functions_in_scope", 0)),
        "functions_executed": int(counts.get("functions_executed", 0)),
        "function_coverage_pct": float(counts.get("function_coverage_pct", 0.0)),
        "line_coverage_pct": float(counts.get("line_coverage_pct", 0.0)),
    }


def parse_executed_surface_manifest(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    data = json.loads(path.read_text())
    priority = data.get("priority_counts") or {}
    gate = data.get("gate") or {}
    return {
        "executed_functions_total": int(data.get("executed_functions_total", 0)),
        "mapped_implemented_total": int(data.get("mapped_implemented_total", 0)),
        "unmapped_total": int(data.get("unmapped_total", 0)),
        "priority_counts": {
            "P0": int(priority.get("P0", 0)),
            "P1": int(priority.get("P1", 0)),
            "P2": int(priority.get("P2", 0)),
        },
        "gate_pass": bool(gate.get("pass", False)),
    }


def yes_no(flag: bool) -> str:
    return "YES" if flag else "NO"


def main() -> int:
    workloads = json.loads(WORKLOADS_PATH.read_text())
    diff_counts, diff_total = parse_diff_results(DIFF_RESULTS_PATH)
    trace_counts, trace_total, top_unresolved = parse_traceability(TRACEABILITY_PATH)
    core_cov = parse_total_coverage(CORE_COVERAGE_RAW_PATH)
    all_cov = parse_total_coverage(ALL_COVERAGE_RAW_PATH)
    p0_gap_count = parse_p0_gap_count(ROOT_TEST_STATUS_PATH)
    ref_cov = parse_reference_coverage_manifest(REF_COVERAGE_MANIFEST_PATH)
    exec_surface = parse_executed_surface_manifest(EXEC_SURFACE_MANIFEST_PATH)

    bench_default = parse_benchmarks(BENCH_DEFAULT_RAW_PATH)
    bench_parallel = parse_benchmarks(BENCH_PARALLEL_RAW_PATH)
    scan_serial_name = "Quadratic 2D: MnScan serial (101 points)"
    scan_parallel_name = "Quadratic 2D: MnScan parallel (101 points)"
    scan_serial_default = bench_default.get(scan_serial_name)
    scan_serial_parallel = bench_parallel.get(scan_serial_name)
    scan_parallel_parallel = bench_parallel.get(scan_parallel_name)
    scan_speedup = None
    if scan_serial_parallel and scan_parallel_parallel and scan_parallel_parallel[0] > 0.0:
        scan_speedup = scan_serial_parallel[0] / scan_parallel_parallel[0]

    differential_gate = diff_counts["fail"] == 0
    traceability_gate = trace_counts["unresolved"] == 0
    root_p0_gate = p0_gap_count == 0
    executed_surface_gate = bool(exec_surface and exec_surface["gate_pass"])
    full_1to1_gate = differential_gate and traceability_gate and root_p0_gate and executed_surface_gate
    full_100_verifiable = full_1to1_gate and diff_counts["warn"] == 0

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z"),
        "project": {
            "name": "minuit2-rs",
            "git_branch": run_cmd(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
            "git_commit": run_cmd(["git", "rev-parse", "HEAD"]),
            "git_dirty": bool(run_cmd(["git", "status", "--porcelain"])),
        },
        "environment": {
            "rustc": run_cmd(["rustc", "--version"]),
            "cargo": run_cmd(["cargo", "--version"]),
            "python3": run_cmd(["python3", "--version"]),
        },
        "reference": workloads["reference"],
        "differential": {
            "workloads_total": diff_total,
            "status_counts": diff_counts,
            "gate_pass": differential_gate,
        },
        "traceability": {
            "symbols_total": trace_total,
            "status_counts": trace_counts,
            "gate_pass": traceability_gate,
        },
        "root_p0_tests": {
            "known_gap_count": p0_gap_count,
            "gate_pass": root_p0_gate,
        },
        "coverage": {
            "core_no_default_features": core_cov,
            "all_features": all_cov,
        },
        "reference_coverage": ref_cov,
        "executed_surface": exec_surface,
        "benchmarks": {
            "scan_serial_default": scan_serial_default[1] if scan_serial_default else None,
            "scan_serial_parallel_feature": scan_serial_parallel[1] if scan_serial_parallel else None,
            "scan_parallel_parallel_feature": scan_parallel_parallel[1] if scan_parallel_parallel else None,
            "scan_speedup_serial_over_parallel": scan_speedup,
        },
        "claims": {
            "numerical_parity_on_covered_workloads": differential_gate,
            "executed_surface_mapping_gate": executed_surface_gate,
            "full_1_to_1_functional_coverage_vs_reference": full_1to1_gate,
            "full_100_percent_verifiable_coverage": full_100_verifiable,
        },
    }

    OUT_MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")

    lines: list[str] = []
    lines.append("# Verification Scorecard (Claim-Oriented)")
    lines.append("")
    lines.append(f"- Generated at: `{manifest['generated_at_utc']}`")
    lines.append(f"- Rust commit: `{manifest['project']['git_commit']}`")
    lines.append(f"- Reference repo: `{manifest['reference']['repo']}`")
    lines.append(f"- Reference subtree: `{manifest['reference']['subdir']}`")
    lines.append(f"- Reference tag: `{manifest['reference']['tag']}`")
    lines.append(f"- Reference commit: `{manifest['reference']['commit']}`")
    lines.append("")
    lines.append("## Evidence Snapshot")
    lines.append("")
    lines.append(
        "- Differential workloads: "
        f"**{diff_total}** (pass={diff_counts['pass']}, warn={diff_counts['warn']}, fail={diff_counts['fail']})"
    )
    lines.append(
        "- Traceability symbols: "
        f"**{trace_total}** (implemented={trace_counts['implemented']}, "
        f"waived={trace_counts['waived']}, unresolved={trace_counts['unresolved']})"
    )
    lines.append(f"- Known ROOT P0 gaps: **{p0_gap_count}**")
    lines.append(
        "- Coverage (line): "
        f"`--no-default-features` **{core_cov['lines_pct']:.2f}%**, "
        f"`--all-features` **{all_cov['lines_pct']:.2f}%**"
    )
    if ref_cov:
        lines.append(
            "- Reference C++ executed-surface: "
            f"{ref_cov['functions_executed']}/{ref_cov['functions_in_scope']} functions "
            f"(**{ref_cov['function_coverage_pct']:.2f}%**) across `math/minuit2`"
        )
    if exec_surface:
        pc = exec_surface["priority_counts"]
        lines.append(
            "- Executed-surface unmapped gaps: "
            f"P0={pc['P0']}, P1={pc['P1']}, P2={pc['P2']} "
            f"(gate={yes_no(bool(exec_surface['gate_pass']))})"
        )
    if scan_serial_default:
        lines.append(f"- Benchmark serial scan (`default`): **{scan_serial_default[1]}**")
    if scan_serial_parallel and scan_parallel_parallel:
        lines.append(
            "- Benchmark scan in `parallel` feature run: "
            f"serial={scan_serial_parallel[1]}, parallel={scan_parallel_parallel[1]}"
        )
    if scan_speedup is not None:
        lines.append(f"- Scan speedup (`serial/parallel`, parallel feature run): **{scan_speedup:.2f}x**")
    lines.append("")
    lines.append("## Claim Gates")
    lines.append("")
    lines.append("| Claim | Gate | Status |")
    lines.append("|---|---|---|")
    lines.append(
        "| Numerical parity on covered workloads | `diff fail == 0` | "
        f"**{yes_no(differential_gate)}** |"
    )
    lines.append(
        "| Full symbol/function traceability | `traceability unresolved == 0` | "
        f"**{yes_no(traceability_gate)}** |"
    )
    lines.append(
        "| ROOT P0 regression parity completeness | `known P0 gaps == 0` | "
        f"**{yes_no(root_p0_gate)}** |"
    )
    lines.append(
        "| Executed-surface mapping completeness | `executed-surface P0/P1 == 0` | "
        f"**{yes_no(executed_surface_gate)}** |"
    )
    lines.append(
        "| Full 1:1 functional coverage claim | all above gates true | "
        f"**{yes_no(full_1to1_gate)}** |"
    )
    lines.append(
        "| Full 100% verifiable coverage claim | 1:1 gate + zero warnings | "
        f"**{yes_no(full_100_verifiable)}** |"
    )
    lines.append("")
    lines.append("## Blocking Gaps")
    lines.append("")
    if not traceability_gate:
        lines.append("- Unresolved traceability still present. Top unresolved upstream files:")
        for upstream_file, count in top_unresolved:
            lines.append(f"  - `{upstream_file}`: {count}")
    if not root_p0_gate:
        lines.append("- Known ROOT P0 test-port gap(s) remain; see `reports/verification/root_test_port_status.md`.")
    if not executed_surface_gate:
        lines.append("- Executed-surface mapping gate fails; see `reports/verification/executed_surface_mapping.md`.")
    if diff_counts["warn"] > 0:
        lines.append("- Differential warnings remain (NFCN divergence); see `reports/verification/diff_summary.md`.")
    if traceability_gate and root_p0_gate and executed_surface_gate and diff_counts["warn"] == 0:
        lines.append("- None.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("scripts/build_root_reference_runner.sh v6-36-08")
    lines.append("python3 scripts/compare_ref_vs_rust.py")
    lines.append("python3 scripts/generate_traceability_matrix.py")
    lines.append("python3 scripts/check_traceability_gate.py --mode non-regression")
    lines.append("python3 scripts/generate_executed_surface_mapping.py")
    lines.append("python3 scripts/check_executed_surface_gate.py --mode non-regression")
    lines.append("cargo test --no-default-features")
    lines.append("PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo test --all-features")
    lines.append("cargo llvm-cov --no-default-features --summary-only > reports/coverage/core_coverage_raw.txt")
    lines.append("PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo llvm-cov --all-features --summary-only > reports/coverage/all_features_coverage_raw.txt")
    lines.append("python3 scripts/generate_coverage_reports.py")
    lines.append("python3 scripts/generate_reference_coverage.py --root-tag v6-36-08")
    lines.append("cargo bench --bench benchmarks -- --noplot > reports/benchmarks/default_raw.txt")
    lines.append("cargo bench --features parallel --bench benchmarks -- --noplot > reports/benchmarks/parallel_raw.txt")
    lines.append("python3 scripts/generate_benchmark_report.py")
    lines.append("python3 scripts/generate_verification_scorecard.py")
    lines.append("```")
    lines.append("")

    OUT_SCORECARD_PATH.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT_MANIFEST_PATH.relative_to(REPO_ROOT)}")
    print(f"Wrote {OUT_SCORECARD_PATH.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
