#!/usr/bin/env python3
"""
Run differential comparisons between ROOT Minuit2 reference runner and Rust runner.

Inputs:
  - verification/workloads/root_minuit2_v6_36_08.json

Outputs:
  - reports/verification/raw/ref/<workload>.json
  - reports/verification/raw/rust/<workload>.json
  - reports/verification/diff_results.csv
  - reports/verification/diff_summary.md
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKLOADS = REPO_ROOT / "verification" / "workloads" / "root_minuit2_v6_36_08.json"
DEFAULT_REF_BIN = REPO_ROOT / "third_party" / "root_ref_build" / "ref_runner" / "ref_runner"
DEFAULT_REPORT_DIR = REPO_ROOT / "reports" / "verification"


@dataclass
class DiffOutcome:
    status: str
    issues: list[str]
    warnings: list[str]
    max_param_abs: float
    max_error_abs: float
    max_cov_abs: float
    minos_abs: float
    fval_abs: float
    edm_abs: float
    nfcn_rel: float


def run_cmd(cmd: list[str]) -> str:
    out = subprocess.check_output(cmd, cwd=REPO_ROOT, text=True)
    return out.strip()


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def max_abs_diff(a: list[float], b: list[float]) -> float:
    if len(a) != len(b):
        return float("inf")
    if not a:
        return 0.0
    return max(abs(x - y) for x, y in zip(a, b))


def max_abs_diff_matrix(a: list[list[float]], b: list[list[float]]) -> float:
    if len(a) != len(b):
        return float("inf")
    if not a:
        return 0.0
    out = 0.0
    for ra, rb in zip(a, b):
        if len(ra) != len(rb):
            return float("inf")
        for xa, xb in zip(ra, rb):
            out = max(out, abs(xa - xb))
    return out


def relative_diff(a: float, b: float) -> float:
    denom = max(abs(a), abs(b), 1e-16)
    return abs(a - b) / denom


def compare(ref: dict[str, Any], rust: dict[str, Any], tol: dict[str, float]) -> DiffOutcome:
    issues: list[str] = []
    warnings: list[str] = []

    ref_valid = bool(ref.get("valid", False))
    rust_valid = bool(rust.get("valid", False))
    if ref_valid != rust_valid:
        issues.append(f"valid mismatch ref={ref_valid} rust={rust_valid}")

    fval_abs = abs(float(ref["fval"]) - float(rust["fval"]))
    if fval_abs > float(tol.get("fval_abs", 0.0)):
        issues.append(f"fval abs diff {fval_abs:.3e} > {tol.get('fval_abs')}")

    edm_abs = abs(float(ref["edm"]) - float(rust["edm"]))
    if edm_abs > float(tol.get("edm_abs", 0.0)):
        issues.append(f"edm abs diff {edm_abs:.3e} > {tol.get('edm_abs')}")

    params_ref = [float(x) for x in ref.get("params", [])]
    params_rust = [float(x) for x in rust.get("params", [])]
    max_param_abs = max_abs_diff(params_ref, params_rust)
    if not math.isfinite(max_param_abs):
        issues.append("parameter vector size mismatch")
    elif max_param_abs > float(tol.get("param_abs", 0.0)):
        issues.append(f"param max abs diff {max_param_abs:.3e} > {tol.get('param_abs')}")

    max_error_abs = 0.0
    if "error_abs" in tol:
        errors_ref = [float(x) for x in ref.get("errors", [])]
        errors_rust = [float(x) for x in rust.get("errors", [])]
        max_error_abs = max_abs_diff(errors_ref, errors_rust)
        if not math.isfinite(max_error_abs):
            issues.append("error vector size mismatch")
        elif max_error_abs > float(tol["error_abs"]):
            issues.append(f"error max abs diff {max_error_abs:.3e} > {tol['error_abs']}")

    max_cov_abs = 0.0
    if "cov_abs" in tol:
        ref_has_cov = bool(ref.get("has_covariance", False))
        rust_has_cov = bool(rust.get("has_covariance", False))
        if ref_has_cov != rust_has_cov:
            issues.append(f"covariance presence mismatch ref={ref_has_cov} rust={rust_has_cov}")
        elif ref_has_cov and rust_has_cov:
            cov_ref = ref.get("covariance", [])
            cov_rust = rust.get("covariance", [])
            max_cov_abs = max_abs_diff_matrix(cov_ref, cov_rust)
            cov_tol = float(tol["cov_abs"])
            if not math.isfinite(max_cov_abs) or max_cov_abs > cov_tol:
                issues.append(f"covariance max abs diff {max_cov_abs:.3e} > {cov_tol}")

    minos_abs = 0.0
    if "minos_abs" in tol:
        ref_has_minos = bool(ref.get("has_minos", False))
        rust_has_minos = bool(rust.get("has_minos", False))
        if ref_has_minos != rust_has_minos:
            issues.append(f"minos presence mismatch ref={ref_has_minos} rust={rust_has_minos}")
        elif ref_has_minos and rust_has_minos:
            minos_tol = float(tol["minos_abs"])
            ref_m = ref.get("minos") or {}
            rust_m = rust.get("minos") or {}
            if bool(ref_m.get("valid", False)) != bool(rust_m.get("valid", False)):
                issues.append("minos validity mismatch")
            lower_abs = abs(float(ref_m.get("lower", 0.0)) - float(rust_m.get("lower", 0.0)))
            upper_abs = abs(float(ref_m.get("upper", 0.0)) - float(rust_m.get("upper", 0.0)))
            minos_abs = max(lower_abs, upper_abs)
            if minos_abs > minos_tol:
                issues.append(f"minos max abs diff {minos_abs:.3e} > {minos_tol}")

    nfcn_ref = float(ref.get("nfcn", 0.0))
    nfcn_rust = float(rust.get("nfcn", 0.0))
    nfcn_rel = relative_diff(nfcn_ref, nfcn_rust)
    nfcn_rel_warn = float(tol.get("nfcn_rel_warn", 1.0))
    if nfcn_rel > nfcn_rel_warn:
        nfcn_waiver = tol.get("nfcn_rel_waiver")
        if isinstance(nfcn_waiver, str) and nfcn_waiver.strip():
            warnings.append(f"waived nfcn relative diff {nfcn_rel:.3f} > {nfcn_rel_warn}: {nfcn_waiver.strip()}")
        else:
            warnings.append(f"nfcn relative diff {nfcn_rel:.3f} > {nfcn_rel_warn}")

    effective_warnings = [w for w in warnings if not w.startswith("waived ")]

    status = "pass"
    if issues:
        status = "fail"
    elif effective_warnings:
        status = "warn"

    return DiffOutcome(
        status=status,
        issues=issues,
        warnings=warnings,
        max_param_abs=max_param_abs if math.isfinite(max_param_abs) else -1.0,
        max_error_abs=max_error_abs if math.isfinite(max_error_abs) else -1.0,
        max_cov_abs=max_cov_abs if math.isfinite(max_cov_abs) else -1.0,
        minos_abs=minos_abs,
        fval_abs=fval_abs,
        edm_abs=edm_abs,
        nfcn_rel=nfcn_rel,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare ROOT Minuit2 reference runner against Rust runner")
    parser.add_argument("--workloads", default=str(DEFAULT_WORKLOADS), help="Path to workload JSON")
    parser.add_argument("--ref-bin", default=str(DEFAULT_REF_BIN), help="Path to reference C++ runner")
    parser.add_argument(
        "--rust-cmd",
        default="cargo run --quiet --bin ref_compare_runner -- --workload",
        help="Rust command prefix; workload id is appended as final argument",
    )
    parser.add_argument("--report-dir", default=str(DEFAULT_REPORT_DIR), help="Output report directory")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workloads_path = Path(args.workloads)
    ref_bin = Path(args.ref_bin)
    report_dir = Path(args.report_dir)

    if not ref_bin.exists():
        raise FileNotFoundError(f"reference runner not found: {ref_bin}")

    spec = json.loads(workloads_path.read_text())
    workloads = spec["workloads"]

    raw_ref_dir = report_dir / "raw" / "ref"
    raw_rust_dir = report_dir / "raw" / "rust"
    raw_ref_dir.mkdir(parents=True, exist_ok=True)
    raw_rust_dir.mkdir(parents=True, exist_ok=True)

    rust_cmd_prefix = args.rust_cmd.split()
    rows: list[dict[str, str]] = []

    for w in workloads:
        wid = w["id"]

        ref_out = run_cmd([str(ref_bin), "--workload", wid])
        rust_out = run_cmd([*rust_cmd_prefix, wid])

        ref_json = json.loads(ref_out)
        rust_json = json.loads(rust_out)

        write_json(raw_ref_dir / f"{wid}.json", ref_json)
        write_json(raw_rust_dir / f"{wid}.json", rust_json)

        outcome = compare(ref_json, rust_json, w.get("tolerances", {}))

        rows.append(
            {
                "workload": wid,
                "status": outcome.status,
                "issues": " | ".join(outcome.issues),
                "warnings": " | ".join(outcome.warnings),
                "fval_abs": f"{outcome.fval_abs:.6e}",
                "edm_abs": f"{outcome.edm_abs:.6e}",
                "max_param_abs": f"{outcome.max_param_abs:.6e}",
                "max_error_abs": f"{outcome.max_error_abs:.6e}",
                "max_cov_abs": f"{outcome.max_cov_abs:.6e}",
                "minos_abs": f"{outcome.minos_abs:.6e}",
                "nfcn_rel": f"{outcome.nfcn_rel:.6e}",
            }
        )

    csv_path = report_dir / "diff_results.csv"
    report_dir.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            lineterminator="\n",
            fieldnames=[
                "workload",
                "status",
                "issues",
                "warnings",
                "fval_abs",
                "edm_abs",
                "max_param_abs",
                "max_error_abs",
                "max_cov_abs",
                "minos_abs",
                "nfcn_rel",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    counts = {"pass": 0, "warn": 0, "fail": 0}
    for row in rows:
        counts[row["status"]] += 1

    lines = [
        "# Differential Verification Summary",
        "",
        f"Reference repo: `{spec['reference']['repo']}`",
        f"Reference subtree: `{spec['reference']['subdir']}`",
        f"Reference tag: `{spec['reference']['tag']}`",
        f"Reference commit: `{spec['reference']['commit']}`",
        "",
        "## Status Counts",
        "",
        f"- pass: **{counts['pass']}**",
        f"- warn: **{counts['warn']}**",
        f"- fail: **{counts['fail']}**",
        "",
        "## Per-Workload Results",
        "",
        "| Workload | Status | Issues | Warnings |",
        "|---|---|---|---|",
    ]

    for row in rows:
        issues = row["issues"] if row["issues"] else "-"
        warnings = row["warnings"] if row["warnings"] else "-"
        lines.append(f"| `{row['workload']}` | `{row['status']}` | {issues} | {warnings} |")

    lines.extend(
        [
            "",
            "## Artifacts",
            "",
            "- `reports/verification/diff_results.csv`",
            "- `reports/verification/raw/ref/*.json`",
            "- `reports/verification/raw/rust/*.json`",
            "",
            "## Notes",
            "",
            "- `fail` means correctness metrics exceeded workload tolerances.",
            "- `warn` means correctness metrics passed, but NFCN divergence exceeded warning threshold.",
        ]
    )

    summary_path = report_dir / "diff_summary.md"
    summary_path.write_text("\n".join(lines) + "\n")

    print(f"Wrote {csv_path.relative_to(REPO_ROOT)}")
    print(f"Wrote {summary_path.relative_to(REPO_ROOT)}")
    print(f"Status counts: pass={counts['pass']} warn={counts['warn']} fail={counts['fail']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
