#!/usr/bin/env python3
"""
Generate executed-surface coverage for ROOT Minuit2 reference runner.

Builds an instrumented reference runner, executes all configured workloads,
and exports function-level coverage for ROOT `math/minuit2`.

Outputs:
  - reports/verification/reference_coverage/executed_functions.csv
  - reports/verification/reference_coverage/unexecuted_functions.csv
  - reports/verification/reference_coverage/summary.md
  - reports/verification/reference_coverage/manifest.json
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKLOADS = REPO_ROOT / "verification" / "workloads" / "root_minuit2_v6_36_08.json"
DEFAULT_REPORT_DIR = REPO_ROOT / "reports" / "verification" / "reference_coverage"
DEFAULT_ROOT_TAG = "v6-36-08"


def run(cmd: list[str], *, env: dict[str, str] | None = None) -> str:
    output = subprocess.check_output(cmd, cwd=REPO_ROOT, env=env, text=True)
    return output


def which_llvm_tool(tool: str) -> list[str]:
    env_override = os.environ.get(tool.upper().replace("-", "_"))
    if env_override:
        return env_override.split()
    direct = shutil.which(tool)
    if direct:
        return [direct]
    xcrun = shutil.which("xcrun")
    if xcrun:
        resolved = subprocess.check_output([xcrun, "--find", tool], text=True).strip()
        if resolved:
            return [resolved]
    raise FileNotFoundError(f"could not resolve tool: {tool}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate C++ reference coverage for ROOT Minuit2")
    parser.add_argument("--root-tag", default=DEFAULT_ROOT_TAG, help="ROOT tag to check out")
    parser.add_argument("--workloads", default=str(DEFAULT_WORKLOADS), help="Workload JSON")
    parser.add_argument("--report-dir", default=str(DEFAULT_REPORT_DIR), help="Output report directory")
    return parser.parse_args()


def build_instrumented_runner(root_tag: str) -> Path:
    env = os.environ.copy()
    cov_flags = "-fprofile-instr-generate -fcoverage-mapping -O0 -g"
    # Let CMake pick these up from the environment.
    env["CFLAGS"] = cov_flags
    env["CXXFLAGS"] = cov_flags
    env["LDFLAGS"] = "-fprofile-instr-generate"
    env["ROOT_BUILD_DIR"] = "third_party/root_ref_build/minuit2_cov"
    env["RUNNER_BUILD_DIR"] = "third_party/root_ref_build/ref_runner_cov"
    env["ROOT_CMAKE_EXTRA_ARGS"] = "-DCMAKE_BUILD_TYPE=Debug"
    env["RUNNER_CMAKE_EXTRA_ARGS"] = "-DCMAKE_BUILD_TYPE=Debug"

    out = run(["bash", "scripts/build_root_reference_runner.sh", root_tag], env=env)
    lines = [line.strip() for line in out.splitlines() if line.strip()]
    if not lines:
        raise RuntimeError("build script produced no output")
    runner = Path(lines[-1])
    if not runner.exists():
        raise FileNotFoundError(f"instrumented runner missing: {runner}")
    return runner


def load_workloads(path: Path) -> list[str]:
    data = json.loads(path.read_text())
    return [w["id"] for w in data["workloads"]]


def normalize_file(path: str) -> str:
    p = Path(path)
    try:
        return str(p.resolve().relative_to(REPO_ROOT.resolve()))
    except Exception:
        return str(p)


def in_scope(filename: str) -> bool:
    return "/math/minuit2/" in filename or filename.startswith("math/minuit2/")


def main() -> int:
    args = parse_args()
    workloads_path = Path(args.workloads)
    report_dir = Path(args.report_dir)
    raw_dir = report_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    runner = build_instrumented_runner(args.root_tag)
    workload_ids = load_workloads(workloads_path)

    profraw_dir = raw_dir / "profraw"
    if profraw_dir.exists():
        for f in profraw_dir.glob("*.profraw"):
            f.unlink()
    profraw_dir.mkdir(parents=True, exist_ok=True)

    for wid in workload_ids:
        env = os.environ.copy()
        env["LLVM_PROFILE_FILE"] = str(profraw_dir / f"ref_runner_{wid}_%p.profraw")
        out = subprocess.check_output([str(runner), "--workload", wid], cwd=REPO_ROOT, env=env, text=True)
        (raw_dir / f"{wid}.json").write_text(out.strip() + "\n")

    profraw_files = sorted(profraw_dir.glob("*.profraw"))
    if not profraw_files:
        raise RuntimeError("no .profraw files produced")

    llvm_profdata = which_llvm_tool("llvm-profdata")
    llvm_cov = which_llvm_tool("llvm-cov")

    merged = raw_dir / "merged.profdata"
    subprocess.check_call([*llvm_profdata, "merge", "-sparse", *[str(p) for p in profraw_files], "-o", str(merged)])

    export_json_path = raw_dir / "llvm_cov_export.json"
    export_json = subprocess.check_output(
        [*llvm_cov, "export", str(runner), f"-instr-profile={merged}"], cwd=REPO_ROOT, text=True
    )
    export_json_path.write_text(export_json)
    payload: dict[str, Any] = json.loads(export_json)

    functions: list[tuple[str, str, int]] = []
    total_in_scope = 0
    executed_in_scope = 0

    file_totals: dict[str, dict[str, int]] = {}
    for data_entry in payload.get("data", []):
        for f in data_entry.get("files", []):
            filename = str(f.get("filename", ""))
            if not in_scope(filename):
                continue
            norm = normalize_file(filename)
            summary = f.get("summary", {}).get("lines", {})
            covered = int(summary.get("covered", 0))
            count = int(summary.get("count", 0))
            file_totals[norm] = {"covered": covered, "count": count}

        for fn in data_entry.get("functions", []):
            name = str(fn.get("name", ""))
            count = int(fn.get("count", 0))
            files = [str(f) for f in fn.get("filenames", [])]
            if not files:
                continue
            primary = files[0]
            if not in_scope(primary):
                continue
            total_in_scope += 1
            if count > 0:
                executed_in_scope += 1
            functions.append((name, normalize_file(primary), count))

    functions.sort(key=lambda x: (x[2] == 0, x[1], x[0]))

    executed_path = report_dir / "executed_functions.csv"
    unexecuted_path = report_dir / "unexecuted_functions.csv"
    report_dir.mkdir(parents=True, exist_ok=True)

    with executed_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["function", "file", "count"])
        for name, file_name, count in functions:
            if count > 0:
                writer.writerow([name, file_name, count])

    with unexecuted_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["function", "file", "count"])
        for name, file_name, count in functions:
            if count == 0:
                writer.writerow([name, file_name, count])

    file_covered = sum(v["covered"] for v in file_totals.values())
    file_count = sum(v["count"] for v in file_totals.values())
    line_coverage = 0.0 if file_count == 0 else (100.0 * file_covered / file_count)
    function_coverage = 0.0 if total_in_scope == 0 else (100.0 * executed_in_scope / total_in_scope)

    summary_md = report_dir / "summary.md"
    summary_md.write_text(
        "\n".join(
            [
                "# Reference Coverage Summary (ROOT Minuit2)",
                "",
                f"- Reference tag: `{args.root_tag}`",
                f"- Workloads executed: `{len(workload_ids)}`",
                f"- In-scope functions (math/minuit2): **{total_in_scope}**",
                f"- Executed functions: **{executed_in_scope}**",
                f"- Function coverage (executed/in-scope): **{function_coverage:.2f}%**",
                f"- File line coverage (math/minuit2 files in export): **{line_coverage:.2f}%**",
                "",
                "## Artifacts",
                "",
                "- `reports/verification/reference_coverage/executed_functions.csv`",
                "- `reports/verification/reference_coverage/unexecuted_functions.csv`",
                "- `reports/verification/reference_coverage/raw/llvm_cov_export.json`",
            ]
        )
        + "\n"
    )

    manifest = {
        "reference_tag": args.root_tag,
        "runner_binary": str(runner),
        "workloads": workload_ids,
        "counts": {
            "functions_in_scope": total_in_scope,
            "functions_executed": executed_in_scope,
            "function_coverage_pct": round(function_coverage, 4),
            "line_coverage_pct": round(line_coverage, 4),
        },
        "artifacts": {
            "executed_functions_csv": str(executed_path),
            "unexecuted_functions_csv": str(unexecuted_path),
            "llvm_cov_export_json": str(export_json_path),
            "summary_md": str(summary_md),
        },
    }
    (report_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"Wrote {summary_md}")
    print(
        f"Function coverage: {executed_in_scope}/{total_in_scope} "
        f"({function_coverage:.2f}%), line coverage {line_coverage:.2f}%"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
