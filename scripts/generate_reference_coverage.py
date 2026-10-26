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
import re
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


BUILD_ROOT_DIR = "third_party/root_ref_build/minuit2_cov"
BUILD_RUNNER_DIR = "third_party/root_ref_build/ref_runner_cov"


def validate_build_paths() -> None:
    """Check every path that the build script deletes or writes, before any write.

    The build script runs in the repo root, so a relative path means a path under it.
    """
    root_ref_dir = os.environ.get("ROOT_REF_DIR", "third_party/root_ref")
    for rel in (
        "third_party",
        root_ref_dir,
        BUILD_ROOT_DIR,
        BUILD_RUNNER_DIR,
        f"{BUILD_RUNNER_DIR}/ref_runner",
    ):
        safe_output(REPO_ROOT / rel)


def build_instrumented_runner(root_tag: str) -> Path:
    validate_build_paths()
    env = os.environ.copy()
    cov_flags = "-fprofile-instr-generate -fcoverage-mapping -O0 -g"
    # Let CMake pick these up from the environment.
    env["CFLAGS"] = cov_flags
    env["CXXFLAGS"] = cov_flags
    env["LDFLAGS"] = "-fprofile-instr-generate"
    env["ROOT_BUILD_DIR"] = BUILD_ROOT_DIR
    env["RUNNER_BUILD_DIR"] = BUILD_RUNNER_DIR
    env["ROOT_CMAKE_EXTRA_ARGS"] = "-DCMAKE_BUILD_TYPE=Debug"
    env["RUNNER_CMAKE_EXTRA_ARGS"] = "-DCMAKE_BUILD_TYPE=Debug"

    out = run(["bash", "scripts/build_root_reference_runner.sh", root_tag], env=env)
    lines = [line.strip() for line in out.splitlines() if line.strip()]
    if not lines:
        raise RuntimeError("build script produced no output")
    runner = safe_output(REPO_ROOT / lines[-1])
    if not runner.exists():
        raise FileNotFoundError(f"instrumented runner missing: {runner}")
    return runner


WORKLOAD_ID_RE = re.compile(r"^[A-Za-z0-9._-]+$")


def load_workloads(path: Path) -> list[str]:
    data = json.loads(path.read_text())
    ids = [str(w["id"]) for w in data["workloads"]]
    for wid in ids:
        if not WORKLOAD_ID_RE.fullmatch(wid) or wid in (".", ".."):
            raise SystemExit(
                f"error: workload id {wid!r} in {path} is not allowed. "
                "Use only letters, digits, '.', '_' and '-', and not '.' or '..'."
            )
    return ids


def safe_output(path: Path) -> Path:
    """Resolve `path` (symlinks followed) and require it to be inside the repo.

    Return the resolved path. Use it for every write, so a symlink or a ".." part
    cannot send a report outside the repo. Exit with an error otherwise.
    """
    resolved = Path(path).resolve(strict=False)
    try:
        resolved.relative_to(REPO_ROOT)
    except ValueError:
        raise SystemExit(
            f"error: {path} resolves to {resolved}, which is outside the repository ({REPO_ROOT}). "
            "Pass a path inside the repository and remove symlinks that point outside it."
        ) from None
    return resolved


def repo_relative(path: Path) -> str:
    """Return the repo-relative POSIX form of `path`, so committed reports hold no local paths."""
    return safe_output(path).relative_to(REPO_ROOT).as_posix()


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
    workload_ids = load_workloads(workloads_path)

    # Validate every output and build path before the first write. Keep the resolved paths.
    validate_build_paths()
    report_dir = safe_output(Path(args.report_dir))
    raw_dir = safe_output(report_dir / "raw")
    profraw_dir = safe_output(raw_dir / "profraw")
    workload_outputs = {wid: safe_output(raw_dir / f"{wid}.json") for wid in workload_ids}
    for wid in workload_ids:
        safe_output(profraw_dir / f"ref_runner_{wid}_0.profraw")
    merged = safe_output(raw_dir / "merged.profdata")
    export_json_path = safe_output(raw_dir / "llvm_cov_export.json")
    executed_path = safe_output(report_dir / "executed_functions.csv")
    unexecuted_path = safe_output(report_dir / "unexecuted_functions.csv")
    summary_md = safe_output(report_dir / "summary.md")
    manifest_path = safe_output(report_dir / "manifest.json")

    raw_dir.mkdir(parents=True, exist_ok=True)

    runner = build_instrumented_runner(args.root_tag)
    safe_output(runner)

    if profraw_dir.exists():
        for f in profraw_dir.glob("*.profraw"):
            f.unlink()
    profraw_dir.mkdir(parents=True, exist_ok=True)

    for wid in workload_ids:
        env = os.environ.copy()
        env["LLVM_PROFILE_FILE"] = str(profraw_dir / f"ref_runner_{wid}_%p.profraw")
        out = subprocess.check_output([str(runner), "--workload", wid], cwd=REPO_ROOT, env=env, text=True)
        workload_outputs[wid].write_text(out.strip() + "\n")

    profraw_files = sorted(profraw_dir.glob("*.profraw"))
    if not profraw_files:
        raise RuntimeError("no .profraw files produced")

    llvm_profdata = which_llvm_tool("llvm-profdata")
    llvm_cov = which_llvm_tool("llvm-cov")

    subprocess.check_call([*llvm_profdata, "merge", "-sparse", *[str(p) for p in profraw_files], "-o", str(merged)])

    export_json = subprocess.check_output(
        [*llvm_cov, "export", str(runner), f"-instr-profile={merged}"], cwd=REPO_ROOT, text=True
    )
    payload: dict[str, Any] = json.loads(export_json)
    # Store repo-relative file names in the committed export.
    export_json_path.write_text(export_json.replace(f"{REPO_ROOT.resolve()}/", ""))

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
        "runner_binary": repo_relative(runner),
        "workloads": workload_ids,
        "counts": {
            "functions_in_scope": total_in_scope,
            "functions_executed": executed_in_scope,
            "function_coverage_pct": round(function_coverage, 4),
            "line_coverage_pct": round(line_coverage, 4),
        },
        "artifacts": {
            "executed_functions_csv": repo_relative(executed_path),
            "unexecuted_functions_csv": repo_relative(unexecuted_path),
            "llvm_cov_export_json": repo_relative(export_json_path),
            "summary_md": repo_relative(summary_md),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"Wrote {summary_md}")
    print(
        f"Function coverage: {executed_in_scope}/{total_in_scope} "
        f"({function_coverage:.2f}%), line coverage {line_coverage:.2f}%"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
