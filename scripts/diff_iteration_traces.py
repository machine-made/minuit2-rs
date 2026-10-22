#!/usr/bin/env python3
"""Generate and diff per-iteration Migrad traces for ROOT vs Rust.

ROOT tracing uses ROOT Minuit2's MnTraceObject hook wired into the local C++
reference runner. MnTraceObject exposes MinimumState at iteration boundaries;
ROOT does not expose the accepted line-search lambda through this hook, so ROOT
records lambda as null while Rust records the accepted lambda from its line
search result.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TRACE_DIR = REPO / "reports" / "verification" / "raw" / "trace"
REF_BIN = REPO / "third_party" / "root_ref_build" / "ref_runner" / "ref_runner"
SUPPORTED = {"rosenbrock2_migrad", "quadratic3_fixx_migrad"}
FIELDS = ("nfcn", "fval", "edm", "grad_norm", "dcovar")


def run(cmd: list[str], env: dict[str, str]) -> None:
    subprocess.run(cmd, cwd=REPO, env=env, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def build_ref() -> Path:
    out = subprocess.run(["bash", "scripts/build_root_reference_runner.sh", "v6-36-08"], cwd=REPO, check=True, stdout=subprocess.PIPE, text=True).stdout.strip().splitlines()
    path = Path(out[-1]) if out else REF_BIN
    return path if path.is_absolute() else REPO / path


def ensure_ref() -> Path:
    return REF_BIN if REF_BIN.exists() else build_ref()


def load_jsonl(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open() as fh:
        for line in fh:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def generate(workload: str) -> tuple[Path, Path]:
    TRACE_DIR.mkdir(parents=True, exist_ok=True)
    root_trace = TRACE_DIR / f"{workload}.root.jsonl"
    rust_trace = TRACE_DIR / f"{workload}.rust.jsonl"
    root_trace.unlink(missing_ok=True)
    rust_trace.unlink(missing_ok=True)

    ref = ensure_ref()
    env = os.environ.copy()
    env["MINUIT2_ROOT_TRACE_JSONL"] = str(root_trace)
    run([str(ref), "--workload", workload], env)
    if not root_trace.exists():
        ref = build_ref()
        run([str(ref), "--workload", workload], env)

    env = os.environ.copy()
    env["MINUIT2_RS_TRACE_JSONL"] = str(rust_trace)
    run(["cargo", "run", "--quiet", "--features", "trace", "--bin", "ref_compare_runner", "--", "--workload", workload], env)
    return root_trace, rust_trace


def differs(a: object, b: object, tol: float) -> bool:
    if a is None or b is None:
        return a != b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return abs(float(a) - float(b)) > tol
    return a != b


def main() -> int:
    ap = argparse.ArgumentParser(description="Diff ROOT and Rust Migrad iteration traces")
    ap.add_argument("--workload", required=True, choices=sorted(SUPPORTED))
    ap.add_argument("--tolerance", type=float, default=1e-12)
    args = ap.parse_args()

    root_path, rust_path = generate(args.workload)
    root = load_jsonl(root_path)
    rust = load_jsonl(rust_path)
    print(f"workload: {args.workload}")
    print(f"root trace: {root_path.relative_to(REPO)} ({len(root)} records)")
    print(f"rust trace: {rust_path.relative_to(REPO)} ({len(rust)} records)")

    for i in range(max(len(root), len(rust))):
        if i >= len(root) or i >= len(rust):
            print(f"first divergent iteration: {i} (length mismatch)")
            return 0
        diffs = [f for f in FIELDS if differs(root[i].get(f), rust[i].get(f), args.tolerance)]
        if diffs:
            print(f"first divergent iteration: {i}")
            print("fields: " + ", ".join(diffs))
            print("root: " + json.dumps(root[i], sort_keys=True))
            print("rust: " + json.dumps(rust[i], sort_keys=True))
            return 0
    print("first divergent iteration: none")
    return 0


if __name__ == "__main__":
    sys.exit(main())
