#!/usr/bin/env python3
"""Measure binary64 distances against the frozen workload tolerances."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import struct
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import cast

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "verification/workloads/root_minuit2_v6_36_08.json"
SIGN = 1 << 63
MASK = (1 << 64) - 1
type Number = int | float | bool


def ulp_distance(reference: float, candidate: float) -> int | None:
    """Complement negative binary64 bits and set the sign bit on nonnegative bits to order their integers.

    Their distance counts representable steps, keeping signed zeros adjacent; equal nonfinite classes return zero, and other nonfinite pairs return None.
    """
    if math.isnan(reference) and math.isnan(candidate):
        return 0
    if not math.isfinite(reference) or not math.isfinite(candidate):
        return 0 if reference == candidate else None

    def ordered(value: float) -> int:
        bits = int.from_bytes(struct.pack(">d", value), "big")
        return (~bits & MASK) if bits & SIGN else bits | SIGN

    return abs(ordered(reference) - ordered(candidate))


def read_object(path: Path) -> dict[str, object]:
    value = cast("object", json.loads(path.read_text()))
    return object_dict(value)


def object_dict(value: object) -> dict[str, object]:
    if not isinstance(value, dict):
        raise TypeError("Expected a JSON object")
    result = cast("dict[object, object]", value)
    if not all(isinstance(key, str) for key in result):
        raise ValueError("Expected string keys")
    return cast("dict[str, object]", result)


def flatten(value: object) -> tuple[dict[str, Number], dict[str, str]]:
    """Collect numeric leaves and node kinds, including nulls and empty containers."""
    numbers: dict[str, Number] = {}
    nodes: dict[str, str] = {}

    def visit(item: object, path: str) -> None:
        if isinstance(item, dict):
            nodes[path] = "object"
            for key, child in sorted(object_dict(cast("object", item)).items()):
                visit(child, f"{path}.{key}" if path else key)
        elif isinstance(item, list):
            children = cast("list[object]", item)
            nodes[path] = f"array[{len(children)}]"
            for index, child in enumerate(children):
                visit(child, f"{path}[{index}]")
        elif isinstance(item, (int, float)):
            numbers[path] = item
            nodes[path] = "boolean" if isinstance(item, bool) else "number"
        elif item is None:
            nodes[path] = "null"
        elif isinstance(item, str):
            nodes[path] = "string"
        else:
            raise ValueError(f"Unsupported JSON value at {path}")

    visit(value, "")
    return numbers, nodes


def tolerance_for(
    path: str, workload: str, tolerances: dict[str, object], manifest_path: str
) -> tuple[dict[str, object] | None, str]:
    base = path.split("[", 1)[0]
    key = {
        "fval": "fval_abs",
        "edm": "edm_abs",
        "params": "param_abs",
        "errors": "error_abs",
        "covariance": "cov_abs",
        "minos.lower": "minos_abs",
        "minos.upper": "minos_abs",
        "nfcn": "nfcn_rel_warn",
    }.get(base)
    value = tolerances.get(key) if key is not None else None
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return None, "no tolerance in the manifest; distance reported only"
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"Invalid tolerance for {workload}.{path}")
    reason = f"tolerance of workload {workload} in {manifest_path}"
    result: dict[str, object] = {
        "atol": float(value) if key != "nfcn_rel_warn" else 0.0,
        "rtol": float(value) if key == "nfcn_rel_warn" else 0.0,
        "precision": "integer" if key == "nfcn_rel_warn" else "f64",
        "signed_zero": "strict",
        "reason": reason,
    }
    if key == "nfcn_rel_warn" and "nfcn_rel_waiver" in tolerances:
        result["waiver"] = tolerances["nfcn_rel_waiver"]
    return result, reason


@dataclass
class FieldResult:
    """Count paired numeric leaves; unmatched leaves and containers carry structure records."""

    name: str
    workload: str
    path: str
    tolerance: dict[str, object] | None
    reason: str
    count: int = 0
    max_abs_err: float | int | None = 0.0
    max_rel_err: float | None = 0.0
    max_ulp: int | None = 0
    worst_abs_index: str | None = None
    nonfinite_mismatch: int = 0
    sign_mismatch: int = 0
    subnormal_in_ref: int = 0
    int_mismatch: int = 0
    structure: list[dict[str, str]] = field(default_factory=list)
    first_failures: list[dict[str, object]] = field(default_factory=list)

    def as_json(self) -> dict[str, object]:
        return {
            "name": self.name,
            "workload": self.workload,
            "path": self.path,
            "pass": not self.first_failures,
            "count": self.count,
            "failures": len(self.first_failures),
            "tolerance": self.tolerance,
            "reason": self.reason,
            "max_abs_err": self.max_abs_err,
            "max_rel_err": self.max_rel_err,
            "max_ulp": self.max_ulp,
            "worst_abs_index": self.worst_abs_index,
            "nonfinite_mismatch": self.nonfinite_mismatch,
            "sign_mismatch": self.sign_mismatch,
            "subnormal_in_ref": self.subnormal_in_ref,
            "int_mismatch": self.int_mismatch,
            "structure": self.structure,
            "first_failures": self.first_failures,
        }


def display_number(value: Number) -> Number | str:
    """Keep report JSON valid when input numbers are nonfinite."""
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def compare_number(result: FieldResult, reference: Number, candidate: Number) -> None:
    """Require exact integers, booleans, zero signs, and nonfinite classes independently of numeric tolerances."""
    result.count = 1
    a, b = float(reference), float(candidate)
    result.max_ulp = (
        0
        if isinstance(reference, int) and isinstance(candidate, int)
        else ulp_distance(a, b)
    )
    result.subnormal_in_ref = int(0 < abs(a) < sys.float_info.min)
    result.sign_mismatch = int(math.copysign(1.0, a) != math.copysign(1.0, b))
    exact = isinstance(reference, int) or isinstance(candidate, int)
    boolean_mismatch = isinstance(reference, bool) != isinstance(candidate, bool)
    result.int_mismatch = int(exact and (reference != candidate or boolean_mismatch))
    failed = bool(result.int_mismatch)
    if not math.isfinite(a) or not math.isfinite(b):
        result.nonfinite_mismatch = int(result.max_ulp is None)
        result.max_abs_err = None if result.nonfinite_mismatch else 0.0
        result.max_rel_err = None if result.nonfinite_mismatch else 0.0
        failed |= bool(result.nonfinite_mismatch)
    else:
        # Subtract integers before converting, so large counters retain exact errors.
        error = abs(reference - candidate)
        scale = max(abs(reference), abs(candidate))
        relative = error / scale if scale else 0.0
        if not math.isfinite(error):
            relative = abs(a / scale - b / scale)
        result.max_abs_err = error if math.isfinite(error) else None
        result.max_rel_err = relative if math.isfinite(relative) else None
        if error:
            result.worst_abs_index = result.path
        failed |= a == b == 0.0 and bool(result.sign_mismatch)
        if not exact and result.tolerance is not None:
            atol = cast("float", result.tolerance["atol"])
            failed |= error > atol
    if failed:
        reason = (
            "integer or boolean mismatch"
            if result.int_mismatch
            else "nonfinite mismatch"
            if result.nonfinite_mismatch
            else "signed zero mismatch"
            if a == b == 0.0 and result.sign_mismatch
            else "absolute tolerance exceeded"
        )
        result.first_failures.append(
            {
                "index": result.path,
                "reference": display_number(reference),
                "candidate": display_number(candidate),
                "reason": reason,
            }
        )


def compare_workload(
    workload: str,
    reference: object,
    candidate: object,
    tolerances: dict[str, object],
    manifest_path: str = "verification/workloads/root_minuit2_v6_36_08.json",
) -> list[FieldResult]:
    left, left_nodes = flatten(reference)
    right, right_nodes = flatten(candidate)
    changed = {
        path
        for path in left_nodes.keys() | right_nodes.keys()
        if left_nodes.get(path) != right_nodes.get(path)
    }
    results: list[FieldResult] = []
    for path in sorted(left.keys() | right.keys() | changed):
        tolerance, reason = tolerance_for(path, workload, tolerances, manifest_path)
        result = FieldResult(f"{workload}.{path}", workload, path, tolerance, reason)
        if path in changed:
            mismatch = {
                "path": path,
                "reference": left_nodes.get(path, "missing"),
                "candidate": right_nodes.get(path, "missing"),
            }
            result.structure.append(mismatch)
            result.first_failures.append(
                {"index": path, "reason": "structure mismatch"}
            )
        if path in left and path in right:
            compare_number(result, left[path], right[path])
        results.append(result)
    return results


def fingerprint(path: Path) -> dict[str, str]:
    return {
        "path": path.resolve().relative_to(ROOT).as_posix(),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def summary(results: list[FieldResult]) -> dict[str, object]:
    return {
        "fields": len(results),
        "fields_failed": sum(bool(result.first_failures) for result in results),
        "values": sum(result.count for result in results),
        "failures": sum(len(result.first_failures) for result in results),
        "max_rel_err": (
            None
            if any(result.max_rel_err is None for result in results)
            else max((result.max_rel_err or 0.0 for result in results), default=0.0)
        ),
        "max_ulp": (
            None
            if any(result.max_ulp is None for result in results)
            else max((result.max_ulp or 0 for result in results), default=0)
        ),
    }


def build_report(ref_dir: Path, rust_dir: Path, manifest: Path) -> dict[str, object]:
    spec = read_object(manifest)
    entries = spec.get("workloads")
    if not isinstance(entries, list) or not entries:
        raise ValueError("Expected a nonempty workload list")
    reference_inputs: list[dict[str, str]] = []
    candidate_inputs: list[dict[str, str]] = []
    workloads: list[dict[str, object]] = []
    results: list[FieldResult] = []
    seen: set[str] = set()
    manifest_path = manifest.resolve().relative_to(ROOT).as_posix()
    for entry in cast("list[object]", entries):
        workload = object_dict(entry)
        name = workload.get("id")
        if not isinstance(name, str) or not name or Path(name).name != name:
            raise ValueError("Expected a workload id without path components")
        if name in seen:
            raise ValueError(f"Duplicate workload: {name}")
        seen.add(name)
        ref_path, rust_path = ref_dir / f"{name}.json", rust_dir / f"{name}.json"
        reference_inputs.append(fingerprint(ref_path))
        candidate_inputs.append(fingerprint(rust_path))
        fields = compare_workload(
            name,
            read_object(ref_path),
            read_object(rust_path),
            object_dict(workload.get("tolerances", {})),
            manifest_path,
        )
        results.extend(fields)
        workloads.append(
            {
                "id": name,
                "pass": not any(result.first_failures for result in fields),
                "summary": summary(fields),
                "fields": [result.name for result in fields],
            }
        )
    return {
        "tool": "parity",
        "version": 1,
        "pass": not any(result.first_failures for result in results),
        "inputs": {
            "reference": reference_inputs,
            "candidate": candidate_inputs,
            "config": fingerprint(manifest),
        },
        "summary": summary(results),
        "workloads": workloads,
        "fields": [result.as_json() for result in results],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    _ = parser.add_argument("--ref-dir", type=Path, required=True)
    _ = parser.add_argument("--rust-dir", type=Path, required=True)
    _ = parser.add_argument("--workloads", type=Path, default=MANIFEST)
    _ = parser.add_argument(
        "--report", type=Path, default=ROOT / "parity/parity_ulp.json"
    )
    args = parser.parse_args(argv)
    try:
        report = build_report(
            cast("Path", args.ref_dir),
            cast("Path", args.rust_dir),
            cast("Path", args.workloads),
        )
        report_path = cast("Path", args.report)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        _ = report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    except (OSError, TypeError, ValueError, OverflowError) as error:
        print(f"parity: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report["summary"]))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
