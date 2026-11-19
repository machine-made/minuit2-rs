#!/usr/bin/env python3
"""Compare every Rust workload with its committed JSON, without writing reports."""

import json
import subprocess
from pathlib import Path
from typing import cast

ROOT = Path(__file__).resolve().parents[1]


def differences(expected: object, actual: object, path: str = "$") -> list[str]:
    """Return field paths whose types or values differ, including missing fields."""
    if type(expected) is not type(actual):
        return [path]
    if isinstance(expected, dict):
        left = cast("dict[str, object]", expected)
        right = cast("dict[str, object]", actual)
        if left.keys() != right.keys():
            return [path]
        return [
            p
            for key in left
            for p in differences(left[key], right[key], f"{path}.{key}")
        ]
    if isinstance(expected, list):
        first = cast("list[object]", expected)
        second = cast("list[object]", actual)
        if len(first) != len(second):
            return [path]
        return [
            p
            for i, (a, b) in enumerate(zip(first, second))
            for p in differences(a, b, f"{path}[{i}]")
        ]
    if isinstance(expected, float):
        return [] if expected.hex() == cast("float", actual).hex() else [path]
    return [] if expected == actual else [path]


def main() -> int:
    spec = cast(
        "dict[str, list[dict[str, str]]]",
        json.loads(
            (ROOT / "verification/workloads/root_minuit2_v6_36_08.json").read_text()
        ),
    )
    passed = 0
    for workload in spec["workloads"]:
        name = workload["id"]
        committed = subprocess.check_output(
            ["git", "show", f"HEAD:reports/verification/raw/rust/{name}.json"],
            cwd=ROOT,
            text=True,
        )
        output = subprocess.check_output(
            [
                "cargo",
                "run",
                "--release",
                "--bin",
                "ref_compare_runner",
                "--",
                "--workload",
                name,
            ],
            cwd=ROOT,
            text=True,
        )
        changed = differences(
            cast("object", json.loads(committed)), cast("object", json.loads(output))
        )
        if changed:
            print(f"FAIL {name}: {', '.join(changed)}")
        else:
            passed += 1
    print(f"{passed} of {len(spec['workloads'])} equal")
    return int(passed != len(spec["workloads"]))


if __name__ == "__main__":
    raise SystemExit(main())
