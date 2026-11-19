#!/usr/bin/env python3
"""Compare every Rust workload with its committed JSON, without writing reports."""

import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]


def differences(expected, actual, path="$"):
    """Return field paths whose types or values differ, including missing fields."""
    if type(expected) is not type(actual):
        return [path]
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            return [path]
        return [p for key in expected for p in differences(expected[key], actual[key], f"{path}.{key}")]
    if isinstance(expected, list):
        if len(expected) != len(actual):
            return [path]
        return [p for i, (a, b) in enumerate(zip(expected, actual)) for p in differences(a, b, f"{path}[{i}]")]
    if isinstance(expected, float):
        return [] if expected.hex() == actual.hex() else [path]
    return [] if expected == actual else [path]


def main():
    spec = json.loads((ROOT / "verification/workloads/root_minuit2_v6_36_08.json").read_text())
    passed = 0
    for workload in spec["workloads"]:
        name = workload["id"]
        committed = subprocess.check_output(
            ["git", "show", f"HEAD:reports/verification/raw/rust/{name}.json"], cwd=ROOT, text=True
        )
        output = subprocess.check_output(
            ["cargo", "run", "--release", "--bin", "ref_compare_runner", "--", "--workload", name],
            cwd=ROOT, text=True
        )
        changed = differences(json.loads(committed), json.loads(output))
        if changed:
            print(f"FAIL {name}: {', '.join(changed)}")
        else:
            passed += 1
    print(f"{passed} of {len(spec['workloads'])} equal")
    return int(passed != len(spec["workloads"]))


if __name__ == "__main__":
    raise SystemExit(main())
