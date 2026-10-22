#!/usr/bin/env python3
"""
Generate benchmark baseline markdown from Criterion raw logs.

Inputs:
  - reports/benchmarks/default_raw.txt
  - reports/benchmarks/parallel_raw.txt

Output:
  - reports/benchmarks/benchmark_baseline.md
"""

from __future__ import annotations

from pathlib import Path
import re


BLOCK_RE = re.compile(
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


def to_us(value: float, unit: str) -> float:
    return value * UNIT_TO_US[unit]


def parse_benchmarks(raw: str) -> dict[str, tuple[float, str]]:
    out: dict[str, tuple[float, str]] = {}
    for m in BLOCK_RE.finditer(raw):
        name = m.group("name").strip()
        mid = float(m.group("mid"))
        mid_unit = m.group("mid_unit")
        out[name] = (to_us(mid, mid_unit), f"{mid:.4g} {mid_unit}")
    return out


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    bench_dir = root / "reports" / "benchmarks"
    bench_dir.mkdir(parents=True, exist_ok=True)

    default_raw = (bench_dir / "default_raw.txt").read_text()
    parallel_raw = (bench_dir / "parallel_raw.txt").read_text()

    default = parse_benchmarks(default_raw)
    parallel = parse_benchmarks(parallel_raw)

    all_names = sorted(set(default) | set(parallel))

    lines = []
    lines.append("# Benchmark Baseline")
    lines.append("")
    lines.append("- Command (default): `cargo bench --bench benchmarks -- --noplot`")
    lines.append("- Command (parallel): `cargo bench --features parallel --bench benchmarks -- --noplot`")
    lines.append("")
    lines.append("## Median Time By Benchmark")
    lines.append("")
    lines.append("| Benchmark | Default Median | Parallel Median | Parallel vs Default |")
    lines.append("|---|---:|---:|---:|")
    for name in all_names:
        d = default.get(name)
        p = parallel.get(name)
        d_txt = d[1] if d else "-"
        p_txt = p[1] if p else "-"
        if d and p:
            ratio = p[0] / d[0] if d[0] > 0.0 else float("inf")
            ratio_txt = f"{ratio:.2f}x"
        else:
            ratio_txt = "-"
        lines.append(f"| {name} | {d_txt} | {p_txt} | {ratio_txt} |")

    serial_name = "Quadratic 2D: MnScan serial (101 points)"
    parallel_name = "Quadratic 2D: MnScan parallel (101 points)"
    if serial_name in parallel and parallel_name in parallel:
        serial_us = parallel[serial_name][0]
        parallel_us = parallel[parallel_name][0]
        speedup = serial_us / parallel_us if parallel_us > 0.0 else 0.0
        lines.append("")
        lines.append("## Scan Comparison (`parallel` run)")
        lines.append("")
        lines.append(f"- Serial median: **{parallel[serial_name][1]}**")
        lines.append(f"- Parallel median: **{parallel[parallel_name][1]}**")
        lines.append(f"- Speedup (`serial / parallel`): **{speedup:.2f}x**")

    (bench_dir / "benchmark_baseline.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
