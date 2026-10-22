#!/usr/bin/env python3
"""
Benchmark C++ (ROOT Minuit2) vs Rust scientific demos with in-process solver timing.

Method:
- each benchmark job loads/prepares data once
- warmups are run but not measured
- measured samples time only solver work (Minuit optimize + uncertainty step)
- jobs are interleaved in random order across impl/case/batch
- optional CPU pinning on Linux via taskset
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import platform
import random
import re
import shutil
import statistics
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CASES: list[tuple[str, str]] = [
    ("noaa_co2", "NOAA CO2"),
    ("nist_strd", "NIST StRD"),
    ("usgs_earthquakes", "USGS EQ"),
    ("cern_dimuon", "CERN dimuon"),
]


def quantile(values: list[float], q: float) -> float:
    if not values:
        return float("nan")
    if len(values) == 1:
        return values[0]
    s = sorted(values)
    idx = q * (len(s) - 1)
    lo = int(math.floor(idx))
    hi = int(math.ceil(idx))
    if lo == hi:
        return s[lo]
    w = idx - lo
    return s[lo] * (1.0 - w) + s[hi] * w


def bootstrap_median_ci(samples: list[float], iters: int, alpha: float, rng: random.Random) -> tuple[float, float]:
    if len(samples) < 2:
        v = statistics.median(samples) if samples else float("nan")
        return (v, v)
    meds: list[float] = []
    n = len(samples)
    for _ in range(iters):
        draw = [samples[rng.randrange(n)] for _ in range(n)]
        meds.append(statistics.median(draw))
    low = quantile(meds, alpha / 2.0)
    high = quantile(meds, 1.0 - alpha / 2.0)
    return (low, high)


def permutation_pvalue_median(a: list[float], b: list[float], iters: int, rng: random.Random) -> float:
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    obs = abs(statistics.median(a) - statistics.median(b))
    combined = a + b
    na = len(a)
    count = 0
    for _ in range(iters):
        rng.shuffle(combined)
        d = abs(statistics.median(combined[:na]) - statistics.median(combined[na:]))
        if d >= obs:
            count += 1
    return (count + 1.0) / (iters + 1.0)


def median_abs_deviation(samples: list[float]) -> float:
    if not samples:
        return float("nan")
    m = statistics.median(samples)
    return statistics.median([abs(v - m) for v in samples])


def trimmed(values: list[float], trim_fraction: float) -> list[float]:
    if not values:
        return []
    if trim_fraction <= 0.0:
        return sorted(values)
    s = sorted(values)
    k = int(len(s) * trim_fraction)
    if 2 * k >= len(s):
        return s
    return s[k : len(s) - k]


def trimmed_median(values: list[float], trim_fraction: float) -> float:
    t = trimmed(values, trim_fraction)
    if not t:
        return float("nan")
    return statistics.median(t)


def parse_bench_times(stdout: str, stderr: str, cmd: list[str]) -> list[float]:
    merged = f"{stdout}\n{stderr}"
    for line in merged.splitlines():
        if line.startswith("BENCH_TIMES_S:"):
            payload = line.split(":", 1)[1].strip()
            if not payload:
                return []
            return [float(x) for x in payload.split(",") if x]
    raise RuntimeError(f"benchmark output missing BENCH_TIMES_S line for command: {' '.join(cmd)}")


def run_capture(cmd: list[str]) -> tuple[int, str, str]:
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def first_non_empty_line(text: str) -> str:
    for line in text.splitlines():
        if line.strip():
            return line.strip()
    return ""


def collect_environment(cpu_core: int | None, strict_env: bool) -> tuple[dict[str, object], list[str], list[str], bool]:
    env: dict[str, object] = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "python": platform.python_version(),
    }
    warnings: list[str] = []
    critical: list[str] = []

    system = platform.system()
    taskset_available = shutil.which("taskset") is not None
    affinity_method = "none"
    if cpu_core is not None:
        if system == "Linux" and taskset_available:
            affinity_method = "taskset"
        else:
            msg = f"CPU affinity requested (core={cpu_core}) but not available on {system}."
            warnings.append(msg)
            if strict_env:
                critical.append(msg)
    env["cpu_affinity"] = {
        "requested_core": cpu_core,
        "method": affinity_method,
        "taskset_available": taskset_available,
    }

    if system == "Linux":
        model = ""
        try:
            with open("/proc/cpuinfo", "r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("model name"):
                        model = line.split(":", 1)[1].strip()
                        break
        except OSError:
            pass

        governors: list[str] = []
        gov_paths = sorted(Path("/sys/devices/system/cpu").glob("cpu*/cpufreq/scaling_governor"))
        for gp in gov_paths:
            try:
                governors.append(gp.read_text(encoding="utf-8").strip())
            except OSError:
                continue
        unique_governors = sorted(set(g for g in governors if g))

        no_turbo = None
        no_turbo_path = Path("/sys/devices/system/cpu/intel_pstate/no_turbo")
        if no_turbo_path.exists():
            try:
                no_turbo = no_turbo_path.read_text(encoding="utf-8").strip()
            except OSError:
                pass

        env["linux"] = {
            "cpu_model": model,
            "governors": unique_governors,
            "intel_pstate_no_turbo": no_turbo,
        }

        if unique_governors and any(g != "performance" for g in unique_governors):
            msg = f"Linux CPU governor is not fully 'performance': {unique_governors}"
            warnings.append(msg)
            if strict_env:
                critical.append(msg)
        if not unique_governors:
            warnings.append("Linux governor not readable; cannot verify fixed-frequency policy.")

    elif system == "Darwin":
        def sysctl_value(key: str) -> str:
            rc, out, _ = run_capture(["sysctl", "-n", key])
            return out.strip() if rc == 0 else ""

        brand = sysctl_value("machdep.cpu.brand_string")
        ncpu = sysctl_value("hw.ncpu")
        perflevel = sysctl_value("hw.perflevel0.physicalcpu")

        rc_custom, out_custom, _ = run_capture(["pmset", "-g", "custom"])
        low_power_mode = None
        if rc_custom == 0:
            section = None
            for raw_line in out_custom.splitlines():
                line = raw_line.strip()
                if line.endswith(":"):
                    section = line[:-1]
                    continue
                if section == "AC Power":
                    m = re.search(r"lowpowermode\s+(\d)", line)
                    if m:
                        low_power_mode = int(m.group(1))
                        break
            if low_power_mode is None:
                m = re.search(r"lowpowermode\s+(\d)", out_custom)
                if m:
                    low_power_mode = int(m.group(1))

        rc_batt, out_batt, _ = run_capture(["pmset", "-g", "batt"])
        power_source = first_non_empty_line(out_batt) if rc_batt == 0 else ""

        env["darwin"] = {
            "cpu_brand": brand,
            "ncpu": ncpu,
            "perflevel0_physicalcpu": perflevel,
            "lowpowermode": low_power_mode,
            "power_source": power_source,
        }

        if low_power_mode == 1:
            msg = "macOS Low Power Mode is enabled; this biases performance measurements."
            warnings.append(msg)
            if strict_env:
                critical.append(msg)
        if power_source and "Battery Power" in power_source:
            warnings.append("Running on battery power; plug in for stable frequency behavior.")
    else:
        warnings.append(f"No OS-specific frequency/governor checks implemented for {system}.")

    try:
        la = os.getloadavg()
        env["loadavg"] = {"1min": la[0], "5min": la[1], "15min": la[2]}
    except OSError:
        env["loadavg"] = None

    return env, warnings, critical, taskset_available


def run_bench(cmd: list[str], cwd: Path, cpu_core: int | None, taskset_available: bool) -> list[float]:
    wrapped = cmd
    if cpu_core is not None and platform.system() == "Linux" and taskset_available:
        wrapped = ["taskset", "-c", str(cpu_core), *cmd]

    p = subprocess.run(wrapped, cwd=cwd, check=True, capture_output=True, text=True)
    return parse_bench_times(p.stdout, p.stderr, wrapped)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate C++ vs Rust solver-only timing comparison figure")
    parser.add_argument("--repeats", type=int, default=41, help="measured repeats per benchmark batch")
    parser.add_argument("--warmups", type=int, default=9, help="warmup repeats per benchmark batch")
    parser.add_argument("--batches", type=int, default=9, help="independent benchmark batches per case/impl")
    parser.add_argument("--bootstrap-iters", type=int, default=5000, help="bootstrap iterations for median CI")
    parser.add_argument(
        "--permutation-iters",
        type=int,
        default=10000,
        help="permutation-test iterations for p-value on median difference",
    )
    parser.add_argument("--trim-fraction", type=float, default=0.10, help="fraction trimmed on each tail for robust median")
    parser.add_argument("--cpu-core", type=int, default=None, help="pin jobs to this CPU core (Linux only via taskset)")
    parser.add_argument(
        "--strict-env",
        action="store_true",
        help="fail if environment checks detect critical benchmark-quality issues",
    )
    parser.add_argument("--seed", type=int, default=20260211, help="random seed")
    args = parser.parse_args()

    repeats = max(3, args.repeats)
    warmups = max(0, args.warmups)
    batches = max(1, args.batches)
    trim_fraction = max(0.0, min(0.24, args.trim_fraction))
    bootstrap_iters = max(200, args.bootstrap_iters)
    permutation_iters = max(500, args.permutation_iters)
    rng = random.Random(args.seed)

    repo = Path(__file__).resolve().parents[1]
    out_dir = repo / "examples" / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "comparison_timings.csv"
    samples_path = out_dir / "comparison_samples.csv"
    env_path = out_dir / "comparison_environment.json"
    fig_path = repo / "examples" / "comparison.png"

    env, warnings, critical, taskset_available = collect_environment(args.cpu_core, args.strict_env)
    env["benchmark_config"] = {
        "repeats_per_batch": repeats,
        "warmups_per_batch": warmups,
        "batches": batches,
        "trim_fraction": trim_fraction,
        "bootstrap_iters": bootstrap_iters,
        "permutation_iters": permutation_iters,
        "seed": args.seed,
    }
    env["warnings"] = warnings
    env["critical_issues"] = critical

    if critical:
        env_path.write_text(json.dumps(env, indent=2) + "\n", encoding="utf-8")
        raise RuntimeError(
            "benchmark environment failed strict checks:\n- " + "\n- ".join(critical) + f"\nSee {env_path}"
        )

    rust_bins = {case: repo / "target" / "release" / "examples" / case for case, _ in CASES}
    cpp_bin = repo / "third_party" / "root_ref_build" / "scientific_runner" / "scientific_runner"

    if not cpp_bin.exists():
        raise FileNotFoundError(f"C++ timing binary not found: {cpp_bin}")
    for case, _ in CASES:
        if not rust_bins[case].exists():
            raise FileNotFoundError(f"Rust timing binary not found: {rust_bins[case]}")

    samples: dict[str, dict[str, list[float]]] = {case: {"rust": [], "cpp": []} for case, _ in CASES}
    raw_rows: list[dict[str, str]] = []

    schedule: list[tuple[str, str, int]] = []
    for case, _ in CASES:
        for b in range(batches):
            schedule.append((case, "rust", b))
            schedule.append((case, "cpp", b))
    rng.shuffle(schedule)

    for case, impl, batch_id in schedule:
        if impl == "rust":
            cmd = [
                str(rust_bins[case]),
                "--mode",
                "solve-only",
                "--bench-repeats",
                str(repeats),
                "--bench-warmups",
                str(warmups),
            ]
        else:
            cmd = [
                str(cpp_bin),
                "--case",
                case,
                "--mode",
                "solve-only",
                "--bench-repeats",
                str(repeats),
                "--bench-warmups",
                str(warmups),
            ]

        vals = run_bench(cmd, repo, args.cpu_core, taskset_available)
        if len(vals) != repeats:
            raise RuntimeError(
                f"expected {repeats} benchmark samples, got {len(vals)} for case={case} impl={impl} batch={batch_id}"
            )

        samples[case][impl].extend(vals)
        for i, v in enumerate(vals):
            raw_rows.append(
                {
                    "case_id": case,
                    "impl": impl,
                    "batch_id": str(batch_id),
                    "sample_index": str(i),
                    "time_s": f"{v:.9f}",
                }
            )

    with samples_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["case_id", "impl", "batch_id", "sample_index", "time_s"])
        writer.writeheader()
        writer.writerows(raw_rows)

    rows: list[dict[str, str]] = []
    labels: list[str] = []
    medians_rust: list[float] = []
    medians_cpp: list[float] = []
    err_low_rust: list[float] = []
    err_high_rust: list[float] = []
    err_low_cpp: list[float] = []
    err_high_cpp: list[float] = []

    for case, label in CASES:
        rust_times = samples[case]["rust"]
        cpp_times = samples[case]["cpp"]

        rust_med = statistics.median(rust_times)
        cpp_med = statistics.median(cpp_times)
        rust_ci_lo, rust_ci_hi = bootstrap_median_ci(rust_times, bootstrap_iters, 0.05, rng)
        cpp_ci_lo, cpp_ci_hi = bootstrap_median_ci(cpp_times, bootstrap_iters, 0.05, rng)
        rust_trim_med = trimmed_median(rust_times, trim_fraction)
        cpp_trim_med = trimmed_median(cpp_times, trim_fraction)
        rust_mad = median_abs_deviation(rust_times)
        cpp_mad = median_abs_deviation(cpp_times)
        p_perm = permutation_pvalue_median(rust_times, cpp_times, permutation_iters, rng)

        rows.append(
            {
                "case_id": case,
                "case_label": label,
                "batches": str(batches),
                "repeats_per_batch": str(repeats),
                "warmups_per_batch": str(warmups),
                "samples_per_impl": str(len(rust_times)),
                "trim_fraction": f"{trim_fraction:.4f}",
                "rust_median_s": f"{rust_med:.9f}",
                "rust_trimmed_median_s": f"{rust_trim_med:.9f}",
                "rust_mad_s": f"{rust_mad:.9f}",
                "rust_ci95_low_s": f"{rust_ci_lo:.9f}",
                "rust_ci95_high_s": f"{rust_ci_hi:.9f}",
                "cpp_median_s": f"{cpp_med:.9f}",
                "cpp_trimmed_median_s": f"{cpp_trim_med:.9f}",
                "cpp_mad_s": f"{cpp_mad:.9f}",
                "cpp_ci95_low_s": f"{cpp_ci_lo:.9f}",
                "cpp_ci95_high_s": f"{cpp_ci_hi:.9f}",
                "speedup_cpp_over_rust": f"{(rust_med / cpp_med):.6f}" if cpp_med > 0 else "inf",
                "speedup_cpp_over_rust_trimmed": f"{(rust_trim_med / cpp_trim_med):.6f}" if cpp_trim_med > 0 else "inf",
                "permutation_pvalue_median": f"{p_perm:.6g}",
            }
        )
        labels.append(label)
        medians_rust.append(rust_med)
        medians_cpp.append(cpp_med)
        err_low_rust.append(max(0.0, rust_med - rust_ci_lo))
        err_high_rust.append(max(0.0, rust_ci_hi - rust_med))
        err_low_cpp.append(max(0.0, cpp_med - cpp_ci_lo))
        err_high_cpp.append(max(0.0, cpp_ci_hi - cpp_med))

    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "case_id",
                "case_label",
                "batches",
                "repeats_per_batch",
                "warmups_per_batch",
                "samples_per_impl",
                "trim_fraction",
                "rust_median_s",
                "rust_trimmed_median_s",
                "rust_mad_s",
                "rust_ci95_low_s",
                "rust_ci95_high_s",
                "cpp_median_s",
                "cpp_trimmed_median_s",
                "cpp_mad_s",
                "cpp_ci95_low_s",
                "cpp_ci95_high_s",
                "speedup_cpp_over_rust",
                "speedup_cpp_over_rust_trimmed",
                "permutation_pvalue_median",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    env_path.write_text(json.dumps(env, indent=2) + "\n", encoding="utf-8")

    x = list(range(len(labels)))
    width = 0.38
    fig, ax = plt.subplots(figsize=(11, 6))
    bars_r = ax.bar(
        [v - width / 2 for v in x],
        medians_rust,
        width=width,
        label="Rust",
        color="#1f77b4",
        yerr=[err_low_rust, err_high_rust],
        capsize=4,
    )
    bars_c = ax.bar(
        [v + width / 2 for v in x],
        medians_cpp,
        width=width,
        label="C++ (ROOT Minuit2)",
        color="#d62728",
        yerr=[err_low_cpp, err_high_cpp],
        capsize=4,
    )

    affinity_note = ""
    if args.cpu_core is not None:
        if platform.system() == "Linux" and taskset_available:
            affinity_note = f", pinned core {args.cpu_core}"
        else:
            affinity_note = f", affinity requested (core {args.cpu_core}) not applied"

    ax.set_title(
        "C++ vs Rust Solver Timing (in-process, median, 95% bootstrap CI)\n"
        f"{batches} interleaved batches x {repeats} runs + {warmups} warmups per impl/case{affinity_note}"
    )
    ax.set_ylabel("Solver compute time (seconds)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.grid(axis="y", alpha=0.3)

    max_r = max(m + e for m, e in zip(medians_rust, err_high_rust, strict=False))
    max_c = max(m + e for m, e in zip(medians_cpp, err_high_cpp, strict=False))
    y_top = max(max_r, max_c) * 1.30
    ax.set_ylim(0.0, y_top)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=2)

    for bars in (bars_r, bars_c):
        for b in bars:
            h = b.get_height()
            ax.annotate(
                f"{h:.3f}s",
                xy=(b.get_x() + b.get_width() / 2, h),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    for i, row in enumerate(rows):
        ax.annotate(
            (
                f"p={row['permutation_pvalue_median']}\n"
                f"R/C={row['speedup_cpp_over_rust']}x\n"
                f"MAD R={float(row['rust_mad_s']):.2e}, C={float(row['cpp_mad_s']):.2e}"
            ),
            xy=(x[i], max(medians_rust[i], medians_cpp[i])),
            xytext=(0, 12),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#888", alpha=0.7),
        )

    fig.text(
        0.01,
        0.01,
        f"Robust stats: trim={trim_fraction:.0%} each tail, MAD around median. Environment log: {env_path.name}",
        fontsize=8,
    )

    fig.tight_layout()
    fig.savefig(fig_path, dpi=240)
    plt.close(fig)

    print(f"Wrote {csv_path}")
    print(f"Wrote {samples_path}")
    print(f"Wrote {env_path}")
    print(f"Wrote {fig_path}")

    if warnings:
        print("Benchmark environment warnings:")
        for w in warnings:
            print(f"- {w}")


if __name__ == "__main__":
    main()
