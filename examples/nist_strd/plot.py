#!/usr/bin/env python3
"""
Create NIST StRD fit figures from examples/nist_strd/output/nist_*.csv.
"""

from pathlib import Path
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_curve(path: Path):
    x, y, yfit, res = [], [], [], []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            x.append(float(row["x"]))
            y.append(float(row["observed"]))
            yfit.append(float(row["fitted"]))
            res.append(float(row["residual"]))
    return x, y, yfit, res


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    out = repo / "examples" / "nist_strd" / "output"
    fig_dir = Path(__file__).resolve().parent / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    datasets = [
        ("Misra1a", out / "nist_misra1a_curve.csv"),
        ("Hahn1", out / "nist_hahn1_curve.csv"),
        ("Rat43", out / "nist_rat43_curve.csv"),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for i, (name, path) in enumerate(datasets):
        x, y, yfit, res = read_curve(path)
        ax_top = axes[0, i]
        ax_bot = axes[1, i]

        ax_top.plot(x, y, "o", ms=3, alpha=0.8, label="Observed")
        ax_top.plot(x, yfit, "-", lw=1.8, color="tab:red", label="Fitted")
        ax_top.set_title(name)
        ax_top.set_xlabel("x")
        ax_top.set_ylabel("y")
        ax_top.grid(alpha=0.3)
        if i == 0:
            ax_top.legend(loc="best")

        ax_bot.plot(x, res, ".", ms=4, color="tab:purple")
        ax_bot.axhline(0.0, color="k", lw=1.0)
        ax_bot.set_xlabel("x")
        ax_bot.set_ylabel("Residual")
        ax_bot.grid(alpha=0.3)

    fig.suptitle("NIST StRD Fits and Residuals")
    fig.tight_layout()
    fig.savefig(fig_dir / "nist_strd_fits.png", dpi=220)
    plt.close(fig)

    summary_path = out / "nist_summary.csv"
    labels, vals = [], []
    with summary_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["metric"] == "param":
                labels.append(f"{row['dataset']}:{row['param']}")
                vals.append(float(row["rel_error"]))

    fig2, ax = plt.subplots(figsize=(13, 5))
    ax.bar(range(len(vals)), vals, color="tab:blue")
    ax.set_yscale("log")
    ax.set_title("NIST Parameter Relative Error vs Certified Values")
    ax.set_ylabel("Relative Error (log scale)")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=75, ha="right")
    ax.grid(axis="y", alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(fig_dir / "nist_strd_rel_error.png", dpi=220)
    plt.close(fig2)


if __name__ == "__main__":
    main()
