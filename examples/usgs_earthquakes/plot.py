#!/usr/bin/env python3
"""
Create USGS Gutenberg-Richter figures from examples/usgs_earthquakes/output/usgs_gutenberg_richter_curve.csv.
"""

from pathlib import Path
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    csv_path = repo / "examples" / "usgs_earthquakes" / "output" / "usgs_gutenberg_richter_curve.csv"
    fig_dir = Path(__file__).resolve().parent / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    m, n, logn, pred, res = [], [], [], [], []
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            m.append(float(row["magnitude_threshold"]))
            n.append(float(row["cumulative_count"]))
            logn.append(float(row["log10_count"]))
            pred.append(float(row["pred_log10"]))
            res.append(float(row["residual"]))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    ax1.plot(m, logn, "o", label="Observed log10 N", ms=4)
    ax1.plot(m, pred, "-", label="Fitted", lw=1.8, color="tab:red")
    ax1.set_ylabel("log10 N(M>=m)")
    ax1.set_title("USGS Gutenberg-Richter Fit")
    ax1.grid(alpha=0.3)
    ax1.legend(loc="best")

    ax2.plot(m, res, "o", ms=4, color="tab:purple")
    ax2.axhline(0.0, color="k", lw=1.0)
    ax2.set_xlabel("Magnitude Threshold m")
    ax2.set_ylabel("Residual")
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(fig_dir / "usgs_gr_fit.png", dpi=220)
    plt.close(fig)

    fig2, ax = plt.subplots(figsize=(9, 5))
    ax.semilogy(m, n, "o-", ms=4)
    ax.set_title("USGS Cumulative Counts")
    ax.set_xlabel("Magnitude Threshold m")
    ax.set_ylabel("N(M>=m)")
    ax.grid(alpha=0.3, which="both")
    fig2.tight_layout()
    fig2.savefig(fig_dir / "usgs_gr_counts_semilogy.png", dpi=220)
    plt.close(fig2)


if __name__ == "__main__":
    main()
