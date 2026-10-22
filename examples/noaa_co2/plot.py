#!/usr/bin/env python3
"""
Create NOAA CO2 fit figures from examples/noaa_co2/output/noaa_co2_curve.csv.
"""

from pathlib import Path
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    csv_path = repo / "examples" / "noaa_co2" / "output" / "noaa_co2_curve.csv"
    fig_dir = Path(__file__).resolve().parent / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    x = []
    y = []
    yfit = []
    r = []
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            x.append(float(row["decimal_date"]))
            y.append(float(row["observed"]))
            yfit.append(float(row["fitted"]))
            r.append(float(row["residual"]))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    ax1.plot(x, y, ".", ms=2, label="Observed", alpha=0.7)
    ax1.plot(x, yfit, "-", lw=1.6, label="Fitted", color="tab:red")
    ax1.set_ylabel("CO2 (ppm)")
    ax1.set_title("NOAA Mauna Loa CO2: Harmonic Trend Fit")
    ax1.grid(alpha=0.3)
    ax1.legend(loc="upper left")

    ax2.plot(x, r, ".", ms=2, color="tab:purple")
    ax2.axhline(0.0, color="k", lw=1.0)
    ax2.set_xlabel("Decimal Year")
    ax2.set_ylabel("Residual (ppm)")
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(fig_dir / "noaa_co2_fit.png", dpi=220)
    plt.close(fig)

    fig2, ax = plt.subplots(figsize=(8, 5))
    ax.hist(r, bins=40, color="tab:blue", alpha=0.8)
    ax.set_title("NOAA CO2 Residual Distribution")
    ax.set_xlabel("Residual (ppm)")
    ax.set_ylabel("Count")
    ax.grid(alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(fig_dir / "noaa_co2_residual_hist.png", dpi=220)
    plt.close(fig2)


if __name__ == "__main__":
    main()
