#!/usr/bin/env python3
"""
Create CERN dimuon fit figures from examples/cern_dimuon/output/cern_*_curve.csv.
"""

from pathlib import Path
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_curve(path: Path):
    x, y, yerr, m, r = [], [], [], [], []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            x.append(float(row["bin_center"]))
            y.append(float(row["count"]))
            yerr.append(float(row["sigma"]))
            m.append(float(row["model"]))
            r.append(float(row["residual"]))
    return x, y, yerr, m, r


def draw_panel(title: str, x, y, yerr, model, residual, out_path: Path) -> None:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    ax1.errorbar(x, y, yerr=yerr, fmt="o", ms=3, lw=0.8, alpha=0.8, label="Data")
    ax1.plot(x, model, "-", lw=1.8, color="tab:red", label="Model")
    ax1.set_ylabel("Counts / bin")
    ax1.set_title(title)
    ax1.grid(alpha=0.3)
    ax1.legend(loc="best")

    ax2.plot(x, residual, "o", ms=3, color="tab:purple")
    ax2.axhline(0.0, color="k", lw=1.0)
    ax2.set_xlabel("Invariant Mass (GeV)")
    ax2.set_ylabel("Residual")
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    out = repo / "examples" / "cern_dimuon" / "output"
    fig_dir = Path(__file__).resolve().parent / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    x1, y1, yerr1, m1, r1 = read_curve(out / "cern_murun2010b0_jpsi_curve.csv")
    draw_panel(
        "CERN MuRun2010B_0: J/psi Region Fit",
        x1,
        y1,
        yerr1,
        m1,
        r1,
        fig_dir / "cern_murun_jpsi_fit.png",
    )

    x2, y2, yerr2, m2, r2 = read_curve(out / "cern_zmumu_zpeak_curve.csv")
    draw_panel(
        "CERN Zmumu: Z-Peak Fit",
        x2,
        y2,
        yerr2,
        m2,
        r2,
        fig_dir / "cern_zmumu_z_fit.png",
    )


if __name__ == "__main__":
    main()
