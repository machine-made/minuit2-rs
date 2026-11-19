"""Check the Rosenbrock covariance error with a user gradient in the comparison report.

Usage: python3 -I scripts/check_gradient_covariance.py
"""

import csv
import sys
from pathlib import Path

LIMIT = 0.0765  # twice the iminuit 2.33.0 value, 0.0382
REPORT = Path(__file__).resolve().parents[1] / "reports/compare/matched_accuracy.csv"

with REPORT.open(newline="") as handle:
    rows = [
        r
        for r in csv.DictReader(handle)
        if (r["group"], r["case"]) == ("gradient", "rosenbrock")
    ]
if len(rows) != 1 or not rows[0]["rust_cov_error"]:
    sys.exit("expected one gradient/rosenbrock row with a rust_cov_error value")
error = float(rows[0]["rust_cov_error"])
print(f"rust_cov_error {error}")
if not error <= LIMIT:
    sys.exit(f"rust_cov_error {error} is above {LIMIT}")
