# Differential Verification Summary

Reference repo: `https://github.com/root-project/root`
Reference subtree: `math/minuit2`
Reference tag: `v6-36-08`
Reference commit: `a8ca1b23e38d7dbe0ff24027894ca0f2ad65f1bd`

## Status Counts

- pass: **12**
- warn: **0**
- fail: **0**

## Per-Workload Results

| Workload | Status | Issues | Warnings |
|---|---|---|---|
| `quadratic3_fixx_migrad` | `pass` | - | - |
| `quadratic3_fixx_hesse` | `pass` | - | waived nfcn relative diff 0.667 > 0.5: D-001: inherits the fixed-parameter Migrad seed NFCN accounting divergence; ROOT v6-36-08 counts the seed central-difference probes before Hesse, while minuit2-rs covariance remains within tolerance. |
| `rosenbrock2_migrad` | `pass` | - | - |
| `quadratic2_minos_p0` | `pass` | - | - |
| `quadratic2_simplex` | `pass` | - | - |
| `rosenbrock2_minimize` | `pass` | - | - |
| `quadratic2_minos_p1` | `pass` | - | - |
| `quadratic2_lower_limited_migrad` | `pass` | - | - |
| `rosenbrock2_migrad_strategy2` | `pass` | - | - |
| `quadratic2_scan_p0` | `pass` | - | - |
| `quadratic2_contours_01` | `pass` | - | - |
| `quadratic2_no_g2_migrad` | `pass` | - | - |

## Artifacts

- `reports/verification/diff_results.csv`
- `reports/verification/raw/ref/*.json`
- `reports/verification/raw/rust/*.json`

## Notes

- `fail` means correctness metrics exceeded workload tolerances.
- `warn` means correctness metrics passed, but NFCN divergence exceeded warning threshold.
