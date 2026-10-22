# Verification Scorecard (Claim-Oriented)

- Generated at: `2026-02-11T13:57:34Z`
- Rust commit: `3de1c79776e9c226606ead74fe556f5e06aa3457`
- Reference repo: `https://github.com/root-project/root`
- Reference subtree: `math/minuit2`
- Reference tag: `v6-36-08`
- Reference commit: `a8ca1b23e38d7dbe0ff24027894ca0f2ad65f1bd`

## Evidence Snapshot

- Differential workloads: **12** (pass=10, warn=2, fail=0)
- Traceability symbols: **415** (implemented=303, waived=112, unresolved=0)
- Known ROOT P0 gaps: **0**
- Coverage (line): `--no-default-features` **73.12%**, `--all-features` **69.51%**
- Reference C++ executed-surface: 618/1840 functions (**33.59%**) across `math/minuit2`
- Executed-surface unmapped gaps: P0=0, P1=48, P2=425 (gate=NO)
- Benchmark serial scan (`default`): **2.441 µs**
- Benchmark scan in `parallel` feature run: serial=2.138 µs, parallel=28.89 µs
- Scan speedup (`serial/parallel`, parallel feature run): **0.07x**

## Claim Gates

| Claim | Gate | Status |
|---|---|---|
| Numerical parity on covered workloads | `diff fail == 0` | **YES** |
| Full symbol/function traceability | `traceability unresolved == 0` | **YES** |
| ROOT P0 regression parity completeness | `known P0 gaps == 0` | **YES** |
| Executed-surface mapping completeness | `executed-surface P0/P1 == 0` | **NO** |
| Full 1:1 functional coverage claim | all above gates true | **NO** |
| Full 100% verifiable coverage claim | 1:1 gate + zero warnings | **NO** |

## Blocking Gaps

- Executed-surface mapping gate fails; see `reports/verification/executed_surface_mapping.md`.
- Differential warnings remain (NFCN divergence); see `reports/verification/diff_summary.md`.

## Reproduce

```bash
scripts/build_root_reference_runner.sh v6-36-08
python3 scripts/compare_ref_vs_rust.py
python3 scripts/generate_traceability_matrix.py
python3 scripts/check_traceability_gate.py --mode non-regression
python3 scripts/generate_executed_surface_mapping.py
python3 scripts/check_executed_surface_gate.py --mode non-regression
cargo test --no-default-features
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo test --all-features
cargo llvm-cov --no-default-features --summary-only > reports/coverage/core_coverage_raw.txt
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo llvm-cov --all-features --summary-only > reports/coverage/all_features_coverage_raw.txt
python3 scripts/generate_coverage_reports.py
python3 scripts/generate_reference_coverage.py --root-tag v6-36-08
cargo bench --bench benchmarks -- --noplot > reports/benchmarks/default_raw.txt
cargo bench --features parallel --bench benchmarks -- --noplot > reports/benchmarks/parallel_raw.txt
python3 scripts/generate_benchmark_report.py
python3 scripts/generate_verification_scorecard.py
```

