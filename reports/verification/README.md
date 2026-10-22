# Cross-Implementation Verification (ROOT vs Rust)

Canonical reference target:
- Repo: `https://github.com/root-project/root`
- Subtree: `math/minuit2`
- Tag: `v6-36-08`
- Commit: `a8ca1b23e38d7dbe0ff24027894ca0f2ad65f1bd`

This directory is a correctness-verification surface. It uses ROOT Minuit2 as a
numerical reference/oracle for behavior comparisons. It does not claim that the
Rust implementation is source-derived from ROOT, and it is not a legal
provenance certificate.

## Reference data regeneration

The `raw/ref/*.json` files are produced by the ROOT C++ reference runner against
the pinned tag/commit above. They are regenerated whenever the reference runner
is rebuilt; rebuilding the runner on a different toolchain or CPU changes only
last-digit floating-point rounding (compiler FMA contraction, `libm`), not the
algorithm or the ROOT version.

Such a rebuild typically shifts the stored values by well under `1e-8` absolute
while `nfcn`, `valid`, and the convergence structure stay identical. The
differential-gate tolerances absorb deltas of this size, so a pure-rounding
refresh is not numerical drift, a regression, or a version change. To avoid
committing rounding-only churn, regenerate these files only alongside a genuine
reference update (a tag/commit bump or an intended behavior change). Regenerate
with `scripts/run_full_verification.sh v6-36-08` (see Reproduce below).

## What is implemented

- ROOT C++ reference runner: `tools/ref_runner_cpp/main.cpp`
- Rust runner: `src/bin/ref_compare_runner.rs`
- Workload set + tolerances: `verification/workloads/root_minuit2_v6_36_08.json`
- Comparator: `scripts/compare_ref_vs_rust.py`
- ROOT build/bootstrap script: `scripts/build_root_reference_runner.sh`
- Traceability matrix generator: `scripts/generate_traceability_matrix.py`
- Traceability gate checker: `scripts/check_traceability_gate.py`
- Executed-surface mapping generator: `scripts/generate_executed_surface_mapping.py`
- Executed-surface gate checker: `scripts/check_executed_surface_gate.py`
- Claim scorecard generator: `scripts/generate_verification_scorecard.py`

## Reproduce

```bash
scripts/run_full_verification.sh v6-36-08

# Or run components manually:
scripts/build_root_reference_runner.sh v6-36-08
python3 scripts/compare_ref_vs_rust.py
python3 scripts/generate_traceability_matrix.py
python3 scripts/check_traceability_gate.py --mode non-regression
python3 scripts/generate_reference_coverage.py --root-tag v6-36-08
python3 scripts/generate_executed_surface_mapping.py
python3 scripts/check_executed_surface_gate.py --mode non-regression
python3 scripts/generate_verification_scorecard.py
```

## Artifacts

- `reports/verification/diff_summary.md`
- `reports/verification/diff_results.csv`
- `reports/verification/raw/ref/*.json`
- `reports/verification/raw/rust/*.json`
- `reports/verification/manifest.json`
- `reports/verification/scorecard.md`
- `reports/verification/known_differences.md`
- `reports/verification/traceability_matrix.csv`
- `reports/verification/traceability_summary.md`
- `reports/verification/reference_coverage/summary.md`
- `reports/verification/executed_surface_mapping.md`
- `reports/verification/executed_surface_gaps.csv`
- `reports/verification/executed_surface_manifest.json`
- `verification/traceability/traceability_baseline.csv`
- `verification/traceability/executed_surface_gaps_baseline.csv`
- `verification/traceability/waivers.csv`

## Last generated snapshot

From `reports/verification/manifest.json`, generated at
`2026-02-11T13:57:34Z`:
- pass: 10
- warn: 2
- fail: 0

Workload count: 12.

Correctness gates passed for all workloads in that snapshot; warnings are NFCN-efficiency deltas only.

This is a real, reproducible differential gate. That snapshot demonstrates correctness parity on the covered workloads, with remaining call-count efficiency differences.

Traceability status snapshot:
- implemented: 303
- waived: 112
- unresolved: 0

Claim snapshot:
- numerical parity on covered workloads: **YES**
- executed-surface mapping completeness (P0/P1): **NO**
- full 1:1 functional coverage claim: **NO**
- full 100% verifiable coverage claim: **NO**
