#!/usr/bin/env bash
set -euo pipefail

ROOT_TAG="${1:-v6-36-08}"
RUN_TESTS="${RUN_TESTS:-1}"
RUN_COVERAGE="${RUN_COVERAGE:-1}"
RUN_REF_COVERAGE="${RUN_REF_COVERAGE:-1}"
RUN_EXEC_SURFACE="${RUN_EXEC_SURFACE:-1}"
RUN_BENCH="${RUN_BENCH:-1}"

echo "[1/9] Building ROOT Minuit2 reference runner (${ROOT_TAG})"
./scripts/build_root_reference_runner.sh "${ROOT_TAG}"

echo "[2/9] Running ROOT-vs-Rust differential workloads"
python3 scripts/compare_ref_vs_rust.py

echo "[3/9] Refreshing traceability matrix and gate"
python3 scripts/generate_traceability_matrix.py
python3 scripts/check_traceability_gate.py --mode non-regression

if [[ "${RUN_TESTS}" == "1" ]]; then
  echo "[4/9] Running test suite"
  cargo test --no-default-features
  PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo test --all-features
else
  echo "[4/9] Skipping tests (RUN_TESTS=${RUN_TESTS})"
fi

if [[ "${RUN_COVERAGE}" == "1" ]]; then
  echo "[5/9] Regenerating coverage reports"
  mkdir -p reports/coverage
  cargo llvm-cov clean --workspace
  cargo llvm-cov --no-default-features --summary-only > reports/coverage/core_coverage_raw.txt
  cargo llvm-cov clean --workspace
  PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 cargo llvm-cov --all-features --summary-only > reports/coverage/all_features_coverage_raw.txt
  python3 scripts/generate_coverage_reports.py
else
  echo "[5/9] Skipping coverage (RUN_COVERAGE=${RUN_COVERAGE})"
fi

if [[ "${RUN_REF_COVERAGE}" == "1" ]]; then
  echo "[6/9] Regenerating C++ reference executed-surface coverage"
  python3 scripts/generate_reference_coverage.py --root-tag "${ROOT_TAG}"
else
  echo "[6/9] Skipping reference coverage (RUN_REF_COVERAGE=${RUN_REF_COVERAGE})"
fi

if [[ "${RUN_EXEC_SURFACE}" == "1" ]]; then
  echo "[7/9] Generating executed-surface mapping and gate"
  python3 scripts/generate_executed_surface_mapping.py
  python3 scripts/check_executed_surface_gate.py --mode non-regression
else
  echo "[7/9] Skipping executed-surface mapping (RUN_EXEC_SURFACE=${RUN_EXEC_SURFACE})"
fi

if [[ "${RUN_BENCH}" == "1" ]]; then
  echo "[8/9] Regenerating benchmark reports"
  mkdir -p reports/benchmarks
  cargo bench --bench benchmarks -- --noplot > reports/benchmarks/default_raw.txt
  cargo bench --features parallel --bench benchmarks -- --noplot > reports/benchmarks/parallel_raw.txt
  python3 scripts/generate_benchmark_report.py
else
  echo "[8/9] Skipping benchmarks (RUN_BENCH=${RUN_BENCH})"
fi

echo "[9/9] Generating claim scorecard"
python3 scripts/generate_verification_scorecard.py

echo "Verification completed."
echo "Scorecard: reports/verification/scorecard.md"
echo "Manifest:  reports/verification/manifest.json"
