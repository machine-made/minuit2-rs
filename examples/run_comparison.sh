#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

BATCH_REPEATS="${1:-41}"
WARMUPS="${2:-9}"
BOOT_ITERS="${3:-5000}"
PERM_ITERS="${4:-10000}"
BATCHES="${5:-9}"
CPU_CORE="${6:-}"
TRIM_FRACTION="${7:-0.10}"
STRICT_ENV="${8:-0}"
SECONDS=0

if [[ ! -f "examples/data/noaa/co2_mm_mlo.csv" ]]; then
  scripts/fetch_scientific_demo_data.sh
fi

# Build Rust examples in release mode for stable timing.
cargo build --release --examples

# Ensure Minuit2 C++ reference library exists.
if [[ ! -f "third_party/root_ref_build/minuit2/Minuit2Config.cmake" ]]; then
  scripts/build_root_reference_runner.sh v6-36-08
fi

MINUIT2_DIR_ABS="$(cd third_party/root_ref_build/minuit2 && pwd)"
cmake -S tools/scientific_runner_cpp \
  -B third_party/root_ref_build/scientific_runner \
  -DCMAKE_BUILD_TYPE=Release \
  -DMinuit2_DIR="${MINUIT2_DIR_ABS}"
cmake --build third_party/root_ref_build/scientific_runner --target scientific_runner -j

PY_ARGS=(
  --repeats "${BATCH_REPEATS}"
  --warmups "${WARMUPS}"
  --batches "${BATCHES}"
  --bootstrap-iters "${BOOT_ITERS}"
  --permutation-iters "${PERM_ITERS}"
  --trim-fraction "${TRIM_FRACTION}"
)

if [[ -n "${CPU_CORE}" ]]; then
  PY_ARGS+=(--cpu-core "${CPU_CORE}")
fi

if [[ "${STRICT_ENV}" == "1" ]]; then
  PY_ARGS+=(--strict-env)
fi

python3 examples/compare_timings.py "${PY_ARGS[@]}"

elapsed=$SECONDS
printf "Total comparison execution time: %02dh:%02dm:%02ds\n" "$((elapsed/3600))" "$(((elapsed%3600)/60))" "$((elapsed%60))"
echo "Generated:"
echo "  examples/comparison.png"
echo "  examples/output/comparison_timings.csv"
echo "  examples/output/comparison_samples.csv"
echo "  examples/output/comparison_environment.json"
