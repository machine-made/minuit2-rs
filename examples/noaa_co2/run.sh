#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"
SECONDS=0

if [[ ! -f "examples/data/noaa/co2_mm_mlo.csv" ]]; then
  scripts/fetch_scientific_demo_data.sh
fi

cargo run --example noaa_co2
python3 "${SCRIPT_DIR}/plot.py"

elapsed=$SECONDS
printf "Total execution time: %02dh:%02dm:%02ds\n" "$((elapsed/3600))" "$(((elapsed%3600)/60))" "$((elapsed%60))"
echo "Generated:"
echo "  examples/noaa_co2/figures/noaa_co2_fit.png"
echo "  examples/noaa_co2/figures/noaa_co2_residual_hist.png"
