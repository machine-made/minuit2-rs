#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"
SECONDS=0

if [[ ! -f "examples/data/cern/MuRun2010B_0.csv" ]]; then
  scripts/fetch_scientific_demo_data.sh
fi

cargo run --example cern_dimuon
python3 "${SCRIPT_DIR}/plot.py"

elapsed=$SECONDS
printf "Total execution time: %02dh:%02dm:%02ds\n" "$((elapsed/3600))" "$(((elapsed%3600)/60))" "$((elapsed%60))"
echo "Generated:"
echo "  examples/cern_dimuon/figures/cern_murun_jpsi_fit.png"
echo "  examples/cern_dimuon/figures/cern_zmumu_z_fit.png"
