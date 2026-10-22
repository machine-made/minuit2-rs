#!/usr/bin/env bash
#
# Fetch and verify scientific demo reference data.
#
# Every file is pinned against examples/data/SHA256SUMS. A run on an intact
# tree performs zero network traffic. Any checksum mismatch, missing manifest
# entry, TLS failure, or partial download aborts with a non-zero exit code.
# Partial downloads are staged to a sibling tempfile so a corrupt fetch can
# never poison examples/data/.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT_DIR}/examples/data"
SUMS_FILE="${DATA_DIR}/SHA256SUMS"

if [[ ! -f "$SUMS_FILE" ]]; then
  echo "error: checksum manifest not found: $SUMS_FILE" >&2
  exit 1
fi

# Portable SHA-256 tool. Both print "<hex-hash>  <path>".
if command -v sha256sum >/dev/null 2>&1; then
  SUM_OF=(sha256sum)
elif command -v shasum >/dev/null 2>&1; then
  SUM_OF=(shasum -a 256)
else
  echo "error: neither sha256sum nor shasum available in PATH" >&2
  exit 1
fi

mkdir -p \
  "${DATA_DIR}/noaa" \
  "${DATA_DIR}/nist" \
  "${DATA_DIR}/usgs" \
  "${DATA_DIR}/cern"

PATHS=(
  "noaa/co2_mm_mlo.csv"
  "nist/Misra1a.dat"
  "nist/Hahn1.dat"
  "nist/Rat43.dat"
  "nist/Misra1b.dat"
  "nist/Chwirut2.dat"
  "nist/Lanczos3.dat"
  "nist/Gauss1.dat"
  "nist/Kirby2.dat"
  "nist/Thurber.dat"
  "nist/Rat42.dat"
  "nist/MGH09.dat"
  "nist/BoxBOD.dat"
  "nist/ENSO.dat"
  "usgs/earthquakes_2025_m4p5.csv"
  "cern/MuRun2010B_0.csv"
  "cern/Zmumu.csv"
)

url_for() {
  case "$1" in
    noaa/co2_mm_mlo.csv)
      printf '%s\n' "https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_mlo.csv" ;;
    nist/Misra1a.dat)
      printf '%s\n' "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Misra1a.dat" ;;
    nist/Hahn1.dat)
      printf '%s\n' "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Hahn1.dat" ;;
    nist/Rat43.dat)
      printf '%s\n' "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Rat43.dat" ;;
    nist/Misra1b.dat)
      printf '%s\n' "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Misra1b.dat" ;;
    nist/Chwirut2.dat)
      printf '%s\n' "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Chwirut2.dat" ;;
    nist/Lanczos3.dat)
      printf '%s\n' "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Lanczos3.dat" ;;
    nist/Gauss1.dat)
      printf '%s\n' "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Gauss1.dat" ;;
    nist/Kirby2.dat)
      printf '%s\n' "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Kirby2.dat" ;;
    nist/Thurber.dat)
      printf '%s\n' "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Thurber.dat" ;;
    nist/Rat42.dat)
      printf '%s\n' "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Rat42.dat" ;;
    nist/MGH09.dat)
      printf '%s\n' "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/MGH09.dat" ;;
    nist/BoxBOD.dat)
      printf '%s\n' "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/BoxBOD.dat" ;;
    nist/ENSO.dat)
      printf '%s\n' "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/ENSO.dat" ;;
    usgs/earthquakes_2025_m4p5.csv)
      printf '%s\n' "https://earthquake.usgs.gov/fdsnws/event/1/query.csv?starttime=2025-01-01&endtime=2025-12-31&minmagnitude=4.5&orderby=time" ;;
    cern/MuRun2010B_0.csv)
      printf '%s\n' "https://eospublic.cern.ch/eos/opendata/cms/Run2010B/Mu/CSV/Apr21ReReco-v1/MuRun2010B_0.csv" ;;
    cern/Zmumu.csv)
      printf '%s\n' "https://eospublic.cern.ch/eos/opendata/cms/Run2011A/DoubleMu/CSV/12Oct2013-v1/Zmumu.csv" ;;
    *)
      echo "error: no URL mapped for $1" >&2; return 1 ;;
  esac
}

# Look up the expected hash for a path from SHA256SUMS. Matches on the last
# whitespace-separated field so it is immune to a single- vs two-space
# sha256sum/shasum format difference.
expected_hash_for() {
  local rel="$1"
  awk -v target="$rel" '
    NF >= 2 && $NF == target { print $1; found = 1; exit }
    END { if (!found) exit 1 }
  ' "$SUMS_FILE"
}

actual_hash_of_path() {
  # Accepts an absolute path to a file (doesn't have to live inside DATA_DIR).
  "${SUM_OF[@]}" "$1" | awk '{print $1}'
}

for rel in "${PATHS[@]}"; do
  abs="${DATA_DIR}/${rel}"

  expected="$(expected_hash_for "$rel")" || {
    echo "error: no checksum entry for $rel in ${SUMS_FILE}" >&2
    exit 2
  }

  if [[ -f "$abs" ]] && [[ "$(actual_hash_of_path "$abs")" == "$expected" ]]; then
    echo "ok    ${rel}  (cached)"
    continue
  fi

  url="$(url_for "$rel")"
  echo "fetch ${rel}"
  echo "      ${url}"

  tmp="$(mktemp "${abs}.XXXXXX")"
  trap 'rm -f "$tmp"' EXIT

  # TLS is verified (no -k). --proto '=https' refuses any non-HTTPS scheme,
  # including on redirects. --tlsv1.2 sets a sane minimum TLS version.
  curl \
    --fail --location --silent --show-error \
    --proto '=https' \
    --tlsv1.2 \
    --retry 3 --retry-delay 2 \
    "$url" -o "$tmp"

  # Verify the download BEFORE publishing it, so a bad fetch can never
  # poison $abs. Only a hash-matching tempfile is ever moved into place.
  actual="$(actual_hash_of_path "$tmp")"
  if [[ "$expected" != "$actual" ]]; then
    echo "error: checksum mismatch for $rel" >&2
    echo "       expected $expected" >&2
    echo "       actual   $actual" >&2
    echo "       (previous ${abs} left untouched)" >&2
    exit 2
  fi

  mv -f "$tmp" "$abs"
  trap - EXIT
  echo "ok    ${rel}  (verified)"
done

echo
echo "All data files present and verified against examples/data/SHA256SUMS."
