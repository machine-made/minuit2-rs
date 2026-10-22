#!/usr/bin/env bash
set -euo pipefail

ROOT_TAG="${1:-v6-36-08}"
ROOT_REF_DIR="${ROOT_REF_DIR:-third_party/root_ref}"
ROOT_BUILD_DIR="${ROOT_BUILD_DIR:-third_party/root_ref_build/minuit2}"
RUNNER_BUILD_DIR="${RUNNER_BUILD_DIR:-third_party/root_ref_build/ref_runner}"
ROOT_CMAKE_EXTRA_ARGS="${ROOT_CMAKE_EXTRA_ARGS:-}"
RUNNER_CMAKE_EXTRA_ARGS="${RUNNER_CMAKE_EXTRA_ARGS:-}"

REPO_URL="https://github.com/root-project/root.git"

mkdir -p "$(dirname "$ROOT_REF_DIR")"

if [[ ! -d "$ROOT_REF_DIR/.git" ]]; then
  git clone --filter=blob:none --no-checkout "$REPO_URL" "$ROOT_REF_DIR"
fi

git -C "$ROOT_REF_DIR" fetch --tags --force origin "$ROOT_TAG"
git -C "$ROOT_REF_DIR" sparse-checkout init --no-cone
# Standalone Minuit2 references these paths in ROOT tree layout.
git -C "$ROOT_REF_DIR" sparse-checkout set \
  /math/minuit2 \
  /math/mathcore/inc/Fit \
  /math/mathcore/inc/Math \
  /math/mathcore/src \
  /core/foundation/inc/ROOT \
  /rootx/src/rootx.cxx \
  /LICENSE \
  /LGPL2_1.txt

git -C "$ROOT_REF_DIR" checkout -f "$ROOT_TAG"

rm -rf "$ROOT_BUILD_DIR"
cmake -S "$ROOT_REF_DIR/math/minuit2" -B "$ROOT_BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -Dminuit2_inroot=ON -Dminuit2_standalone=ON ${ROOT_CMAKE_EXTRA_ARGS}
cmake --build "$ROOT_BUILD_DIR" --target Minuit2 -j

ROOT_BUILD_DIR_ABS="$(cd "$ROOT_BUILD_DIR" && pwd)"

rm -rf "$RUNNER_BUILD_DIR"
cmake -S tools/ref_runner_cpp -B "$RUNNER_BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -DMinuit2_DIR="$ROOT_BUILD_DIR_ABS" ${RUNNER_CMAKE_EXTRA_ARGS}
cmake --build "$RUNNER_BUILD_DIR" --target ref_runner -j

echo "$RUNNER_BUILD_DIR/ref_runner"
