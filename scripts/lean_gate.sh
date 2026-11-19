#!/usr/bin/env bash
set -euo pipefail

repo=$(realpath "$(dirname "${BASH_SOURCE[0]}")/..")
project="$repo/proofs/lean"
scratch="$repo/target/lean-scratch"
mkdir -p "$scratch"

version=$(lean --version)
printf '%s\n' "$version"
if [[ ! "$version" =~ ^Lean\ \(version\ 4\.35\.0, ]]; then
    printf '%s\n' 'Lean 4.35.0 is required.' >&2
    exit 1
fi

# Keep downloads inside the worktree and run heavy steps in sequence.
export MATHLIB_CACHE_DIR="$project/.lake/cache"
export RAYON_NUM_THREADS=2
export CARGO_BUILD_JOBS=2
env -C "$project" jobscope lake exe cache get | tee "$scratch/cache.log"
env -C "$project" jobscope lake build | tee "$scratch/build.log"

if grep -rnE '\b(sorry|admit|axiom)\b' "$project" --include='*.lean' --exclude-dir=.lake; then
    printf '%s\n' 'Proof integrity check failed.' >&2
    exit 1
else
    grep_status=$?
    if [[ "$grep_status" -ne 1 ]]; then
        exit "$grep_status"
    fi
fi

env -C "$project" jobscope lake env lean "$project/CheckAxioms.lean" | tee "$scratch/axioms.log"
python3 - "$scratch/axioms.log" <<'PY'
import re
import sys
from pathlib import Path

output = Path(sys.argv[1]).read_text()
expected = {
    'Minuit2Proofs.dfpUpdate_posDef',
    'Minuit2Proofs.dfpUpdate_secant',
    'Minuit2Proofs.dfpUpdate_neg_delgam_not_posDef',
}
allowed = {'propext', 'Classical.choice', 'Quot.sound'}
entries = re.findall(r"'([^']+)' depends on axioms:\s*\[([^\]]*)\]", output)
if len(entries) != 3 or {name for name, _ in entries} != expected:
    sys.exit('Missing or unexpected theorem dependency output.')
for name, dependencies in entries:
    found = {part.strip() for part in dependencies.split(',') if part.strip()}
    if not found <= allowed:
        sys.exit(f'{name}: unexpected dependencies: {found - allowed}')
    print('axioms ok')
PY
