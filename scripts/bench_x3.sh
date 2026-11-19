#!/usr/bin/env bash
# bench_x3.sh BASE HEAD
#
# Compares the criterion benchmarks of two git revisions of this clone.
# Run it on a quiet machine. BASE and HEAD are git revisions of the clone,
# for example "origin/main HEAD". Only committed revisions are measured.
#
# Steps:
#   1. Create two git worktrees, bench-work/base and bench-work/head.
#      Build each once for the default features and once for --features parallel.
#   2. Run ROUNDS rounds (the variable below, set to 3). Each round runs the
#      benchmark program in this order: base, head, base again.
#      The second base run is the control. Each run covers both feature sets.
#   3. Read the median point estimate of every benchmark id from
#      target/criterion/<id>/new/estimates.json of the worktree.
#      The id comes from benchmark.json in the same folder.
#   4. Print one table row per benchmark id and feature set, then a RESULT line.
#      The full log is written to bench-logs/bench_x3_<UTC date>.log.
#   5. Remove the two worktrees. Only worktrees created by this run are removed.
#
# Run time: a run with the default criterion timings takes about 30 minutes,
# plus the builds. BENCH_EXTRA_ARGS below shortens it for a check of the script.
#
# Formulas. Round r gives the median times base1_r, head_r and base2_r of one
# id and feature set. The table shows a mean over the rounds:
#   base1, head, base2 = mean over r of base1_r, head_r, base2_r
#   ratio              = mean over r of (head_r / base1_r)
#   control            = mean over r of (base2_r / base1_r)
#   max_ratio          = largest ratio over all rows
#   control_max_drift  = largest |control - 1| over all rows
# The number of rows N must equal the number of ids listed by
# "cargo bench --bench benchmarks -- --list" for the two feature sets together.
#
# RESULT line and exit code:
#   0  RESULT: OK rows=N max_ratio=R control_max_drift=D
#        max_ratio is at most 1.10 and control_max_drift is at most 0.05.
#   1  RESULT: FAIL ...     max_ratio is above 1.10. The line names the slowest benchmark.
#   2  usage error         wrong argument count or a revision that does not exist.
#   3  RESULT: INVALID ...  control_max_drift is above 0.05, so the machine was busy.
#                           This is checked before FAIL.
#   4  RESULT: ERROR ...    a build or run failed, or the ids do not match the listing.
#
# Environment variable:
#   BENCH_EXTRA_ARGS  extra arguments for criterion, appended after --noplot.
#                     Default: empty. Example for a quick check of the script:
#                     BENCH_EXTRA_ARGS="--warm-up-time 1 --measurement-time 2 --sample-size 10"
#
# Cargo runs with CARGO_BUILD_JOBS=2 under "nice -n 19". The script also uses
# jobscope when it is installed. It loads ~/.cargo/env when cargo is not on PATH.
set -euo pipefail

ROUNDS=3

USAGE="usage: bash scripts/bench_x3.sh BASE HEAD   (two git revisions of this clone, for example: origin/main HEAD)"

usage_error() {
    echo "$USAGE" >&2
    [ "$#" -eq 0 ] || echo "bench_x3.sh: $*" >&2
    exit 2
}

[ "$#" -eq 2 ] || usage_error "need exactly two arguments, got $#"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
git -C "$ROOT" rev-parse --git-dir >/dev/null 2>&1 || usage_error "$ROOT is not a git clone"

resolve() {
    case "$1" in
        -*) usage_error "not a revision: $1" ;;
    esac
    git -C "$ROOT" rev-parse --verify -q "$1^{commit}" || usage_error "not a revision of this clone: $1"
}
BASE_REV="$1"
HEAD_REV="$2"
BASE_SHA="$(resolve "$BASE_REV")"
HEAD_SHA="$(resolve "$HEAD_REV")"

WORK="$ROOT/bench-work"
BASE_WT="$WORK/base"
HEAD_WT="$WORK/head"
for wt in "$BASE_WT" "$HEAD_WT"; do
    if [ -e "$wt" ]; then
        echo "bench_x3.sh: $wt exists from an earlier run. Remove it with: git worktree remove --force $wt" >&2
        exit 4
    fi
done

if ! command -v cargo >/dev/null 2>&1 && [ -f "$HOME/.cargo/env" ]; then
    # shellcheck source=/dev/null
    . "$HOME/.cargo/env"
fi
command -v cargo >/dev/null 2>&1 || { echo "bench_x3.sh: cargo not found" >&2; exit 4; }
command -v python3 >/dev/null 2>&1 || { echo "bench_x3.sh: python3 not found" >&2; exit 4; }

export CARGO_BUILD_JOBS=2
unset CRITERION_HOME CARGO_TARGET_DIR

RUN=(nice -n 19)
if command -v jobscope >/dev/null 2>&1; then
    RUN=(jobscope nice -n 19)
fi

read -r -a EXTRA <<<"${BENCH_EXTRA_ARGS:-}"

STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOGDIR="$ROOT/bench-logs"
LOG="$LOGDIR/bench_x3_$STAMP.log"
DATA="$LOGDIR/bench_x3_$STAMP.data"
mkdir -p "$LOGDIR" "$DATA"
: >"$LOG"

say() { printf '%s\n' "$*" | tee -a "$LOG"; }

CREATED=()
# shellcheck disable=SC2329  # runs through the EXIT trap
cleanup() {
    local rc=$?
    trap - EXIT
    local wt
    for wt in ${CREATED[@]+"${CREATED[@]}"}; do
        git -C "$ROOT" worktree remove --force "$wt" >>"$LOG" 2>&1 || echo "bench_x3.sh: could not remove $wt" >&2
    done
    rmdir "$WORK" 2>/dev/null || true
    exit "$rc"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

cargo_in() { # worktree, then the cargo arguments
    local wt=$1
    shift
    (cd "$wt" && CARGO_TARGET_DIR="$wt/target" "${RUN[@]}" cargo "$@")
}

step() { # label, then a command whose output goes to the log
    local label=$1
    shift
    say "$(date -u +%H:%M:%S) $label"
    if ! "$@" >>"$LOG" 2>&1; then
        say "RESULT: ERROR step failed: $label (see $LOG)"
        exit 4
    fi
}

set_features() { # sets FEAT for the feature set named by $1
    if [ "$1" = parallel ]; then
        FEAT=(--locked --features parallel)
    else
        FEAT=(--locked)
    fi
}

say "bench_x3 $STAMP"
say "host: $(hostname -s 2>/dev/null || hostname)"
say "cpu: $(awk -F': *' '/model name/ {print $2; exit}' /proc/cpuinfo 2>/dev/null || true)"
say "rustc: $(rustc --version)"
say "load at start: $(cut -d' ' -f1-3 /proc/loadavg 2>/dev/null || true)"
say "base: $BASE_REV = $BASE_SHA"
say "head: $HEAD_REV = $HEAD_SHA"
say "criterion extra arguments: ${BENCH_EXTRA_ARGS:-none}"
say "log: $LOG"

step "worktree base" git -C "$ROOT" worktree add --detach "$BASE_WT" "$BASE_SHA"
CREATED+=("$BASE_WT")
step "worktree head" git -C "$ROOT" worktree add --detach "$HEAD_WT" "$HEAD_SHA"
CREATED+=("$HEAD_WT")

for who in base head; do
    wt="$WORK/$who"
    for feat in default parallel; do
        set_features "$feat"
        step "build $who $feat" cargo_in "$wt" build --release --benches "${FEAT[@]}"
    done
done

for feat in default parallel; do
    set_features "$feat"
    say "list ids $feat"
    if ! cargo_in "$HEAD_WT" bench --bench benchmarks "${FEAT[@]}" -- --list >"$DATA/list_$feat.txt" 2>>"$LOG"; then
        say "RESULT: ERROR could not list the benchmark ids for the $feat features (see $LOG)"
        exit 4
    fi
done

for round in $(seq 1 "$ROUNDS"); do
    for pass in "base1 $BASE_WT" "head $HEAD_WT" "base2 $BASE_WT"; do
        name=${pass%% *}
        wt=${pass#* }
        for feat in default parallel; do
            set_features "$feat"
            step "bench round $round $name $feat" cargo_in "$wt" bench --bench benchmarks "${FEAT[@]}" -- --noplot ${EXTRA[@]+"${EXTRA[@]}"}
            # Keep this pass apart: the next pass in the same worktree writes target/criterion again.
            mv "$wt/target/criterion" "$wt/target/criterion.r$round.$name.$feat"
        done
    done
done

say "load at end: $(cut -d' ' -f1-3 /proc/loadavg 2>/dev/null || true)"
say "elapsed: ${SECONDS} s"

rc=0
python3 -I - "$BASE_WT" "$HEAD_WT" "$DATA/list_default.txt" "$DATA/list_parallel.txt" "$DATA/rows.json" "$ROUNDS" <<'PY' | tee -a "$LOG" || rc=$?
import json
import os
import sys

base_wt, head_wt, list_default, list_parallel, out_json = sys.argv[1:6]
ROUNDS = int(sys.argv[6])
MAX_RATIO = 1.10
MAX_DRIFT = 0.05
RUNS = (("base1", base_wt), ("head", head_wt), ("base2", base_wt))
FEATURES = (("default", list_default), ("parallel", list_parallel))
SUFFIX = ": benchmark"


def listed(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.endswith(SUFFIX):
                ids.append(line[: -len(SUFFIX)])
    return ids


def medians(snapshot):
    found = {}
    for name in sorted(os.listdir(snapshot)):
        folder = os.path.join(snapshot, name, "new")
        estimates = os.path.join(folder, "estimates.json")
        if not os.path.isfile(estimates):
            continue
        bench_id = name
        meta = os.path.join(folder, "benchmark.json")
        if os.path.isfile(meta):
            with open(meta) as fh:
                bench_id = json.load(fh).get("full_id") or name
        with open(estimates) as fh:
            found[bench_id] = float(json.load(fh)["median"]["point_estimate"])
    return found


def mean(values):
    return sum(values) / len(values)


rows = []
problems = []
listed_total = 0
for feature, list_path in FEATURES:
    ids = listed(list_path)
    listed_total += len(ids)
    data = {}
    for rnd in range(1, ROUNDS + 1):
        for run, wt in RUNS:
            snapshot = os.path.join(wt, "target", "criterion.r%d.%s.%s" % (rnd, run, feature))
            try:
                data[(rnd, run)] = medians(snapshot)
            except (OSError, ValueError, KeyError) as err:
                problems.append("%s round %d %s: cannot read %s (%s)" % (feature, rnd, run, snapshot, err))
                data[(rnd, run)] = {}
            if data[(rnd, run)] and set(data[(rnd, run)]) != set(ids):
                problems.append("%s round %d %s: measured ids differ from the listed ids" % (feature, rnd, run))
    for bench_id in ids:
        base1s, heads, base2s = [], [], []
        for rnd in range(1, ROUNDS + 1):
            values = [data[(rnd, run)].get(bench_id) for run, _ in RUNS]
            if any(v is None or not v > 0.0 for v in values):
                break
            base1s.append(values[0])
            heads.append(values[1])
            base2s.append(values[2])
        if len(heads) != ROUNDS:
            problems.append("%s: no usable median for %s" % (feature, bench_id))
            continue
        ratios = [h / b for h, b in zip(heads, base1s)]
        controls = [b2 / b1 for b2, b1 in zip(base2s, base1s)]
        rows.append(
            {
                "features": feature,
                "id": bench_id,
                "base1_ns": mean(base1s),
                "head_ns": mean(heads),
                "base2_ns": mean(base2s),
                "ratio": mean(ratios),
                "control": mean(controls),
                "ratio_rounds": ratios,
                "control_rounds": controls,
            }
        )

if len(rows) != listed_total and not problems:
    problems.append("rows=%d but the listings hold %d ids" % (len(rows), listed_total))

if problems:
    for text in problems:
        print(text)
    print("RESULT: ERROR rows=%d listed=%d (see the lines above)" % (len(rows), listed_total))
    sys.exit(4)

width = max(len(r["id"]) for r in rows)
print("%-8s  %-*s  %10s  %10s  %10s  %7s  %8s" % ("features", width, "id", "base1 ns", "head ns", "base2 ns", "ratio", "control"))
for r in rows:
    print(
        "%-8s  %-*s  %10.1f  %10.1f  %10.1f  %7.3f  %8.3f"
        % (r["features"], width, r["id"], r["base1_ns"], r["head_ns"], r["base2_ns"], r["ratio"], r["control"])
    )

slowest = max(rows, key=lambda r: r["ratio"])
drifting = max(rows, key=lambda r: abs(r["control"] - 1.0))
max_ratio = slowest["ratio"]
max_drift = abs(drifting["control"] - 1.0)
head = "rows=%d max_ratio=%.3f control_max_drift=%.3f" % (len(rows), max_ratio, max_drift)
tail = 'slowest="%s" (%s features)' % (slowest["id"], slowest["features"])

if max_drift > MAX_DRIFT:
    status, code = "INVALID", 3
    line = 'RESULT: INVALID %s %s; the machine was busy, control drift at "%s" (%s features)' % (
        head, tail, drifting["id"], drifting["features"])
elif max_ratio > MAX_RATIO:
    status, code = "FAIL", 1
    line = "RESULT: FAIL %s %s" % (head, tail)
else:
    status, code = "OK", 0
    line = "RESULT: OK %s" % head

with open(out_json, "w") as fh:
    json.dump({"status": status, "rows": rows, "max_ratio": max_ratio, "control_max_drift": max_drift}, fh, indent=1)
print(line)
sys.exit(code)
PY
exit "$rc"
