# NegativeG2LineSearch parity audit (ROOT v6-36-08)

Reference checkout: `third_party/root_ref` at `a8ca1b23e38d7dbe0ff24027894ca0f2ad65f1bd` (`v6-36-08`).

## Findings table

| Call site | ROOT behavior | Rust behavior | Gap? |
|---|---|---|---|
| `MnSeedGenerator::operator()` numerical seed (`math/minuit2/src/MnSeedGenerator.cxx:43-125`, calls at `:99-105`) | Builds initial state from numerical gradient, then, when no user covariance exists, calls `NegativeG2LineSearch` if any `G2 <= 0`. The line search repairs before optional strategy-2 Hesse (`:110-119`). | `MigradSeedGenerator::generate` computes numerical seed and calls `escape_negative_curvature` before building `V0`. | Implemented gap: old Rust repaired all negative-G2 coordinates in one vector step per pass. It now repairs one coordinate, recomputes full gradient, then scans again, matching ROOT's per-coordinate control flow (`NegativeG2LineSearch.cxx:62-123`). Direction split: ROOT-reachable `grad < 0` uses `+gstep`; ROOT-reachable `grad >= 0` with `|g2| > eps` uses `-gstep`. Intentional micro-waiver: for the start-at-limit/zero-Jacobian case that ROOT skips (`NegativeG2LineSearch.cxx:68-71`) and therefore defines no direction for, Rust keeps the positive escape direction. At an exact transform singularity, the positive probe is floored to a finite internal step so `tests/limit_boundary.rs` does not reconverge at the bound before the variable-metric builder has usable curvature. |
| `MnSeedGenerator::CallWithAnalyticalGradientCalculator` analytical seed (`math/minuit2/src/MnSeedGenerator.cxx:127-248`, calls at `:228-234`) | If analytical calculator can provide G2/Hessian, builds seed and still runs `NegativeG2LineSearch` when any `G2 <= 0` and no covariance exists. If it cannot provide G2, falls back to numerical seed (`:138-142`). | `generate_with_gradient` has the same negative-G2 guard and analytical line-search path, but current `AnalyticalGradientCalculator` synthesizes positive heuristic G2 rather than consuming user-supplied external G2/Hessian in seed generation. | Intentional documented gap for this bead: no observable negative analytical G2 reaches seed.rs today without changing `src/gradient/analytical.rs` (outside the allow-list). The guard/path is in place for traceability. |
| `NegativeG2LineSearch` algorithm (`math/minuit2/src/NegativeG2LineSearch.cxx:28-140`; declaration `inc/Minuit2/NegativeG2LineSearch.h:29-35`) | No-op if no `G2 <= 0` (`:42-44`). For each pass, finds the first non-positive G2 coordinate, optionally skips near-zero grad+G2, line-searches only that coordinate with numerical `Gstep` or analytical unit step (`:62-93`), updates params, recomputes gradients/G2 (`:97-114`), and repeats. Rebuilds diagonal covariance from G2 and EDM (`:125-140`). | `escape_negative_curvature_with` now mirrors first-offender coordinate repair and full-gradient recomputation; seed construction then rebuilds diagonal covariance and EDM from repaired G2. | Implemented, with the direction split and start-at-limit micro-waiver above. |
| `VariableMetricBuilder::Minimum` iteration loop (`math/minuit2/src/VariableMetricBuilder.cxx:205-375`) | No include or call to `NegativeG2LineSearch`; after line search it recomputes gradient (`:297-300`), estimates EDM (`:301`), fixes non-posdef matrix if needed (`:309-321`), then applies error update (`:323-330`). Top-level builder may Hesse/reiterate based on EDM/matrix quality (`:134-181`), but not negative-G2-specific. | `VariableMetricBuilder::iterate` recomputes gradients after line search and handles non-descent/non-posdef via `make_pos_def`, but has no negative-G2 line search. | No gap: ROOT has no in-iteration NegativeG2LineSearch call in v6-36-08. |
| `MnApplication` / Migrad / Minimize strategy retries (`src/MnApplication.cxx:27-53`; `src/ModularFunctionMinimizer.cxx:33-59`; `src/CombinedMinimumBuilder.cxx:27-41`; `inc/Minuit2/VariableMetricMinimizer.h:41-52`; `inc/Minuit2/CombinedMinimizer.h:30-39`) | `MnApplication` delegates to the selected minimizer. `ModularFunctionMinimizer` constructs the gradient calculator, calls the minimizer's `SeedGenerator` once (`ModularFunctionMinimizer.cxx:57-59`), then builder. `CombinedMinimumBuilder` retries Migrad after Simplex by creating a new VM seed (`CombinedMinimumBuilder.cxx:39-41`), so any negative-G2 handling happens through `MnSeedGenerator`, not in the application/retry layer. | `MnMigrad` and `MnMinimize` construct seeds before builder passes; negative-G2 handling is seed-local. | No additional gap for in-iteration/application retry handling. |

## Regression test

Added integration regression `negative_g2_seed_line_search_repairs_one_coordinate_before_recomputing_gradient` in `tests/negative_g2.rs`.

Pre-change equivalent failure (old all-negative-components vector repair moved both coordinates before recomputing):

```text
thread 'negative_g2_seed_line_search_repairs_one_coordinate_before_recomputing_gradient' panicked at tests/negative_g2.rs:...:
assertion `left == right` failed: ROOT repairs only one coordinate before recomputing the full gradient
  left: 0.1
 right: 0.0
```

Current run:

```text
cargo test --all-features --test negative_g2
running 1 test
test negative_g2_seed_line_search_repairs_one_coordinate_before_recomputing_gradient ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
```

## Verification log

Required gates:

```text
bash .sc/minuit2-rs-k5h.gate.sh
== targeted tests ==

test result: ok. 10 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s

     Running tests/robustness.rs (target/debug/deps/robustness-d6f3a296f5f7c4a2)

running 6 tests
test boundary_edge_case ... ok
test goldstein_price ... ok
test inf_resilience ... ok
test nan_resilience ... ok
test rosenbrock_hard_start ... ok
test high_dim_stress ... ok

test result: ok. 6 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s

== negative_g2 tests ==
running 1 test
test negative_g2_seed_line_search_repairs_one_coordinate_before_recomputing_gradient ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s

== clippy ==
    Checking minuit2 v0.5.1 (/home/rfrantz/Projects/minuit2-rs)
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 1.51s
== diff_results.csv status ==
```

```text
cargo test --all-features
...
     Running tests/limit_boundary.rs (target/debug/deps/limit_boundary-fd9c2ea6456b18c3)

running 3 tests
test lower_limited_started_at_lower_bound_escapes ... ok
test doubly_limited_started_at_lower_bound_escapes ... ok
test upper_limited_started_at_upper_bound_escapes ... ok

test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
...
     Running tests/negative_g2.rs (target/debug/deps/negative_g2-8b1dc33806096911)

running 1 test
test negative_g2_seed_line_search_repairs_one_coordinate_before_recomputing_gradient ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
...
Doc-tests minuit2

running 2 tests
test src/lib.rs - (line 20) ... ok
test src/lib.rs - (line 5) ... ok

test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
```

```text
python3 scripts/compare_ref_vs_rust.py
Wrote reports/verification/diff_results.csv
Wrote reports/verification/diff_summary.md
Status counts: pass=12 warn=0 fail=0
```

```text
python3 scripts/check_executed_surface_gate.py --mode non-regression
Current summary: P0=0 P1=48 P2=425
Baseline summary: P0=0 P1=48 P2=425
PASS: non-regression executed-surface gate
```
