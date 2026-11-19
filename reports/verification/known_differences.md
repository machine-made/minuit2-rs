# Known Differences Ledger (ROOT `v6-36-08` vs `minuit2-rs`)

This ledger records currently observed differences from differential verification runs.

## D-001: Waived NFCN divergence on quadratic fixed-parameter workloads

- Workloads:
  - `quadratic3_fixx_migrad`
  - `quadratic3_fixx_hesse`
- Classification: `algorithmic/path-efficiency`
- Severity: `low` (correctness metrics pass; NFCN-only waiver)
- Expected: `yes`
- Evidence:
  - `reports/verification/diff_summary.md`
  - `reports/verification/diff_results.csv`
  - `verification/workloads/root_minuit2_v6_36_08.json` (`nfcn_rel_waiver` on the two workloads only)
- ROOT v6-36-08 source cited:
  - `third_party/root_ref/math/minuit2/src/MnSeedGenerator.cxx`: numerical seed generation first evaluates the starting point with `MnFcnCaller{fcn}(x)` and then calls `gc(pa)`.
  - `third_party/root_ref/math/minuit2/src/Numerical2PGradientCalculator.cxx`: `operator()(par)` constructs an `InitialGradientCalculator` and then central-difference probes each variable parameter via `mfcnCaller(x)` at `xtf + step` and `xtf - step` for up to `GradientNCycles()` cycles.
  - `third_party/root_ref/math/minuit2/src/InitialGradientCalculator.cxx`: the initial rough gradient itself is heuristic and does not call the FCN.
- Details:
  - The counter test in `tests/nfcn_accounting.rs` checks running totals on a synthetic fixed-parameter function.
    It does not check this runner's 24-call and 34-call totals.
  - Migrad uses 24 calls against ROOT's 29. Rust stops at EDM `1.492603e-05`; ROOT reaches `1.072811e-25`.
  - Hesse uses 34 calls against ROOT's 39. The Rust count includes 24 Migrad calls and 10 Hesse calls.
  - Analytical gradients now receive numerical curvature and steps before numerical Hesse, following `MnHesse.cxx:202-210`.
  - The runner uses analytical gradients in Rust. The remaining gap includes different seed work and stopping EDM; the old counter reset does not explain it.
  - Correctness metrics pass tolerances, including `quadratic3_fixx_migrad` parameter max abs diff `3.886512e-04` and `quadratic3_fixx_hesse` covariance max abs diff `5.869853e-10` in the current differential report.
- Disposition: explicit per-workload waiver; no global NFCN threshold change.
- Next step:
  - Retain the existing manifest waivers. Isolate the remaining seed and stopping differences before further changes.

## D-003: Residual Rosenbrock covariance differences

- Workloads and measured maximum absolute covariance differences:
  - `rosenbrock2_migrad`: `1.771473e-09`
  - `rosenbrock2_minimize`: `1.771473e-09`
  - `rosenbrock2_migrad_strategy2`: `1.106493e-11`
- Classification: `numerical/residual`
- Severity: `low`; all covariance differences remain below the configured limit of 0.5.
- Evidence: `reports/verification/diff_results.csv` and the stored outputs in `reports/verification/raw/ref/` and `reports/verification/raw/rust/`.
- Details: the strategy-2 seed Hesse, triangular covariance-update sum, and line-search corrections closed the former 0.13–0.18 gap.
- Disposition: retain this entry for the measured residuals; their precise rounding source remains unisolated.

## Resolved since previous snapshot

- Upper-limited Migrad now follows ROOT's negative-curvature seed step at a zero gradient.
  `NegativeG2LineSearch.cxx:80-83` uses the signed gradient step when curvature is negative.
  The first differing evaluation was call 10, during seed generation.
  Rust formerly probed `[0.3999999999999999, -1.9142135623730951]`; ROOT probed `[0.4, -1.5000001261539815]`.
  The upper transform now retains rounding residuals through ROOT's expression, evaluated in `long double` at `SqrtUpParameterTransformation.cxx:22-26`.
  Both runs use 61 calls, including 39 seed calls and 10 Hesse calls.
  All call vectors agree within relative `1e-9`; the measured maximum difference is `2.3e-16`.
  The sine derivative retains its sign, following `SinParameterTransformation.cxx:53-56`, fixing covariance signs on the falling branch.
  The manifest now compares the two-sided, upper-limited, and limited-scan workloads, including errors and covariance.
  Every added numerical field differs from the stored ROOT output by less than `1e-9` absolutely.
  The live reference binary differs from those stored files in parameters, errors, covariance, function value, and EDM.
  Those residuals remain below `1e-9` absolutely; the stored reference files remain the oracle.

- D-005: the seed gradient now keeps its curvature floor fixed across each parameter's cycles, following `Numerical2PGradientCalculator.cxx:115-151`.
  Hesse refinement starts from the diagonal sweep and follows `HessianGradientCalculator.cxx:68-132` for steps and stopping.
  Rosenbrock now uses 140 calls on both sides, including four refinement calls.
  `tests/call_trace_parity.rs` compares all 140 vectors with the stored ROOT trace at relative tolerance `1e-9`.
  Hahn1 now uses 456 calls on both sides, including fourteen refinement calls.
  Its Hesse pass starts after 294 calls; both continuations take 112 calls.
  The Rust traces under `raw/trace/` were regenerated twice with identical results and verified user call totals.
- D-006: both Minos rows now report 38 calls before Minos runs.
  Both runners make 28 Migrad calls and 10 standalone Hesse calls.
  Each Hesse pass uses one centre, four diagonal, four refinement and one mixed evaluation.
  Evidence: `tools/root_evidence/nfcn_accounting.rs` and the matching stored ROOT and Rust traces under `raw/trace/`.

- Previous `NoG2CallsWhenFCHasNoG2` P0 gap resolved by adding optional FCN `has_g2()/g2()` and `has_hessian()/hessian()` contract and porting the regression as `tests/root_reference_minuit2.rs::root_no_g2_calls_when_fcn_has_no_g2`.
- Previous `quadratic2_simplex` correctness fail resolved after aligning Simplex EDM tolerance semantics with ROOT (`minedm = tolerance * Up`).
- Previous `rosenbrock2_minimize` correctness fail resolved after aligning `MnMinimize` flow with ROOT combined minimizer logic and propagating covariance from final minimum state.

## Reference trace evidence

- `nfcn_accounting_quadratic2_upper_limited_migrad.{root,rust}.json` records every upper-limited evaluation and its phase.
  The headers of `tools/root_evidence/nfcn_accounting.{cpp,rs}` give the build and run commands.
  `tests/limited_call_trace_parity.rs` checks all 61 calls against the stored ROOT trace.

- The stored ROOT trace of the strategy-2 Rosenbrock run is generated by `tools/ref_runner_cpp/main.cpp` using `MnTraceObject`, attached through `Builder().SetTraceObject`.
  The program's header gives build and run commands for `reports/verification/raw/trace/rosenbrock2_migrad_strategy2.root.jsonl`.
  It uses the reference workload's objective, start, steps, and strategy.
  Iteration 1 gives `dcovar = 0.16098749348302552`; successful seed Hesse gives exactly zero.
- The stored direct ROOT v6-36-08 line-search evidence is generated by `tools/root_evidence/linesearch_open_step.cpp`.
  Its header gives build and run commands for `reports/verification/raw/trace/linesearch_open_step.root.json`.
  It calls `MnLineSearch` on the quartic in `tests/linesearch_open_step.rs`.
  The stored start, direction, and directional derivative reproduce lambda `10.39095385821868` and fval `-2191.6554284992781`.
  Both records use the existing reference Minuit2 library.

## Strict ULP report

The report covers 15 workloads; all three added workloads pass.
Four strict mismatches also appear in the previous report, which covered 12 workloads.
The fixed-parameter call counts are ROOT/Rust `29/24` and `39/34`, within their existing manifest warning thresholds.
The other two mismatches are covariance entries `0` versus `-0.0` in `quadratic2_no_g2_migrad`.
The report requires exact integer equality and matching zero signs independently of numerical tolerances.
The differential comparison passes all 15 workloads; the strict ULP command still exits with status 1.
