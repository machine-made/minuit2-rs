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
    It does not check this runner's 22-call and 36-call totals.
  - Migrad uses 22 calls against ROOT's 29. Rust stops at EDM `1.492603e-05`; ROOT reaches `1.072811e-25`.
  - Hesse uses 36 calls against ROOT's 39. The Rust count includes 22 Migrad calls and 14 Hesse calls.
  - The runner uses analytical gradients in Rust. The remaining gap includes different seed work and stopping EDM; the old counter reset does not explain it.
  - Correctness metrics pass tolerances, including `quadratic3_fixx_migrad` parameter max abs diff `3.886512e-04` and `quadratic3_fixx_hesse` covariance max abs diff `1.006964e-09` in the current differential report.
- Disposition: explicit per-workload waiver; no global NFCN threshold change.
- Next step:
  - Retain the existing manifest waivers. Isolate the remaining seed and stopping differences before changing the numerical path.

## D-003: Residual Rosenbrock covariance differences

- Workloads and measured maximum absolute covariance differences:
  - `rosenbrock2_migrad`: `1.432156e-08`
  - `rosenbrock2_minimize`: `1.432156e-08`
  - `rosenbrock2_migrad_strategy2`: `7.781219e-09`
- Classification: `numerical/residual`
- Severity: `low`; all covariance differences remain below the configured limit of 0.5.
- Evidence: `reports/verification/diff_results.csv` and the stored outputs in `reports/verification/raw/ref/` and `reports/verification/raw/rust/`.
- Details: the strategy-2 seed Hesse, triangular covariance-update sum, and line-search corrections closed the former 0.13–0.18 gap.
- Disposition: retain this entry for the measured residuals; their precise rounding source remains unisolated.

## D-004: No-g2 warning reflects numerical Hesse work, not a ROOT undercount

- Workload: `quadratic2_no_g2_migrad`.
- Classification: `algorithmic/hessian-path`.
- Severity: `low` (numerical metrics pass; call-count warning remains).
- Evidence programs: `tools/root_evidence/nfcn_accounting.cpp` and `tools/root_evidence/nfcn_accounting.rs` reuse the runners' objectives.
- Stored output:
  - [nfcn_accounting_quadratic2_no_g2_migrad.root.json](raw/trace/nfcn_accounting_quadratic2_no_g2_migrad.root.json)
  - [nfcn_accounting_quadratic2_no_g2_migrad.rust.json](raw/trace/nfcn_accounting_quadratic2_no_g2_migrad.rust.json)
- Measurements:

  | Runner | Reported `nfcn` | User-side function calls | Analytical Hessian calls | State call totals |
  |---|---:|---:|---:|---|
  | ROOT `v6-36-08` | 2 | 2 | 2 | 1, 2, 2 |
  | Rust | 12 | 12 | 1 | 2, 12 |

- ROOT's count is accurate. Rust's ten extra evaluations occur during numerical Hesse verification after iteration reaches two calls.
- ROOT `MnHesse.cxx:87-95,100-168` uses the supplied analytical Hessian without objective evaluations during verification.
  Rust `src/hesse/calculator.rs::calculate` uses finite differences, even when the objective supplies a Hessian.
- Rust's 12 calls comprise 1 seed, 1 Migrad iteration, 1 Hesse centre, 8 diagonal, 0 refinement, and 1 mixed-derivative evaluation.
  ROOT makes 1 seed and 1 iteration evaluation, with no objective calls during analytical Hesse.
  Rust adds the centre evaluation to match ROOT's numerical Hesse.
  The old Rust count omitted its seed and centre evaluations.
- The relative difference is `10/12 = 0.833333`, above the unchanged threshold of `0.8`.
- Disposition: leave the warning open. Using the supplied analytical Hessian in Hesse changes the numerical path.
  Matching ROOT's analytical path is a separate change. No waiver is added.

## D-005: Hesse gradient refinement and seed gradient steps differ from ROOT

- Follow-up: aligning the steps with ROOT's would tighten the Rosenbrock and Hahn1 call-count bounds to zero.
- Rosenbrock strategy 1 uses 144 calls against ROOT's 140. Both enter Hesse verification at 130.
  Rust then makes 1 centre, 4 diagonal, 8 refinement and 1 mixed-derivative evaluation.
  ROOT makes 1 centre, 4 diagonal, 4 refinement and 1 mixed-derivative evaluation.
  The four extra refinement calls explain the entire current gap.
- The Rosenbrock traces first differ at call 4, during seed generation on both sides:
  ROOT evaluates `(0.00069049922586661124, 0)`; Rust evaluates `(0.00069049941105437588, 0)`.
  This change leaves the numerical refinement and seed-step rules unchanged.
- Rust's refinement starts from zero and needs at least two cycles per parameter.
  ROOT starts from the diagonal sweep's gradient and takes one cycle on these fits
  (`HessianGradientCalculator.cxx:68,112-127`; Rust `src/hesse/gradient.rs`).
- Hahn1 uses 470 calls against ROOT's 456. The stored counter harnesses measure these phases:

  | Phase | ROOT calls | Rust calls |
  |---|---:|---:|
  | Before Hesse | 294 | 294 |
  | Centre evaluation | 1 | 1 |
  | Diagonal evaluations | 14 | 14 |
  | Gradient refinement in that pass | 14 | 28 |
  | Mixed derivatives | 21 | 21 |
  | Continuation after Hesse | 112 | 112 |
  | Final total | 456 | 470 |

  Each side makes one Hesse pass.
  ROOT's pass spans calls 295–344, with refinement at 310–323.
  Rust's pass spans calls 295–358, with refinement at 310–337.
  Refinement explains all fourteen extra calls; zero calls remain unexplained.
- Evidence programs: `tools/root_evidence/nfcn_accounting.cpp` and `.rs` identify the Rosenbrock seed and refinement differences.
  `tools/root_evidence/hahn1_phase_counts.cpp` and `.rs` count Hahn1 calls using the probe's objective and the parity test's objective.
  Their shared phase counters classify stack frames and the number of coordinates changed from the Hesse centre.
- Stored output:
  - [nfcn_accounting_rosenbrock2_maxfcn140.root.json](raw/trace/nfcn_accounting_rosenbrock2_maxfcn140.root.json)
  - [nfcn_accounting_rosenbrock2_maxfcn140.rust.json](raw/trace/nfcn_accounting_rosenbrock2_maxfcn140.rust.json)
  - [hahn1_phase_counts_hahn1_start2.root.json](raw/trace/hahn1_phase_counts_hahn1_start2.root.json)
  - [hahn1_phase_counts_hahn1_start2.rust.json](raw/trace/hahn1_phase_counts_hahn1_start2.rust.json)
- Disposition: retain bounds of four calls for Rosenbrock and fourteen for Hahn1 until the steps are aligned.
  These totals include Rust's centre evaluation, matching ROOT's numerical Hesse.
  No manifest tolerance or waiver changes.

## D-006: Minos row counts before Minos

- Workloads: `quadratic2_minos_p0` and `quadratic2_minos_p1`.
- Classification: `algorithmic/hessian-path`.
- Both runners store the minimum's count after Migrad and standalone Hesse, before Minos runs.
  ROOT reports 28 Migrad calls plus 10 standalone Hesse calls, totalling 38.
  Rust reports 32 Migrad calls plus 14 standalone Hesse calls, totalling 46.
- Calls by phase:

  | Phase | ROOT calls | Rust calls |
  |---|---:|---:|
  | Seed | 7 | 7 |
  | Migrad iteration | 11 | 11 |
  | Verification Hesse centre | 1 | 1 |
  | Verification Hesse diagonal | 4 | 4 |
  | Verification Hesse refinement | 4 | 8 |
  | Verification Hesse mixed derivatives | 1 | 1 |
  | Standalone Hesse centre | 1 | 1 |
  | Standalone Hesse diagonal | 4 | 4 |
  | Standalone Hesse refinement | 4 | 8 |
  | Standalone Hesse mixed derivatives | 1 | 1 |
  | Final total before Minos | 38 | 46 |

  Each side makes two Hesse passes: one during Migrad verification, then one standalone pass.
  Each Rust pass adds four refinement calls, explaining all eight extra calls.
  This is the Hesse gradient refinement difference described above; no calls remain unexplained.
- Evidence programs: `tools/root_evidence/nfcn_accounting.cpp` and `.rs` repeat the runners' shared setup and stop before Minos.
- Stored output:
  - [nfcn_accounting_quadratic2_minos_before_minos.root.json](raw/trace/nfcn_accounting_quadratic2_minos_before_minos.root.json)
  - [nfcn_accounting_quadratic2_minos_before_minos.rust.json](raw/trace/nfcn_accounting_quadratic2_minos_before_minos.rust.json)
- Disposition: retain the measured counts until the refinement steps are aligned.

## Resolved since previous snapshot

- Previous `NoG2CallsWhenFCHasNoG2` P0 gap resolved by adding optional FCN `has_g2()/g2()` and `has_hessian()/hessian()` contract and porting the regression as `tests/root_reference_minuit2.rs::root_no_g2_calls_when_fcn_has_no_g2`.
- Previous `quadratic2_simplex` correctness fail resolved after aligning Simplex EDM tolerance semantics with ROOT (`minedm = tolerance * Up`).
- Previous `rosenbrock2_minimize` correctness fail resolved after aligning `MnMinimize` flow with ROOT combined minimizer logic and propagating covariance from final minimum state.

## Reference trace evidence

- The stored ROOT trace of the strategy-2 Rosenbrock run is generated by `tools/ref_runner_cpp/main.cpp` using `MnTraceObject`, attached through `Builder().SetTraceObject`.
  The program's header gives build and run commands for `reports/verification/raw/trace/rosenbrock2_migrad_strategy2.root.jsonl`.
  It uses the reference workload's objective, start, steps, and strategy.
  Iteration 1 gives `dcovar = 0.16098749348302552`; successful seed Hesse gives exactly zero.
- The stored direct ROOT v6-36-08 line-search evidence is generated by `tools/root_evidence/linesearch_open_step.cpp`.
  Its header gives build and run commands for `reports/verification/raw/trace/linesearch_open_step.root.json`.
  It calls `MnLineSearch` on the quartic in `tests/linesearch_open_step.rs`.
  The stored start, direction, and directional derivative reproduce lambda `10.39095385821868` and fval `-2191.6554284992781`.
  Both records use the existing reference Minuit2 library.
