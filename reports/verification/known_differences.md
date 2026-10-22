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
  - `math/minuit2/src/MnSeedGenerator.cxx`: numerical seed generation first evaluates the starting point with `MnFcnCaller{fcn}(x)` and then calls `gc(pa)`.
  - `math/minuit2/src/Numerical2PGradientCalculator.cxx`: `operator()(par)` constructs an `InitialGradientCalculator` and then central-difference probes each variable parameter via `mfcnCaller(x)` at `xtf + step` and `xtf - step` for up to `GradientNCycles()` cycles.
  - `math/minuit2/src/InitialGradientCalculator.cxx`: the initial rough gradient itself is heuristic and does not call the FCN.
- Details:
  - The extra C++ evaluations are seed-phase numerical-gradient probes counted by ROOT's `MnFcn`; they are not Hesse sweeps or covariance-squeeze work.
  - The fixed parameter is not differentiated as a variable slot in either implementation; the divergence appears before Migrad iteration 0 because ROOT counts the seed central-difference probes while minuit2-rs reaches the same fixed-parameter minimum without reproducing that extra seed work.
  - Correctness metrics pass tolerances, including `quadratic3_fixx_migrad` parameter max abs diff `3.886512e-04` and `quadratic3_fixx_hesse` covariance max abs diff `9.940541e-10` in the current differential report.
- Disposition: explicit per-workload waiver; no global NFCN threshold change.
- Next step:
  - If future fixed-parameter workloads show correctness degradation, revisit route (a) and reproduce the ROOT seed probe accounting in code.

## D-003: Rosenbrock covariance differs but remains within configured tolerance

- Workloads:
  - `rosenbrock2_migrad`
  - `rosenbrock2_minimize`
- Classification: `algorithmic/path-efficiency`
- Severity: `low` (correctness gate passes with explicit covariance tolerance)
- Expected: `yes`
- Details:
  - Covariance presence now matches and is enforced in the differential gate.
  - Covariance element max abs difference is currently ~0.39 (tolerance: 0.5), consistent with stopping at slightly different nearby points.
- Disposition: `document`
- Next step:
  - Add stricter Rosenbrock convergence workloads (lower tolerance and/or larger call budget) if tighter covariance parity is required.

## Resolved since previous snapshot

- Previous `NoG2CallsWhenFCHasNoG2` P0 gap resolved by adding optional FCN `has_g2()/g2()` and `has_hessian()/hessian()` contract and porting the regression as `tests/root_reference_minuit2.rs::root_no_g2_calls_when_fcn_has_no_g2`.
- Previous `quadratic2_simplex` correctness fail resolved after aligning Simplex EDM tolerance semantics with ROOT (`minedm = tolerance * Up`).
- Previous `rosenbrock2_minimize` correctness fail resolved after aligning `MnMinimize` flow with ROOT combined minimizer logic and propagating covariance from final minimum state.
