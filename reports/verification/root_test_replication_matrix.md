# ROOT Minuit2 Test Replication Matrix

Source of truth for C++ reference tests:
- Repository: `https://github.com/root-project/root`
- Path: `math/minuit2/test`
- Baseline tag: `v6-36-08`

This matrix defines which ROOT tests we can replicate directly in pure Rust (`minuit2-rs`) and which require adaptation because they depend on ROOT fitting/histogram infrastructure outside Minuit2 core.

## Test Inventory

Top-level tests discovered in `math/minuit2/test`:
- `testMinuit2.cxx` (gtest regression tests)
- `testMinimizer.cxx`
- `testCovariance.cxx`
- `testADMinim.cxx` (conditional on autodiff/clad setup)
- `testNdimFit.cxx` (ROOT Fit framework)
- `testUnbinGausFit.cxx` (ROOT Fit framework)
- `testUserFunc.cxx` (ROOT hist/TF1/UI path)
- `MnTutorial/*` (Quad1F/4F/8F/12F)
- `MnSim/*` demos and long tests

## Replication Plan

| ROOT test | Core Minuit2 only? | Rust replication strategy | Priority |
|---|---|---|---|
| `testMinuit2.cxx` | Yes | Replicate as regression tests in `tests/root_reference_minuit2.rs`; include Hessian external-indexing and no-`G2` call behavior checks | P0 |
| `testCovariance.cxx` | Mostly yes | Replicate covariance transform round-trip/consistency checks into `tests/root_reference_covariance.rs` | P0 |
| `testMinimizer.cxx` | Mostly yes | Extract deterministic objective families (Rosenbrock, ChebyQuad, TrigoFletcher) into fixture-based differential tests | P0 |
| `MnTutorial/Quad*` | Yes | Replicate as dimensional scaling correctness/robustness tests | P1 |
| `MnSim/ParallelTest.cxx` | Mostly yes | Replicate as controlled performance/consistency test with deterministic data | P1 |
| `MnSim/PaulTest*.cxx`, `ReneTest.cxx` | Mixed | Replicate algorithmic core portions; replace ROOT I/O with checked-in fixture data | P2 |
| `testADMinim.cxx` | Mixed | Replicate only if Rust-side AD path is in scope; otherwise track as out-of-scope with rationale | P2 |
| `testNdimFit.cxx` | No (ROOT Fit layer) | Recreate equivalent statistical objective directly in Rust without ROOT `Fit::Fitter` | P2 |
| `testUnbinGausFit.cxx` | No (ROOT Fit layer) | Recreate likelihood workload directly using Rust FCN interface | P2 |
| `testUserFunc.cxx` | No (hist/TF1/UI) | Treat as integration/UI test not applicable to pure Minuit2 core | P3 |

## Current Replication Status (Implemented)

| ROOT source | Rust test | Status | Notes |
|---|---|---|---|
| `testMinuit2.cxx` `HessianExternalIndexing_Numeric` | `tests/root_reference_minuit2.rs::root_hessian_external_indexing_numeric` | Implemented | Passes |
| `testMinuit2.cxx` `HessianExternalIndexing_Analytic` | `tests/root_reference_minuit2.rs::root_hessian_external_indexing_with_gradient` | Implemented | Passes |
| `testCovariance.cxx` upper/lower/double/unbounded covariance mapping checks | `tests/root_reference_covariance.rs::{root_covariance_upper,root_covariance_lower,root_covariance_double,root_covariance_unbounded}` | Implemented | Passes |
| `testMinuit2.cxx` `NoG2CallsWhenFCHasNoG2` | `tests/root_reference_minuit2.rs::root_no_g2_calls_when_fcn_has_no_g2` | Implemented | Passes |

## Rules for Replication

1. Preserve objective definitions and parameter initialization exactly where feasible.
2. Remove ROOT framework dependencies only at the harness layer, never by changing the optimization math.
3. Keep deterministic seeds/datasets and commit them under `verification/workloads/`.
4. For each replicated test, store:
   - source ROOT file reference,
   - mapped Rust test name,
   - tolerance contract,
   - differential result status.

## Deliverables

- `tests/root_reference_minuit2.rs`
- `tests/root_reference_covariance.rs`
- `verification/workloads/root_reference_*.json`
- `reports/verification/root_test_port_status.md` (coverage of ROOT tests replicated)

## Acceptance for "ROOT tests replicated"

Minimum bar:
- All P0 tests implemented and passing against both reference runner and Rust runner.
- P1/P2 tests either implemented or explicitly marked with rationale and ETA.
- No silent omissions.
