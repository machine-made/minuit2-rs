# ROOT Test Replication Status (`v6-36-08`)

Reference:
- Repository: `https://github.com/root-project/root`
- Subtree: `math/minuit2`
- Tag: `v6-36-08`

## P0 Regression Replications

| ROOT test intent | Rust evidence | Status |
|---|---|---|
| Hessian external indexing (numeric) | `tests/root_reference_minuit2.rs::root_hessian_external_indexing_numeric` | Passing |
| Hessian external indexing (analytical gradient path) | `tests/root_reference_minuit2.rs::root_hessian_external_indexing_with_gradient` | Passing |
| No G2 calls when FCN reports `has_g2() == false` | `tests/root_reference_minuit2.rs::root_no_g2_calls_when_fcn_has_no_g2` | Passing |
| Covariance transform consistency (unbounded/upper/lower/double) | `tests/root_reference_covariance.rs::{root_covariance_unbounded,root_covariance_upper,root_covariance_lower,root_covariance_double}` | Passing |

## Last Validation Commands

```bash
cargo test --test root_reference_minuit2
cargo test --test root_reference_covariance
```

## Differential Harness Snapshot

- Summary: `reports/verification/diff_summary.md`
- CSV: `reports/verification/diff_results.csv`
- Snapshot counts: pass=10, warn=2, fail=0
- Warning-only cases: NFCN divergence on `quadratic3_fixx_migrad` and `quadratic3_fixx_hesse`
- Failing cases: none
