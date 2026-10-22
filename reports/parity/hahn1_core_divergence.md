# Hahn1 core divergence diagnosis and fix

## Summary

Hahn1 failed because minuit2-rs stopped after a nominal low-EDM MIGRAD pass while the covariance was still uncertain. ROOT Minuit2 verifies that case with MnHesse and, if Hesse raises EDM above tolerance, continues MIGRAD from the Hesse state with an enlarged call budget. Implementing that ROOT-cited verified-continuation path, together with ROOT's per-iteration EDM/update ordering and Davidon negative-curvature update behavior, makes minuit2-rs reach the certified Hahn1 basin from NIST Start 2.

## Why `valid=False`

Pre-fix baseline reproduction:

```text
. .venv-maturin/bin/activate && python scripts/nist_hard_baseline.py --dataset Hahn1
Hahn1 | OK (valid=True, fval=1.53244, params=True, nfcn=581) | FAIL (valid=True, fval=44997.1, params=False, nfcn=580) | FAIL (valid=False, fval=51623.2, params=False, nfcn=72) | FAIL (valid=False, fval=46949.9, params=False, nfcn=260)
```

Direct `fmin` inspection showed the invalid flag was `above_max_edm`, not bad parameters or call limit:

```text
start2 strategy 1: valid=False fval=51623.23745281273 nfcn=72 edm=1.0088355348927458e+16 is_above_max_edm=True has_reached_call_limit=False has_valid_parameters=True
start2 strategy 2: valid=False fval=46949.88221182479 nfcn=260 edm=2.395921281972186e+17 is_above_max_edm=True has_reached_call_limit=False has_valid_parameters=True
```

## Trajectory comparison and divergence point

Using a local ROOT Minuit2 trace runner built against `third_party/root_ref` and minuit2-rs compiled with the `trace` feature, both with NIST Start 2 and user errors of `0.1`, accepted-step traces differ numerically from the first accepted line-search step:

| accepted iteration | ROOT Minuit2 nfcn | ROOT fval | ROOT edm | minuit2-rs nfcn | minuit2-rs fval | minuit2-rs edm |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 76 | 1412.2155975159774 | 1840139.4174936966 | 97 | 946.5561824770050 | 2318.203029072101 |

Trace artifacts:

- `reports/verification/raw/trace/hahn1.root.jsonl`
- `reports/verification/raw/trace/hahn1.rs.jsonl`

The decisive missing ROOT behavior occurred later at the first nominal convergence point. minuit2-rs reached fval=3.232059 with EDM=1.3596e-4 and `Dcovar=0.05159`; pre-fix it accepted that low-EDM/local state. ROOT's `VariableMetricBuilder` checks `strategy == 1 && min.Error().Dcovar() > 0.05`, runs MnHesse, and continues MIGRAD if Hesse recomputes EDM above tolerance and above the machine-accuracy floor. With this behavior, minuit2-rs continues through fval=1.56445, Hesse-verifies again, and reaches fval=1.532439 with certified parameters.

## ROOT-cited mechanisms implemented

- **Hesse-verified continuation**: ROOT runs MnHesse for strategy >= 2, or strategy 1 with `Dcovar() > 0.05`; if Hesse EDM remains above tolerance and above machine accuracy, it continues minimization with 1.3x call budget. Citation: `third_party/root_ref/math/minuit2/src/VariableMetricBuilder.cxx:112-151`.
- **EDM/update ordering**: ROOT calls `Estimator().Estimate(g, s0.Error())` before `ErrorUpdator().Update(s0, p, g)` and only corrects the loop-local EDM by `(1 + 3 * e.Dcovar())` after storing the state. Citation: `third_party/root_ref/math/minuit2/src/VariableMetricBuilder.cxx:291-312`.
- **Davidon negative-curvature update**: ROOT only returns without update for `delgam == 0` or `gvg <= 0`; for `delgam < 0`, it warns but still applies the update. Citation: `third_party/root_ref/math/minuit2/src/DavidonErrorUpdator.cxx:42-78`.

## Regression test

Added `tests/hahn1_core_parity.rs`, which fits unscaled Hahn1 from NIST Start 2 with strategy 1 and asserts validity, certified RSS basin, and certified parameters at 1e-2 relative tolerance.

Pre-change failing command:

```text
cargo test --test hahn1_core_parity -- --nocapture
```

Pre-change failure mode was the baseline failure (`valid=False`, fval ≈ 51623, above-max EDM, and parameters outside certified tolerance).

Post-fix focused result:

```text
cargo test --test hahn1_core_parity -- --nocapture
running 1 test
test hahn1_core_start2_reaches_certified_solution ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.01s
```

Post-fix hard baseline:

```text
. .venv-maturin/bin/activate && python scripts/nist_hard_baseline.py
Hahn1 | OK (valid=True, fval=1.53244, params=True, nfcn=581) | FAIL (valid=True, fval=44997.1, params=False, nfcn=580) | OK (valid=True, fval=1.53244, params=True, nfcn=721) | FAIL (valid=False, fval=12076.9, params=False, nfcn=514)
```
