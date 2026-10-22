# Function Parity Gaps

Upstream repo: `root-project/root`
Upstream subdir: `math/minuit2`
Upstream ref: `v6-36-08`
Upstream commit: `a8ca1b23e38d7dbe0ff24027894ca0f2ad65f1bd`

## Summary

- Total upstream symbols in scope: **415**
- `implemented`: **217**
- `missing`: **0**
- `needs-review`: **170**
- `intentionally-skipped`: **28**

## Top Files by `missing` Symbols

- None

## Top Files by `needs-review` Symbols

- `inc/Minuit2/MnUserParameterState.h`: 27
- `inc/Minuit2/MnApplication.h`: 21
- `inc/Minuit2/MnPrint.h`: 20
- `inc/Minuit2/MnUserTransformation.h`: 20
- `inc/Minuit2/MnStrategy.h`: 10
- `inc/Minuit2/MinimumBuilder.h`: 9
- `inc/Minuit2/FunctionMinimum.h`: 7
- `src/MnPrint.cxx`: 5
- `inc/Minuit2/CombinedMinimumBuilder.h`: 4
- `inc/Minuit2/FCNBase.h`: 4
- `inc/Minuit2/MnTraceObject.h`: 4
- `inc/Minuit2/GradientCalculator.h`: 2
- `inc/Minuit2/MinuitParameter.h`: 2
- `inc/Minuit2/MnParameterScan.h`: 2
- `inc/Minuit2/BasicFunctionGradient.h`: 1
- `inc/Minuit2/BasicFunctionMinimum.h`: 1
- `inc/Minuit2/BasicMinimumError.h;src/BasicMinimumError.cxx`: 1
- `inc/Minuit2/BasicMinimumParameters.h`: 1
- `inc/Minuit2/BasicMinimumSeed.h`: 1
- `inc/Minuit2/BasicMinimumState.h`: 1

## Notes

- Symbol extraction is heuristic (regex-based), not a full C++ parser.
- `intentionally-skipped` currently captures constructor/destructor/operator-style symbols that map to Rust idioms.
- `needs-review` includes architectural refactors where strict 1:1 symbol naming is not expected.
- Use `reports/parity/functions.csv` as the source of truth for triage and manual confirmation.
