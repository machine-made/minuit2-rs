# Similarity Audit

This is a mechanical provenance triage report, not a legal conclusion.

## Inputs

- ROOT reference: `third_party/root_ref/math/minuit2`
- Parity mapping: `reports/parity/functions.csv`
- Rust implementation scope: `src/**/*.rs`
- ROOT scope: `third_party/root_ref/math/minuit2/**/*.{h,cxx,cpp,cc,c}`

## Method

- Strip comments and string literals before code-token comparison.
- Compare retained lexical tokens and generic tokens, where identifiers become `ID` and numbers become `NUM`.
- Compare exact code token shingles, exact generic code shingles, exact comment word shingles, and string literals.
- Mark ROOT/Rust pairs from `reports/parity/functions.csv` as mapped; also scan all unmapped pairs for stronger accidental matches.

Risk labels are triage labels:

- `high`: copied comments/strings or exact retained code-token shingles need immediate human review.
- `medium`: mapped pair with enough lexical/structural overlap to inspect manually.
- `low-review`: mapped pair with weak signals worth keeping in the inventory.
- `low`: no mechanical similarity signal beyond normal algorithm/API vocabulary.

## Summary

- Pairs scanned: **10438**
- High risk: **11**
- Medium risk: **58**
- Low-review: **19**
- Low: **10350**
- Mapped high risk: **1**
- Mapped medium risk: **58**

## High-Risk Findings

| Rust file | ROOT file | Signals |
|---|---|---|
| `src/bin/ref_compare_runner.rs` | `test/testMinuit2.cxx` | comments=0, code=22, strings=0 |
| `src/covariance_squeeze.rs` | `test/testMinimizer.cxx` | comments=0, code=12, strings=0 |
| `src/covariance_squeeze.rs` | `test/testADMinim.cxx` | comments=0, code=12, strings=0 |
| `src/covariance_squeeze.rs` | `test/testMinuit2.cxx` | comments=0, code=6, strings=0 |
| `src/bin/ref_compare_runner.rs` | `test/testUserFunc.cxx` | comments=0, code=4, strings=0 |
| `src/bin/ref_compare_runner.rs` | `test/MnTutorial/Quad4F.h` | comments=0, code=4, strings=0 |
| `src/bin/ref_compare_runner.rs` | `test/MnTutorial/Quad8F.h` | comments=0, code=4, strings=0 |
| `src/bin/ref_compare_runner.rs` | `test/MnTutorial/Quad12F.h` | comments=0, code=4, strings=0 |
| `src/user_parameters.rs` | `src/MnUserTransformation.cxx` | comments=3, code=0, strings=0 |
| `src/gradient/analytical.rs` | `test/testMinuit2.cxx` | comments=0, code=3, strings=0 |
| `src/parameter.rs` | `test/testMinimizer.cxx` | comments=0, code=3, strings=0 |

## Medium-Risk Triage Queue

| Rust file | ROOT file | token cosine | generic cosine | generic shingles |
|---|---|---:|---:|---:|
| `src/minimum/mod.rs` | `inc/Minuit2/FunctionMinimum.h` | 0.8653 | 0.9565 | 1 |
| `src/linesearch.rs` | `src/MnLineSearch.cxx` | 0.8621 | 0.9753 | 12 |
| `src/migrad/builder.rs` | `src/BFGSErrorUpdator.cxx` | 0.8254 | 0.9566 | 6 |
| `src/user_covariance.rs` | `inc/Minuit2/MnUserCovariance.h` | 0.8079 | 0.926 | 0 |
| `src/user_parameter_state.rs` | `inc/Minuit2/FunctionMinimum.h` | 0.7736 | 0.9345 | 0 |
| `src/parameter.rs` | `inc/Minuit2/MinuitParameter.h` | 0.7561 | 0.9338 | 0 |
| `src/minimum/state.rs` | `inc/Minuit2/MinimumState.h` | 0.7313 | 0.9129 | 0 |
| `src/parabola.rs` | `inc/Minuit2/MnParabola.h` | 0.7218 | 0.8814 | 0 |
| `src/user_parameters.rs` | `inc/Minuit2/MnUserTransformation.h` | 0.7132 | 0.9146 | 0 |
| `src/hesse/mod.rs` | `inc/Minuit2/MnHesse.h` | 0.7031 | 0.9274 | 0 |
| `src/user_transformation.rs` | `inc/Minuit2/MnUserTransformation.h` | 0.6931 | 0.9077 | 0 |
| `src/migrad/builder.rs` | `inc/Minuit2/VariableMetricBuilder.h` | 0.6898 | 0.9454 | 0 |
| `src/migrad/builder.rs` | `inc/Minuit2/MinimumBuilder.h` | 0.6825 | 0.9216 | 0 |
| `src/minos/minos_error.rs` | `inc/Minuit2/MinosError.h` | 0.675 | 0.8596 | 0 |
| `src/minimum/seed.rs` | `inc/Minuit2/MinimumSeed.h` | 0.6698 | 0.8704 | 0 |
| `src/minimize/mod.rs` | `inc/Minuit2/MnMinimize.h` | 0.6576 | 0.9385 | 15 |
| `src/simplex/parameters.rs` | `inc/Minuit2/SimplexParameters.h` | 0.655 | 0.8903 | 0 |
| `src/minimum/error.rs` | `inc/Minuit2/MinimumError.h` | 0.6521 | 0.8885 | 0 |
| `src/precision.rs` | `inc/Minuit2/MnMachinePrecision.h` | 0.6459 | 0.872 | 0 |
| `src/user_parameters.rs` | `inc/Minuit2/MnUserParameters.h` | 0.6438 | 0.8866 | 0 |
| `src/minimum/parameters.rs` | `inc/Minuit2/MinimumParameters.h` | 0.6403 | 0.9191 | 0 |
| `src/minimum/error.rs` | `inc/Minuit2/FunctionMinimum.h` | 0.6371 | 0.8756 | 0 |
| `src/minos/cross.rs` | `inc/Minuit2/MnCross.h` | 0.6298 | 0.8683 | 0 |
| `src/contours/mod.rs` | `inc/Minuit2/MnContours.h` | 0.622 | 0.924 | 0 |
| `src/simplex/mod.rs` | `inc/Minuit2/MnSimplex.h` | 0.6182 | 0.9204 | 0 |

## Top Mechanical Similarity Rows

| Risk | Rust file | ROOT file | mapped | token cosine | generic shingles | comments | strings |
|---|---|---|---|---:|---:|---:|---:|
| high | `src/bin/ref_compare_runner.rs` | `test/testMinuit2.cxx` | no | 0.7567 | 22 | 0 | 0 |
| high | `src/covariance_squeeze.rs` | `test/testMinimizer.cxx` | no | 0.8177 | 12 | 0 | 0 |
| high | `src/covariance_squeeze.rs` | `test/testADMinim.cxx` | no | 0.8169 | 12 | 0 | 0 |
| high | `src/covariance_squeeze.rs` | `test/testMinuit2.cxx` | no | 0.8198 | 6 | 0 | 0 |
| high | `src/bin/ref_compare_runner.rs` | `test/testUserFunc.cxx` | no | 0.8693 | 4 | 0 | 0 |
| high | `src/bin/ref_compare_runner.rs` | `test/MnTutorial/Quad4F.h` | no | 0.5888 | 4 | 0 | 0 |
| high | `src/bin/ref_compare_runner.rs` | `test/MnTutorial/Quad8F.h` | no | 0.5038 | 5 | 0 | 0 |
| high | `src/bin/ref_compare_runner.rs` | `test/MnTutorial/Quad12F.h` | no | 0.4774 | 5 | 0 | 0 |
| high | `src/user_parameters.rs` | `src/MnUserTransformation.cxx` | yes | 0.827 | 0 | 3 | 0 |
| high | `src/gradient/analytical.rs` | `test/testMinuit2.cxx` | no | 0.8141 | 4 | 0 | 0 |
| high | `src/parameter.rs` | `test/testMinimizer.cxx` | no | 0.6108 | 5 | 0 | 0 |
| medium | `src/minimum/mod.rs` | `inc/Minuit2/FunctionMinimum.h` | yes | 0.8653 | 1 | 0 | 0 |
| medium | `src/linesearch.rs` | `src/MnLineSearch.cxx` | yes | 0.8621 | 12 | 0 | 0 |
| medium | `src/migrad/builder.rs` | `src/BFGSErrorUpdator.cxx` | yes | 0.8254 | 6 | 0 | 0 |
| medium | `src/user_covariance.rs` | `inc/Minuit2/MnUserCovariance.h` | yes | 0.8079 | 0 | 0 | 0 |
| medium | `src/user_parameter_state.rs` | `inc/Minuit2/FunctionMinimum.h` | yes | 0.7736 | 0 | 0 | 0 |
| medium | `src/parameter.rs` | `inc/Minuit2/MinuitParameter.h` | yes | 0.7561 | 0 | 0 | 0 |
| medium | `src/minimum/state.rs` | `inc/Minuit2/MinimumState.h` | yes | 0.7313 | 0 | 0 | 0 |
| medium | `src/parabola.rs` | `inc/Minuit2/MnParabola.h` | yes | 0.7218 | 0 | 0 | 0 |
| medium | `src/user_parameters.rs` | `inc/Minuit2/MnUserTransformation.h` | yes | 0.7132 | 0 | 0 | 0 |
| medium | `src/hesse/mod.rs` | `inc/Minuit2/MnHesse.h` | yes | 0.7031 | 0 | 0 | 0 |
| medium | `src/user_transformation.rs` | `inc/Minuit2/MnUserTransformation.h` | yes | 0.6931 | 0 | 0 | 0 |
| medium | `src/migrad/builder.rs` | `inc/Minuit2/VariableMetricBuilder.h` | yes | 0.6898 | 0 | 0 | 0 |
| medium | `src/migrad/builder.rs` | `inc/Minuit2/MinimumBuilder.h` | yes | 0.6825 | 0 | 0 | 0 |
| medium | `src/minos/minos_error.rs` | `inc/Minuit2/MinosError.h` | yes | 0.675 | 0 | 0 | 0 |

## Next Manual Review

1. Inspect every `high` row first, if any.
2. Inspect every `medium` row and decide `keep`, `document`, or `rewrite`.
3. For algorithm hot spots, compare against papers/manuals rather than ROOT source when rewriting.
4. Rerun this script after wording cleanup or rewrites.

Current human triage notes, if present: `reports/provenance/manual_review.md`.

## Reproduce

```bash
python3 scripts/generate_similarity_audit.py
```
