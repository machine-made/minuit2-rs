# Executed Surface Mapping

Join of reference executed C++ functions with traceability matrix mappings.

## Summary

- Executed C++ functions: **618**
- Mapped to implemented Rust symbols: **145**
- Unmapped executed functions: **473**
- Unmapped priority split: P0=0, P1=48, P2=425
- Gate (`P0 == 0 and P1 == 0`): **FAIL**
- Coverage workloads used: **12**

## Artifacts

- `reports/verification/executed_surface_mapping.md`
- `reports/verification/executed_surface_gaps.csv`
- `reports/verification/executed_surface_manifest.json`

## Top Gap Files

- `inc/Minuit2/MnPrint.h`: 220
- `inc/Minuit2/MnMatrix.h`: 75
- `inc/ROOT/span.hxx`: 22
- `src/MnMatrix.cxx`: 14
- `src/MnPrint.cxx`: 9
- `inc/Minuit2/FunctionGradient.h`: 8
- `inc/Minuit2/MinimumSeed.h`: 6
- `src/MnHesse.cxx`: 6
- `inc/Minuit2/MinimumParameters.h`: 5
- `src/MPIProcess.h`: 5
- `src/MnApplication.cxx`: 5
- `inc/Minuit2/FunctionMinimum.h`: 4
- `inc/Minuit2/MinimumState.h`: 4
- `inc/Minuit2/MnUserParameterState.h`: 4
- `src/MPIProcess.cxx`: 4

## Top P0/P1 Gaps

| Priority | Upstream file | Symbol | Mapping status | Call count |
|---|---|---|---|---|
| P1 | `inc/Minuit2/FunctionGradient.h` | `G2` | `missing` | 206 |
| P1 | `inc/Minuit2/FunctionGradient.h` | `Grad` | `missing` | 992 |
| P1 | `inc/Minuit2/FunctionGradient.h` | `Gstep` | `missing` | 169 |
| P1 | `inc/Minuit2/FunctionGradient.h` | `Vec` | `missing` | 600 |
| P1 | `inc/Minuit2/FunctionMinimum.h` | `Error` | `missing` | 29 |
| P1 | `inc/Minuit2/FunctionMinimum.h` | `Parameters` | `missing` | 8 |
| P1 | `inc/Minuit2/HessianGradientCalculator.h` | `Fcn` | `missing` | 15 |
| P1 | `inc/Minuit2/HessianGradientCalculator.h` | `Strategy` | `missing` | 30 |
| P1 | `inc/Minuit2/MinimumError.h` | `InvHessian` | `missing` | 672 |
| P1 | `inc/Minuit2/MinimumParameters.h` | `Dirin` | `missing` | 2 |
| P1 | `inc/Minuit2/MinimumParameters.h` | `Vec` | `missing` | 1474 |
| P1 | `inc/Minuit2/MinimumSeed.h` | `Error` | `missing` | 59 |
| P1 | `inc/Minuit2/MinimumSeed.h` | `Gradient` | `missing` | 60 |
| P1 | `inc/Minuit2/MinimumSeed.h` | `Parameters` | `missing` | 146 |
| P1 | `inc/Minuit2/MinimumSeed.h` | `Precision` | `missing` | 28 |
| P1 | `inc/Minuit2/MinimumSeed.h` | `State` | `missing` | 381 |
| P1 | `inc/Minuit2/MinimumState.h` | `Gradient` | `missing` | 735 |
| P1 | `inc/Minuit2/MinimumState.h` | `Vec` | `missing` | 366 |
| P1 | `inc/Minuit2/MinuitParameter.h` | `GetName` | `missing` | 484 |
| P1 | `inc/Minuit2/MnFcn.h` | `Fcn` | `missing` | 853 |
| P1 | `inc/Minuit2/MnFcn.h` | `Trafo` | `missing` | 1334 |
| P1 | `inc/Minuit2/MnUserParameterState.h` | `IntCovariance` | `missing` | 32 |
| P1 | `inc/Minuit2/MnUserParameterState.h` | `IntParameters` | `missing` | 38 |
| P1 | `inc/Minuit2/MnUserParameterState.h` | `Parameters` | `missing` | 45 |
| P1 | `inc/Minuit2/MnUserParameterState.h` | `Trafo` | `missing` | 128 |
| P1 | `inc/Minuit2/MnUserParameters.h` | `Trafo` | `missing` | 417 |
| P1 | `inc/Minuit2/MnUserTransformation.h` | `InitialParValues` | `missing` | 488 |
| P1 | `inc/Minuit2/MnUserTransformation.h` | `Parameters` | `missing` | 94 |
| P1 | `inc/Minuit2/MnUserTransformation.h` | `Precision` | `missing` | 553 |
| P1 | `inc/Minuit2/VariableMetricBuilder.h` | `ErrorUpdator` | `missing` | 81 |
| P1 | `inc/Minuit2/VariableMetricBuilder.h` | `Estimator` | `missing` | 81 |
| P1 | `src/HessianGradientCalculator.cxx` | `DeltaGradient` | `missing` | 15 |
| P1 | `src/HessianGradientCalculator.cxx` | `Precision` | `missing` | 165 |
| P1 | `src/MnHesse.cxx` | `ComputeAnalytical` | `missing` | 1 |
| P1 | `src/MnHesse.cxx` | `ComputeNumerical` | `missing` | 15 |
| P1 | `src/MnScan.cxx` | `Scan` | `missing` | 1 |
| P1 | `src/MnUserParameterState.cxx` | `MinuitParameters` | `missing` | 12 |
| P1 | `src/MnUserParameterState.cxx` | `Name` | `missing` | 16 |
| P1 | `src/MnUserParameterState.cxx` | `Parameter` | `missing` | 380 |
| P1 | `src/MnUserParameterState.cxx` | `Precision` | `missing` | 44 |

