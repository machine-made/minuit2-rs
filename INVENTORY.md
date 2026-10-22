# Minuit2 File Inventory

Inventory map for Minuit2, categorized by phase and priority.

Important:
- Reference baseline for this Rust implementation is ROOT Minuit2:
  - `https://github.com/root-project/root/tree/master/math/minuit2`
  - target release baseline: `v6-36-08`
- Paths used here (`inc/Minuit2/*`, `src/*`) are compatible with ROOT `math/minuit2` and are used by parity tooling.

**Raw URL base**: `https://raw.githubusercontent.com/GooFit/Minuit2/master/`

## Legend

- **Implement** = provide Rust behavior-compatible implementation
- **Skip** = ROOT/Fumili/MPI adapter, not needed
- **Replace** = use `nalgebra` instead

---

## Phase 1 — Core Types & Infrastructure

### Function interfaces
| File | Action | Notes |
|------|--------|-------|
| [inc/Minuit2/FCNBase.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/FCNBase.h) | Implement | → `trait FCN` |
| [inc/Minuit2/FCNGradientBase.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/FCNGradientBase.h) | Implement | → `trait FCNGradient: FCN` |
| [inc/Minuit2/GenericFunction.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/GenericFunction.h) | Skip | Rust traits replace C++ base class |
| [inc/Minuit2/FCNAdapter.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/FCNAdapter.h) | Skip | ROOT Math adapter |
| [inc/Minuit2/FCNGradAdapter.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/FCNGradAdapter.h) | Skip | ROOT Math adapter |

### Parameter management
| File | Action |
|------|--------|
| [inc/Minuit2/MinuitParameter.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MinuitParameter.h) | Implement |
| [inc/Minuit2/MnUserParameters.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnUserParameters.h) | Implement |
| [src/MnUserParameters.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnUserParameters.cxx) | Implement |
| [inc/Minuit2/MnUserParameterState.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnUserParameterState.h) | Implement |
| [src/MnUserParameterState.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnUserParameterState.cxx) | Implement |
| [inc/Minuit2/MnUserCovariance.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnUserCovariance.h) | Implement |
| [inc/Minuit2/MnUserTransformation.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnUserTransformation.h) | Implement |
| [src/MnUserTransformation.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnUserTransformation.cxx) | Implement |
| [inc/Minuit2/MnUserFcn.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnUserFcn.h) | Implement |
| [src/MnUserFcn.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnUserFcn.cxx) | Implement |

### Parameter transforms (bounded → unbounded)
| File | Action |
|------|--------|
| [inc/Minuit2/SinParameterTransformation.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/SinParameterTransformation.h) | Implement |
| [src/SinParameterTransformation.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/SinParameterTransformation.cxx) | Implement |
| [inc/Minuit2/SqrtLowParameterTransformation.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/SqrtLowParameterTransformation.h) | Implement |
| [src/SqrtLowParameterTransformation.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/SqrtLowParameterTransformation.cxx) | Implement |
| [inc/Minuit2/SqrtUpParameterTransformation.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/SqrtUpParameterTransformation.h) | Implement |
| [src/SqrtUpParameterTransformation.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/SqrtUpParameterTransformation.cxx) | Implement |
| [inc/Minuit2/MnVectorTransform.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnVectorTransform.h) | Implement |

### Result types
| File | Action |
|------|--------|
| [inc/Minuit2/FunctionMinimum.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/FunctionMinimum.h) | Implement |
| [inc/Minuit2/BasicFunctionMinimum.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/BasicFunctionMinimum.h) | Implement |
| [inc/Minuit2/MinimumState.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MinimumState.h) | Implement |
| [inc/Minuit2/BasicMinimumState.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/BasicMinimumState.h) | Implement |
| [inc/Minuit2/MinimumParameters.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MinimumParameters.h) | Implement |
| [inc/Minuit2/BasicMinimumParameters.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/BasicMinimumParameters.h) | Implement |
| [inc/Minuit2/MinimumError.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MinimumError.h) | Implement |
| [inc/Minuit2/BasicMinimumError.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/BasicMinimumError.h) | Implement |
| [src/BasicMinimumError.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/BasicMinimumError.cxx) | Implement |
| [inc/Minuit2/MinimumSeed.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MinimumSeed.h) | Implement |
| [inc/Minuit2/BasicMinimumSeed.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/BasicMinimumSeed.h) | Implement |
| [inc/Minuit2/FunctionGradient.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/FunctionGradient.h) | Implement |
| [inc/Minuit2/BasicFunctionGradient.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/BasicFunctionGradient.h) | Implement |

### Strategy & config
| File | Action |
|------|--------|
| [inc/Minuit2/MnStrategy.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnStrategy.h) | Implement |
| [src/MnStrategy.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnStrategy.cxx) | Implement |
| [inc/Minuit2/MnMachinePrecision.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnMachinePrecision.h) | Implement |
| [src/MnMachinePrecision.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnMachinePrecision.cxx) | Implement |
| [inc/Minuit2/MnTiny.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnTiny.h) | Implement |
| [src/MnTiny.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnTiny.cxx) | Implement |
| [inc/Minuit2/MnConfig.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnConfig.h) | Implement |

### Minimizer framework (abstract base)
| File | Action |
|------|--------|
| [inc/Minuit2/MnApplication.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnApplication.h) | Implement |
| [src/MnApplication.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnApplication.cxx) | Implement |
| [inc/Minuit2/FunctionMinimizer.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/FunctionMinimizer.h) | Implement |
| [inc/Minuit2/ModularFunctionMinimizer.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/ModularFunctionMinimizer.h) | Implement |
| [src/ModularFunctionMinimizer.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/ModularFunctionMinimizer.cxx) | Implement |
| [inc/Minuit2/MinimumBuilder.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MinimumBuilder.h) | Implement |
| [src/MinimumBuilder.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MinimumBuilder.cxx) | Implement |
| [inc/Minuit2/MinimumSeedGenerator.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MinimumSeedGenerator.h) | Implement |
| [inc/Minuit2/MinimumErrorUpdator.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MinimumErrorUpdator.h) | Implement |
| [inc/Minuit2/GradientCalculator.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/GradientCalculator.h) | Implement |
| [inc/Minuit2/MnFcn.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnFcn.h) | Implement |
| [src/MnFcn.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnFcn.cxx) | Implement |

### Printing/tracing
| File | Action |
|------|--------|
| [inc/Minuit2/MnPrint.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnPrint.h) | Implement |
| [src/MnPrint.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnPrint.cxx) | Implement |
| [src/MnPrintImpl.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnPrintImpl.cxx) | Implement |
| [inc/Minuit2/MnTraceObject.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnTraceObject.h) | Implement |
| [src/MnTraceObject.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnTraceObject.cxx) | Implement |

### Linear algebra (replace with `nalgebra`)
| File | Action | Notes |
|------|--------|-------|
| [inc/Minuit2/LASymMatrix.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/LASymMatrix.h) | Replace | → `nalgebra::DMatrix<f64>` |
| [inc/Minuit2/LAVector.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/LAVector.h) | Replace | → `nalgebra::DVector<f64>` |
| [inc/Minuit2/MnMatrix.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnMatrix.h) | Replace | Type aliases |
| [inc/Minuit2/ABObj.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/ABObj.h) | Replace | Expression templates → nalgebra ops |
| [inc/Minuit2/ABProd.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/ABProd.h) | Replace | |
| [inc/Minuit2/ABSum.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/ABSum.h) | Replace | |
| [inc/Minuit2/ABTypes.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/ABTypes.h) | Replace | |
| [inc/Minuit2/LaInverse.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/LaInverse.h) | Replace | |
| [inc/Minuit2/LaOuterProduct.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/LaOuterProduct.h) | Replace | |
| [inc/Minuit2/LaProd.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/LaProd.h) | Replace | |
| [inc/Minuit2/LaSum.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/LaSum.h) | Replace | |
| [inc/Minuit2/MatrixInverse.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MatrixInverse.h) | Replace | |
| [inc/Minuit2/VectorOuterProduct.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/VectorOuterProduct.h) | Replace | |
| [inc/Minuit2/StackAllocator.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/StackAllocator.h) | Skip | Custom allocator, not needed in Rust |
| [inc/Minuit2/MnRefCountedPointer.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnRefCountedPointer.h) | Skip | → `Arc`/`Rc` |
| [inc/Minuit2/MnReferenceCounter.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnReferenceCounter.h) | Skip | → `Arc`/`Rc` |

### BLAS-like routines (replace with `nalgebra`)
| File | Action |
|------|--------|
| [src/LaEigenValues.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/LaEigenValues.cxx) | Replace |
| [src/LaInnerProduct.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/LaInnerProduct.cxx) | Replace |
| [src/LaInverse.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/LaInverse.cxx) | Replace |
| [src/LaOuterProduct.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/LaOuterProduct.cxx) | Replace |
| [src/LaSumOfElements.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/LaSumOfElements.cxx) | Replace |
| [src/LaVtMVSimilarity.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/LaVtMVSimilarity.cxx) | Replace |
| [src/mndasum.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/mndasum.cxx) | Replace |
| [src/mndaxpy.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/mndaxpy.cxx) | Replace |
| [src/mnddot.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/mnddot.cxx) | Replace |
| [src/mndscal.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/mndscal.cxx) | Replace |
| [src/mndspmv.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/mndspmv.cxx) | Replace |
| [src/mndspr.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/mndspr.cxx) | Replace |
| [src/mnlsame.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/mnlsame.cxx) | Replace |
| [src/mnvert.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/mnvert.cxx) | Replace |
| [src/mnxerbla.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/mnxerbla.cxx) | Replace |

---

## Phase 2 — Simplex

| File | Action |
|------|--------|
| [inc/Minuit2/MnSimplex.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnSimplex.h) | Implement |
| [inc/Minuit2/SimplexMinimizer.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/SimplexMinimizer.h) | Implement |
| [inc/Minuit2/SimplexBuilder.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/SimplexBuilder.h) | Implement |
| [src/SimplexBuilder.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/SimplexBuilder.cxx) | Implement |
| [inc/Minuit2/SimplexParameters.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/SimplexParameters.h) | Implement |
| [src/SimplexParameters.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/SimplexParameters.cxx) | Implement |
| [inc/Minuit2/SimplexSeedGenerator.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/SimplexSeedGenerator.h) | Implement |
| [src/SimplexSeedGenerator.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/SimplexSeedGenerator.cxx) | Implement |

---

## Phase 3 — Migrad

### Core iteration
| File | Action |
|------|--------|
| [inc/Minuit2/MnMigrad.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnMigrad.h) | Implement |
| [inc/Minuit2/VariableMetricMinimizer.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/VariableMetricMinimizer.h) | Implement |
| [inc/Minuit2/VariableMetricBuilder.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/VariableMetricBuilder.h) | Implement |
| [src/VariableMetricBuilder.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/VariableMetricBuilder.cxx) | Implement |
| [inc/Minuit2/VariableMetricEDMEstimator.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/VariableMetricEDMEstimator.h) | Implement |
| [src/VariableMetricEDMEstimator.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/VariableMetricEDMEstimator.cxx) | Implement |

### Gradient calculators
| File | Action |
|------|--------|
| [inc/Minuit2/Numerical2PGradientCalculator.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/Numerical2PGradientCalculator.h) | Implement |
| [src/Numerical2PGradientCalculator.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/Numerical2PGradientCalculator.cxx) | Implement |
| [inc/Minuit2/InitialGradientCalculator.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/InitialGradientCalculator.h) | Implement |
| [src/InitialGradientCalculator.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/InitialGradientCalculator.cxx) | Implement |
| [inc/Minuit2/HessianGradientCalculator.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/HessianGradientCalculator.h) | Implement |
| [src/HessianGradientCalculator.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/HessianGradientCalculator.cxx) | Implement |
| [inc/Minuit2/AnalyticalGradientCalculator.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/AnalyticalGradientCalculator.h) | Implement |
| [src/AnalyticalGradientCalculator.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/AnalyticalGradientCalculator.cxx) | Implement |

### Line search
| File | Action |
|------|--------|
| [inc/Minuit2/MnLineSearch.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnLineSearch.h) | Implement |
| [src/MnLineSearch.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnLineSearch.cxx) | Implement |
| [inc/Minuit2/MnParabola.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnParabola.h) | Implement |
| [inc/Minuit2/MnParabolaFactory.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnParabolaFactory.h) | Implement |
| [src/MnParabolaFactory.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnParabolaFactory.cxx) | Implement |
| [inc/Minuit2/MnParabolaPoint.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnParabolaPoint.h) | Implement |

### Hessian updators
| File | Action |
|------|--------|
| [inc/Minuit2/DavidonErrorUpdator.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/DavidonErrorUpdator.h) | Implement |
| [src/DavidonErrorUpdator.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/DavidonErrorUpdator.cxx) | Implement |
| [inc/Minuit2/BFGSErrorUpdator.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/BFGSErrorUpdator.h) | Implement |
| [src/BFGSErrorUpdator.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/BFGSErrorUpdator.cxx) | Implement |

### Seed generator & helpers
| File | Action |
|------|--------|
| [inc/Minuit2/MnSeedGenerator.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnSeedGenerator.h) | Implement |
| [src/MnSeedGenerator.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnSeedGenerator.cxx) | Implement |
| [inc/Minuit2/NegativeG2LineSearch.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/NegativeG2LineSearch.h) | Implement |
| [src/NegativeG2LineSearch.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/NegativeG2LineSearch.cxx) | Implement |
| [inc/Minuit2/MnPosDef.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnPosDef.h) | Implement |
| [src/MnPosDef.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnPosDef.cxx) | Implement |
| [inc/Minuit2/MnEigen.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnEigen.h) | Implement |
| [src/MnEigen.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnEigen.cxx) | Implement |
| [src/mnteigen.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/mnteigen.cxx) | Implement |

---

## Phase 4 — Error Analysis

### Hesse
| File | Action |
|------|--------|
| [inc/Minuit2/MnHesse.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnHesse.h) | Implement |
| [src/MnHesse.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnHesse.cxx) | Implement |
| [inc/Minuit2/MnCovarianceSqueeze.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnCovarianceSqueeze.h) | Implement |
| [src/MnCovarianceSqueeze.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnCovarianceSqueeze.cxx) | Implement |
| [inc/Minuit2/MnGlobalCorrelationCoeff.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnGlobalCorrelationCoeff.h) | Implement |
| [src/MnGlobalCorrelationCoeff.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnGlobalCorrelationCoeff.cxx) | Implement |

### Minos
| File | Action |
|------|--------|
| [inc/Minuit2/MnMinos.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnMinos.h) | Implement |
| [src/MnMinos.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnMinos.cxx) | Implement |
| [inc/Minuit2/MinosError.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MinosError.h) | Implement |
| [inc/Minuit2/MnCross.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnCross.h) | Implement |
| [inc/Minuit2/MnFunctionCross.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnFunctionCross.h) | Implement |
| [src/MnFunctionCross.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnFunctionCross.cxx) | Implement |

### Scan & Contours
| File | Action |
|------|--------|
| [inc/Minuit2/MnScan.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnScan.h) | Implement |
| [src/MnScan.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnScan.cxx) | Implement |
| [inc/Minuit2/MnParameterScan.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnParameterScan.h) | Implement |
| [src/MnParameterScan.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnParameterScan.cxx) | Implement |
| [inc/Minuit2/ScanBuilder.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/ScanBuilder.h) | Implement |
| [src/ScanBuilder.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/ScanBuilder.cxx) | Implement |
| [inc/Minuit2/ScanMinimizer.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/ScanMinimizer.h) | Implement |
| [inc/Minuit2/MnContours.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnContours.h) | Implement |
| [src/MnContours.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnContours.cxx) | Implement |
| [inc/Minuit2/ContoursError.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/ContoursError.h) | Implement |

---

## Phase 5 — Extended & Stretch

### Combined minimizer (Simplex → Migrad)
| File | Action |
|------|--------|
| [inc/Minuit2/MnMinimize.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnMinimize.h) | Implement |
| [inc/Minuit2/CombinedMinimizer.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/CombinedMinimizer.h) | Implement |
| [inc/Minuit2/CombinedMinimumBuilder.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/CombinedMinimumBuilder.h) | Implement |
| [src/CombinedMinimumBuilder.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/CombinedMinimumBuilder.cxx) | Implement |

### Plotting (optional)
| File | Action |
|------|--------|
| [inc/Minuit2/MnPlot.h](https://raw.githubusercontent.com/GooFit/Minuit2/master/inc/Minuit2/MnPlot.h) | Optional |
| [src/MnPlot.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/MnPlot.cxx) | Optional |
| [src/mnbins.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/mnbins.cxx) | Optional |
| [src/mntplot.cxx](https://raw.githubusercontent.com/GooFit/Minuit2/master/src/mntplot.cxx) | Optional |

---

## Skip — ROOT/Fumili/MPI (not needed)

| File | Reason |
|------|--------|
| `inc/Math/*.h` (11 files) | ROOT Math interface |
| `inc/Fit/ParameterSettings.h` | ROOT Fit interface |
| `inc/LinkDef.h` | ROOT dictionary generator |
| `inc/TMinuit2TraceObject.h`, `src/TMinuit2TraceObject.cxx` | ROOT TObject |
| `inc/Minuit2/Minuit2Minimizer.h`, `src/Minuit2Minimizer.cxx` | ROOT Math::Minimizer adapter |
| `inc/Minuit2/Fumili*.h` (10 files), `src/Fumili*.cxx` (5 files) | Fumili method (niche, add later) |
| `inc/Minuit2/MPIProcess.h`, `src/MPIProcess.cxx` | MPI parallelism |
| `inc/Minuit2/ParametricFunction.h`, `src/ParametricFunction.cxx` | Fumili dependency |
| `src/math/*.cxx` (2 files) | ROOT Math options |
| `src/FitterUtil.h` | Internal ROOT util |

---

## Summary

| Category | Files | Action |
|----------|-------|--------|
| Core implementation scope (Phase 1-4) | ~95 | Implement |
| Replace with nalgebra | ~28 | Use crate-specific linear algebra |
| Skip (ROOT/Fumili/MPI) | ~34 | Ignore |
| Optional (Phase 5) | ~8 | Low priority |
| **Total** | **~165** | |
