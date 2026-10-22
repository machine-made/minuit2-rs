# NIST StRD Example

## Goal
Validate nonlinear regression fits against NIST Statistical Reference Datasets (StRD) certified values.

## Data
- `examples/data/nist/Misra1a.dat`
- `examples/data/nist/Hahn1.dat`
- `examples/data/nist/Rat43.dat`

Primary references:
- `https://www.nist.gov/itl/sed/products-services/statistical-reference-data-sets-strd`
- `https://www.itl.nist.gov/div898/strd/nls/data/misra1a.shtml`
- `https://www.itl.nist.gov/div898/strd/nls/data/ratkowsky3.shtml`

## Run
```bash
examples/nist_strd/run.sh
```

`run.sh` prints the total execution time and generates figures.

## Fits
- Misra1a: exponential model (2 parameters)
- Hahn1: rational cubic/cubic model (7 parameters)
- Rat43: sigmoidal growth model (4 parameters)

The example prints fitted vs certified parameter values and RSS differences.
`Misra1a` and `Rat43` match certified values closely with default settings; `Hahn1`
is stiffer and can require tighter tuning/restarts for certified-level agreement.

## Output
- `examples/nist_strd/output/nist_summary.csv`
- `examples/nist_strd/output/nist_misra1a_curve.csv`
- `examples/nist_strd/output/nist_hahn1_curve.csv`
- `examples/nist_strd/output/nist_rat43_curve.csv`
- Figures:
  - `examples/nist_strd/figures/nist_strd_fits.png`
  - `examples/nist_strd/figures/nist_strd_rel_error.png`
