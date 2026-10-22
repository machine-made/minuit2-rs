# CERN Dimuon Example

## Goal
Fit the Z-peak region of dimuon mass spectra using two CERN Open Data CSV datasets.

## Data
- `examples/data/cern/MuRun2010B_0.csv`
- `examples/data/cern/Zmumu.csv`

Sources:
- `https://opendata.cern.ch/api/records/700`
- `https://opendata.cern.ch/api/records/545`

Downloaded from CERN EOS public endpoint by:
- `scripts/fetch_scientific_demo_data.sh`

## Run
```bash
examples/cern_dimuon/run.sh
```

`run.sh` prints the total execution time and generates figures.

## Method
- Dataset 1 (`MuRun2010B_0.csv`):
  - uses direct invariant mass column `M`
  - fits the J/psi region (2-5 GeV) with Gaussian + linear background
- Dataset 2 (`Zmumu.csv`):
  - reconstructs mass from `(pt, eta, phi)` pairs with the massless approximation
  - fits the Z region (60-120 GeV) with Gaussian + linear background

## Output
- `examples/cern_dimuon/output/cern_murun2010b0_jpsi_curve.csv`
- `examples/cern_dimuon/output/cern_zmumu_zpeak_curve.csv`
- Figures:
  - `examples/cern_dimuon/figures/cern_murun_jpsi_fit.png`
  - `examples/cern_dimuon/figures/cern_zmumu_z_fit.png`
