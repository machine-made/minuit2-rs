# NOAA CO2 Example

## Goal
Fit a harmonic trend model to Mauna Loa monthly atmospheric CO2 data and export curve/residual data.

## Data
- Local file: `examples/data/noaa/co2_mm_mlo.csv`
- Source: NOAA GML
  - `https://gml.noaa.gov/ccgg/trends/`
  - `https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_mlo.csv`

## Run
```bash
examples/noaa_co2/run.sh
```

`run.sh` prints the total execution time and generates figures.

## Model
The fitted model is:

`C(t)=a0+a1*t+a2*t^2+b1*sin(2*pi*t)+c1*cos(2*pi*t)+b2*sin(4*pi*t)+c2*cos(4*pi*t)+d1*t*sin(2*pi*t)`

where `t` is years from the first observation.

## Output
- `examples/noaa_co2/output/noaa_co2_curve.csv`
  - columns: `decimal_date, observed, uncertainty, fitted, residual`
- Figures:
  - `examples/noaa_co2/figures/noaa_co2_fit.png`
  - `examples/noaa_co2/figures/noaa_co2_residual_hist.png`
