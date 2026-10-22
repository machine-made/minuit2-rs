# USGS Earthquake Example

## Goal
Fit a Gutenberg-Richter relationship to a USGS earthquake catalog snapshot and export fit/residual curves.

## Data
- Local file: `examples/data/usgs/earthquakes_2025_m4p5.csv`
- Source query:
  - `https://earthquake.usgs.gov/fdsnws/event/1/query.csv?starttime=2025-01-01&endtime=2025-12-31&minmagnitude=4.5&orderby=time`
- API schema:
  - `https://earthquake.usgs.gov/fdsnws/event/1/schemas`

## Run
```bash
examples/usgs_earthquakes/run.sh
```

`run.sh` prints the total execution time and generates figures.

## Model
The fit uses:

`log10 N(M>=m) = a - b*m`

with Poisson-derived uncertainty on `log10 N`.

## Output
- `examples/usgs_earthquakes/output/usgs_gutenberg_richter_curve.csv`
  - columns: `magnitude_threshold, cumulative_count, log10_count, sigma_log10, pred_log10, residual`
- Figures:
  - `examples/usgs_earthquakes/figures/usgs_gr_fit.png`
  - `examples/usgs_earthquakes/figures/usgs_gr_counts_semilogy.png`
