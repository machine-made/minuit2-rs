# Scientific Demo Examples

This repository now includes four real-data scientific demos:

1. `examples/noaa_co2/run.sh`
2. `examples/nist_strd/run.sh`
3. `examples/usgs_earthquakes/run.sh`
4. `examples/cern_dimuon/run.sh`

Before running examples, download the reference data:

```bash
scripts/fetch_scientific_demo_data.sh
```

Each demo writes curve/residual CSV output into its own folder:
- `examples/noaa_co2/output/`
- `examples/nist_strd/output/`
- `examples/usgs_earthquakes/output/`
- `examples/cern_dimuon/output/`
Each `run.sh` prints total execution time and generates at least one figure.

Per-example documentation:
- `examples/noaa_co2/README.md`
- `examples/nist_strd/README.md`
- `examples/usgs_earthquakes/README.md`
- `examples/cern_dimuon/README.md`

Aggregate C++ vs Rust timing comparison:

```bash
examples/run_comparison.sh 41 9 5000 10000 9
# optional:
# examples/run_comparison.sh 41 9 5000 10000 9 2 0.10 1
#                                      batches core trim strict_env
```

Artifacts:
- `examples/comparison.png`
- `examples/output/comparison_timings.csv`
- `examples/output/comparison_samples.csv`
- `examples/output/comparison_environment.json`

The comparison now uses in-process solver timing:
- data are loaded/prepared once per benchmark job
- only solver computation is timed (warmups excluded)
- jobs are interleaved across C++/Rust and all cases
- robust stats are included (`trimmed median` + `MAD`)
- Linux CPU pinning is supported via `taskset` when a core is provided
- environment checks (governor/power mode) are logged

GitHub Actions workflow:
- `.github/workflows/scientific-demo.yml` runs tests, regenerates benchmark figures/artifacts, uploads them, and can optionally auto-commit updated benchmark outputs.
- `.github/workflows/scientific-demo-scheduled.yml` can be run manually to refresh and auto-commit benchmark artifacts on `main`.

![C++ vs Rust timing comparison](comparison.png)
