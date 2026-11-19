# Matched-accuracy measurements

Run this harness from a trusted workspace with Rust, uv, and the existing ROOT reference build available.
The entry point installs its Python packages in an isolated uv environment.
It pins iminuit to 2.33.0 and leaves the repository's Python package unbuilt.

Set `compare_repo` to the absolute repository path, then run:

```sh
env -C "$compare_repo" RAYON_NUM_THREADS=2 CARGO_BUILD_JOBS=2 jobscope python3 -I "$compare_repo/verification/compare/run.py"
```

The reports are `reports/compare/matched_accuracy.csv` and `reports/compare/matched_accuracy.md`.
Neither includes timing, a host name, or absolute paths.
The elapsed time is printed to the terminal.
The default uv cache is `target/compare-scratch/uv-cache`.

Run the reference checks and verdict checks with:

```sh
env -C "$compare_repo" UV_CACHE_DIR="$compare_repo/target/compare-scratch/uv-cache" RAYON_NUM_THREADS=2 CARGO_BUILD_JOBS=2 jobscope uv run --no-project --with mpmath --with numpy --with pytest pytest "$compare_repo/verification/compare/test_cases.py" "$compare_repo/verification/compare/test_run.py"
```

`cases.py` owns the objective definitions, analytic gradients, data, limits, and steps.
It computes reference minima and Hessians with mpmath at 40 decimal digits.
The Gaussian data have seeded noise. The Poisson data contain expected bin counts without sampling noise.

`compare_runner.rs` reads one JSON object from standard input.
Required fields are `function`, `start`, and `steps`.
Optional fields are `limits`, `fixed`, `strategy`, `tolerance`, `gradient`, and `errordef`.
Dense quadratics also need `center` and `hessian`. Peak fits need `data`. Bound quadratics need `target`.
The case generator supplies every field for measurements.
The runner prints snapshots after Migrad and after Hesse, with cumulative value and gradient counters.
It also prints the difference between those counters for Hesse alone.

The Python worker counts inside the objective and gradient, using the same case data.
Migrad runs once with wrapper retries and Simplex disabled.
For strategy 3, the worker uses the bundled constructor because the public strategy setter rejects that level.
The worker checks that the constructed level is 3.
The gradient group uses tolerance 0.0001 to reach matched accuracy in both libraries.

The existing ROOT runner supplies the user-Hessian row.
That imported row retains its reported historical count and is excluded from paired verdicts.
The ROOT runner lacks the unbounded reference workload.
`root_counts.cpp` measures that workload using counters inside the function and links the existing ROOT libraries.
It does not rebuild ROOT.

The CSV stores both libraries' validity, objective, true distance, EDM, parameters, covariance, and counters.
`value_calls` and `gradient_calls` refer to the group's selected phase.
Hessian rows select the Hesse snapshot. Other groups select the Migrad snapshot.
`migrad_value_calls` always selects the earlier reading; `total_value_calls` always selects the later reading.
`hesse_calls` is the number of additional value calls. Gradient fields follow the same convention.

Failed fits retain available measurements and messages. Unsupported Rust strategy 3 rows say `not implemented`.
The report states the verdict thresholds, eligible rows, covariance convention, versions, and reference counts.
An empty matched set says `not reached` for its medians.
