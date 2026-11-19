# Changelog

## Unreleased

- Fixed Numerical gradient steps now keep ROOT's curvature floor, and Hesse gradient refinement follows ROOT's steps and stopping rule. Call counts and some covariance digits change
- Fixed Hesse now uses a supplied analytical Hessian, changing covariance values and removing extra function calls from `nfcn()`. This needs two things: the saved gradient must be analytical, and the function must report both a gradient and a Hessian. A function that reports a Hessian without a gradient still uses finite differences
- Fixed Hesse marks the diagonal covariance returned at the call limit as invalid and not accurate
- Fixed `nfcn()` after Hesse and strategy-2 Migrad now reports the running total, including seed evaluations. Analytical Migrad also keeps the seed count for other strategies. Call-limit checks include that total and stop before iteration when the seed exhausts the budget, as ROOT does. Seed generation can still exceed the requested limit, as in ROOT. Numerical Hesse now evaluates the centre once, as ROOT does
- Fixed Migrad can leave a parameter limit when the line-search parabola has no minimum. On `quadratic2`, the former upper-limit result was 3.6775; the minimum is -0.69956
- Fixed Strategy 2 now computes the full inverse Hessian before its first step, including the analytical-gradient path
- Fixed The covariance update now sums one triangle, counting each off-diagonal element once
- Fixed The line search follows ROOT Minuit2 v6-36-08 when shortening trials, opening steps, bounding parabola vertices, and halving poor trials. Zero-coordinate step bounds, adaptive point spacing, and the shared twelve-trial counter also follow ROOT
- Fixed With tolerance 0.1, the Hahn1 fit from start 2 now stops where ROOT Minuit2 stops. Relative errors against certified values are 1.347e-2 for b1 and 1.678e-2 for b4. Releases before this change ended within 1e-2 of the certified values because their line search differed from ROOT's
- Fixed Analytical seeds use the Hessian diagonal when supplied, otherwise the function's finite g2, otherwise the step-size estimate. After a negative-curvature step, a finite user g2 takes precedence over the Hessian diagonal. Limited parameters scale second derivatives by the squared transformation derivative
- Fixed Against released 0.5.5, the Rosenbrock strategy-1 run takes 140 calls, as ROOT does, where 0.5.5 took 152. The strategy-2 run takes 162 calls, as ROOT does, where 0.5.5 took 160
- Fixed The Rosenbrock function values now agree with ROOT Minuit2 v6-36-08: 3.459747538517e-05 against 3.459747567501e-05 for strategy 1, and 2.97582264946e-06 against 2.975822646186e-06 for strategy 2
- Fixed Maximum Rosenbrock covariance differences from ROOT are 1.8e-9 for strategy 1 and 1.1e-11 for strategy 2, down from 0.0981 and 0.0491 in 0.5.5
- Fixed Numerical results remain within their tolerances for all twelve reference workloads, and no workload warns about call counts. Full parameters, errors, covariance, and smaller changes on the other workloads are recorded in [reports/verification/](reports/verification/)
- Fixed The errors of parameters with limits now follow ROOT, which averages the two one-sided changes of the external value. The covariance was already equal to ROOT's. On the three limited `quadratic2` test fits, five of the six errors are smaller than before and one is larger
- Changed The minimum supported Rust version is 1.86, the oldest version that compiles the crate. Releases up to 0.5.5 did not build on 1.85 either
- Changed Machine precision now matches ROOT: `eps` is four times `f64::EPSILON`, and `eps2` is twice its square root. `MnMachinePrecision::compute_precision` returns the same precision as `new()`
- Changed Migrad traces start with one seed record at `iter: 0`. Accepted iterations start at 1 and keep increasing across continuation passes. Analytical trace call counts add subsequent evaluations to the stored seed count, so they cannot jump backward. Analytical seed calls now also appear in the reported minimum count
- Changed The comparison report in `reports/compare/` is regenerated. Its Rosenbrock covariance error with a user gradient is now 0.0399, against 0.0382 for iminuit and 0.662 before the Hesse changes

## 0.5.5

- Changed The repository, homepage and docs links now point to github.com/machine-made/minuit2-rs. The code does not change
- Changed The crate is now published to crates.io with trusted publishing instead of a stored token
- Changed The benchmark and test tooling (criterion 0.7) and the build workflows use newer versions. The library code and its dependencies do not change

## 0.5.4

- Fixed The 0.5.3 README told readers to run `pip install minuit2`, but no `minuit2` package exists on PyPI, so the command would install whatever someone else registers under that name. The README now gives only the maturin build steps. The code does not change
- Removed `.github/workflows/wheels.yml` went. minuit2 has no PyPI release, and the Python package stays a build-from-source option with maturin

## 0.5.3

- Added `README.md` gained a Disclaimer section. It states that the software comes "as is", without warranty, and that it is not validated for engineering design, certification, or safety-critical use. This release puts it on the crates.io page
- Added `.github/workflows/wheels.yml` built the `minuit2` Python package (abi3, Python 3.9 and newer) for Linux, Windows, and macOS. Its PyPI upload never succeeded, and 0.5.4 removed the workflow: minuit2 is not on PyPI
- Fixed The optional `python` feature now uses `pyo3` 0.29 and `numpy` 0.29 (was 0.28). `pyo3` 0.28 had RUSTSEC-2026-0176 (out-of-bounds read in `nth` and `nth_back` for `PyList` and `PyTuple` iterators) and RUSTSEC-2026-0177 (missing `Sync` bound on `PyCFunction::new_closure` closures). The Python tests pass on a wheel built with 0.29
- Fixed `Cargo.lock` now pins `crossbeam-epoch` 0.9.21 (was 0.9.18) for RUSTSEC-2026-0204, reached through the optional `rayon` feature, and `anyhow` 1.0.104 (was 1.0.102) for RUSTSEC-2026-0190, a development dependency only. The `Security audit` and `cargo-deny` CI jobs pass again
- Improved The README now reports measured wall-clock results against C++ Minuit2 instead of "not benchmarked": Rust is faster on USGS EQ (0.79x) and CERN dimuon (0.78x), and slower on NOAA CO2 (1.6x) and NIST StRD (2.6x), with runs pinned to one core. `examples/comparison.png` and the CERN dimuon fits were regenerated against v0.5.2
- Improved The README verification section now matches v0.5.2: the differential suite gives 12 pass, 0 warn, 0 fail. Lanczos3 and MGH09 remain the out-of-the-box failures, shared with C++ Minuit2
- Improved The README banner now loads from an absolute URL, so it also shows on crates.io
- Changed The project entered maintenance mode. The README states that the crate is feature-complete. Dependabot now opens security updates only; the `cargo-audit` and `cargo-deny` CI gates stay

## 0.5.2

- Fixed Migrad now performs ROOT's Hesse-verified convergence check (strategy ≥ 2, or strategy 1 with `dcovar > 0.05`) and continues minimization from the Hesse state with an enlarged call budget when Hesse raises EDM above tolerance, mirroring `VariableMetricBuilder.cxx`. The NIST Hahn1 dataset now reaches the certified solution from the standard starting values with plain strategy-1 Migrad — previously it returned `valid=False` (`above_max_edm`). Includes MnHesse parity fixes (cross-term `yy` caching, sag-failure failed state), a NaN/negative-EDM guard in the inner loop, and correct call-limit reporting for fits that converge inside the extended verification budget.

## 0.5.1

- Added NIST StRD certified-oracle tests: 8 datasets (Misra1a/b, Chwirut2, Rat42, Kirby2, Thurber, Gauss1, ENSO) fit to the NIST-certified parameter values
- Added Property-based metamorphic tests (`proptest`): translation, scaling, permutation, and start-point invariances over randomized inputs
- Added Expanded the `iminuit` differential drop-in harness from 20 to 27 checks
- Added README: an honest "Pure Rust vs C++ Minuit2" comparison (head-to-head function-evaluation counts from the ROOT harness) and a Testing section
- Fixed macOS CI: `cargo test --all-features` now links correctly when the pyo3 `extension-module` + `abi3` features are enabled (added `-undefined dynamic_lookup` via `.cargo/config.toml`). Dev/CI only; the published library is unaffected
- Improved Reduced allocation churn in the hot minimizer paths, with **bit-identical** numerical results (NFCN, fval, parameters, and covariance unchanged — verified by byte-identical regeneration of the ROOT differential outputs): `Numerical2PGradientCalculator` reuses two scratch vectors across the 2-point central-difference loop instead of cloning the parameter vector on every function evaluation. `mn_linesearch` reuses one scratch vector (`copy_from` + `axpy`) instead of allocating on each of its up-to-14 step evaluations. These source-level changes benefit all downstream builds. Criterion on this crate's bench suite shows ~13–29% wall-clock improvement (Rosenbrock Migrad -15%, Migrad+Hesse -16%, Minos -16%, Scan serial -29%)
- Improved Enabled `lto = true` + `codegen-units = 1` for this crate's own `release` builds (benches/examples). No effect on numerical output

## 0.5.0

- Added Python packaging: `pyproject.toml` (maturin backend, mixed layout under `python/`), `cdylib` crate-type, and an abi3 wheel. The bindings are now installable/importable (`from minuit2 import Minuit`)
- Added measured iminuit-compatible Python API subset (phase 1): method chaining (`migrad`/`hesse`/`simplex`/`minos` return `self`); `minos` populates `merrors`; new `errordef`/`strategy`/`tol` (read/write) and `fmin`/`params`/`merrors`/`nfcn`/`npar`/`nfit`/`parameters`/`accurate` properties; new `FMin`/`Param`/`MError` classes
- Added measured iminuit-compatible scan subset (phase 2): `mnprofile` (re-minimized 1D profile), `contour` (fixed-others value grid). Scan/profile/contour methods return numpy arrays
- Fixed Migrad: a free parameter started *exactly* on a limit no longer gets stuck on the bound. The bound transform's zero Jacobian at the boundary produced a zero seed gradient (EDM=0) and a premature "converged" result; the seed now performs a negative-curvature escape (`escape_negative_curvature` in `src/migrad/seed.rs`) when `g2 <= 0`, matching ROOT/iminuit. Off-boundary fits are unchanged. Regression test: `tests/limit_boundary.rs`
- Fixed all generated-CSV writers (`generate_parity_report.py`, `compare_ref_vs_rust.py`, `generate_similarity_audit.py`) to emit `LF` line terminators instead of the `csv` module's default `CRLF`, which had tripped `git diff --check`
- Changed Renamed `Minuit.contour(par_x, par_y, npoints)` -> `Minuit.mncontour(x, y, *, cl=None, size=100)` (Minos contour boundary points, returned as a closed `(size+1, 2)` numpy array). `cl` scaling is not yet implemented
- Changed Renamed `Minuit.scan(param, nsteps, low, high)` -> `Minuit.profile(vname, *, size=100, bound=2.0, subtract_min=False)` (1D profile with other params held fixed, returns `(x, fvals)` numpy arrays)
- Changed `Minuit.scan(...)` now refers to iminuit's brute-force global minimizer and currently raises `NotImplementedError` (deferred), rather than silently doing a 1D profile under the old name
- Changed Slimmed the published crate from 8.4 MiB (237 files) to ~97 KiB by switching `Cargo.toml` from an `exclude` blocklist to an anchored `include` allowlist. The crate now ships only the library/example/bench source, `README`, `CHANGELOG`, and the two license files; example datasets, verification reports, coverage dumps, figures, and dev tooling are no longer published
- Changed Licensed the crate as `LGPL-2.1-or-later`, matching the upstream ROOT Minuit2 that it reimplements. (Earlier 0.5.0 development drafts had used `MIT OR Apache-2.0`; since this is a reimplementation of LGPL code rather than a clean-room work, the honest and safe license is the upstream one. No permissively-licensed version was ever published.) Added the `LGPL-2.1` text as `LICENSE` and a `NOTICE` crediting the original CERN/PH-SFT authors
- Changed Clarified public provenance wording: this crate is an independent Rust reimplementation validated against ROOT Minuit2 as a numerical reference; ROOT Minuit2 is not a source dependency. Removed loose "clean-room" phrasing
- Changed Updated the release workflow gate to require release tags to be reachable from `origin/main` and to run `cargo package --locked` before `cargo publish --locked`
- Changed Retargeted the scheduled scientific-demo artifact commit gate from `master` to `main`

## 0.4.2

- Added a `cargo-deny` CI gate alongside the existing `cargo-audit` check. `deny.toml` enforces a permissive license allowlist, denies unknown registries and git sources, and gates the advisory database against the entire dependency graph
- Added `examples/data/SHA256SUMS` pinning the SHA-256 of every scientific reference data file to the blob tracked in git
- Fixed Upgraded `pyo3` 0.28.0 → 0.28.3. The 0.28.0 lockfile entry was yanked on crates.io and carried RUSTSEC-2026-0013 (type confusion in `#[pyclass(extends=<native>)]` with the `abi3` feature on Python 3.12+). The Python bindings in this crate do not use `extends=`, so the library was not materially exposed, but the yanked version was replaced anyway
- Fixed Hardened all four GitHub Actions workflows against shell injection: every `${{ inputs.* }}` and `${{ github.event.inputs.* }}` used in a `run:` block is now routed through a step-level `env:` map so user input lands in bash as data, never as program source. Numeric inputs are validated at the top of the step
- Fixed Hardened the release workflow against awk code injection through a maliciously-named tag. `VERSION` is now passed to awk via `-v` (data binding) and matched with `index()` (literal prefix), not spliced into the awk program
- Fixed Pinned every third-party GitHub Action by commit SHA across `ci.yml`, `release.yml`, `scientific-demo.yml`, and `scientific-demo-scheduled.yml` (8 distinct actions). Floating refs like `@v4` and `@stable` can be force-moved upstream; SHA pins close that supply-chain vector
- Fixed Split scientific-demo workflows so the benchmark job runs with `permissions: contents: read`. Only the follow-up `commit-refresh` / `commit-artifacts` job carries `contents: write`, and only when gated on `master` (scheduled) or the `auto_commit` input (manual)
- Fixed Release job now rejects any tag whose commit is not an ancestor of `origin/master`: a release cannot ship from an unreviewed branch
- Changed Security hardening and supply-chain release. No runtime API or behavior changes; the Rust library surface is identical to 0.4.1
- Changed Rewrote `scripts/fetch_scientific_demo_data.sh` to verify content against the manifest, stage downloads to a sibling tempfile, and only move into place after a successful checksum match. Restored TLS verification on every download (no more `curl -k`); `--proto '=https'` refuses any non-HTTPS scheme including on redirects. Any mismatch or missing manifest entry aborts with a clear diagnostic
- Changed Refreshed the NOAA CO₂ reference snapshot to the current upstream (additive February 2026 row and small retroactive uncertainty revisions, consistent with NOAA post-analysis practice)
- Removed Miri CI job (nightly toolchain incompatibility with `quote` crate build script)

## 0.4.1

- Added ROOT v6-36-08 verification harness with parity baseline for cross-validation
- Added Executed-surface gate with expanded verification workloads
- Added Scientific demonstration examples (Rosenbrock, Himmelbayas, Rastrigin, sphere functions)
- Added Fair solver benchmark workflow with automated artifact refresh
- Added scheduled scientific benchmark workflow with auto-commit to master
- Added verification snapshot and developer guide sections
- Fixed Security: upgraded PyO3 0.23 → 0.28 to resolve RUSTSEC-2025-0020
- Fixed Resolved all clippy warnings for clean CI lint gate
- Fixed Corrected two API errors in README examples
- Fixed Applied cargo fmt across entire codebase
- Improved `MnScan`, `MnContours`, and `MnMinos` internal implementations
- Improved Expanded README with comprehensive user guide integrating DOC.md content
- Changed Internal refactoring: tightened executed-surface mapping and reduced false positives
- Changed Enhanced `MnMinimize` internal logic and state handling
- Changed Expanded test coverage (93 tests, up from 92)
- Changed Finalized examples cleanup for cleaner CI workflow
- Removed main branch alias, targeting master only

## 0.4.0

- Added Analytical gradient support via `AnalyticalGradientCalculator` trait — users can provide exact gradients for faster convergence
- Added Parallel parameter scans with rayon (opt-in via `parallel` feature) — simultaneous evaluation across scan ranges
- Added Python bindings via PyO3 (opt-in via `python` feature) — call minuit2 minimizers from Python
- Added `MnMinimize` combined minimizer — automatic Simplex → Migrad fallback strategy
- Added Criterion benchmarks for performance regression tracking
- Added Robustness test suite — NaN/Inf resilience, boundary edge cases, high-dimensional stress tests
- Added error analysis examples
- Improved Streamlined README and DOC.md for v0.4.0 features
- Improved CHANGELOG.md now used for GitHub release notes
- Changed Enhanced panic safety in core minimizers with comprehensive edge-case handling
- Changed CI/CD improvements: security audit, code coverage, benchmark jobs

## 0.3.0

- Added `MnHesse` — Full Hessian computation for accurate parameter errors and covariance
- Added `MnMinos` — Profile-likelihood asymmetric confidence intervals
- Added `MnScan` — 1D parameter scans with minimum tracking
- Added `MnContours` — 2D confidence contour computation
- Added Global correlation coefficients and covariance squeeze utilities
- Added 76 tests (48 unit + 8 Migrad + 6 Simplex + 4 Hesse + 3 Minos + 3 Scan + 2 Contours + 2 doctests)

## 0.2.0

- Added `MnMigrad` — quasi-Newton minimizer with DFP/BFGS inverse Hessian update
- Added Quadratic interpolation (`MnParabola`, `MnParabolaPoint`)
- Added Parabolic line search (`mn_linesearch`)
- Added Positive-definite matrix correction (`make_pos_def`)
- Added Two-point central difference gradient (`Numerical2PGradientCalculator`)
- Added Full Migrad implementation: seed, builder, minimizer, public API
- Added 58 tests (42 unit + 8 Migrad + 6 Simplex + 2 doctests)

## 0.1.0

- Added `FCN` trait with blanket impl for closures
- Added Parameter management: free, bounded, fixed, constant
- Added Sin/sqrt parameter transforms for bounded optimization
- Added `MnSimplex` minimizer (Minuit variant of Nelder-Mead)
- Added `FunctionMinimum` result type with Display impl
- Added `MnStrategy` (low/medium/high optimization presets)
- Added `nalgebra` for all linear algebra
