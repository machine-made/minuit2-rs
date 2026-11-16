# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Fixed

- Migrad no longer ends on a limit with `valid=True` in a fit of two or more
  parameters. This happened when a doubly limited parameter started on its
  upper limit and the line search met a parabola without a minimum. On
  `quadratic2` the fit gave fval 3.6775 with `y` stuck on its limit. The minimum
  is -0.69956. The lower-limit start followed in the precision change below.
- The line search changes in three places. A negative parabola vertex now goes
  to the step limit (`MnLineSearch.cxx:137` in ROOT Minuit2 v6.36.08). A
  parabola without a minimum now opens the step from the best point
  (`MnLineSearch.cxx:224-229`). After an open step, the window for the next
  vertex follows the best point (`MnLineSearch.cxx:219` and `:234-245`).
- The negative-vertex change also applies to fits that do not start on a limit.
- The twelve reference workloads stayed inside the stated limits. Changed
  lists the values that moved.

### Changed

- Strategy 2 replaced the diagonal seed with the inverse Hessian before the first step.
- The covariance update summed one triangle. Each off-diagonal was counted once.
- The line search kept a finite vertex inside the window around the best step.
- A trial worse than the three stored points was halved toward the best step.
- The line search shortened the next trial while the first steps did not lower the function.
- `rosenbrock2_migrad` and `rosenbrock2_minimize` fval moved from 3.294567072291e-05 to 3.459740876221e-05. ROOT fval was 3.459747567500577e-05.
- edm moved from 3.423047369382e-05 to 3.451709943513e-05. ROOT edm was 3.451850739606217e-05.
- nfcn moved from 142 to 143. ROOT nfcn was 140.
- Parameters moved from 0.994468418655 and 0.989120644009 to 0.994149191534 and 0.988272162645. ROOT parameters were 0.994149185750 and 0.988272151209.
- `rosenbrock2_migrad_strategy2` fval moved from 3.294658278025e-05 to 2.975823963030e-06. ROOT fval was 2.975822646186e-06.
- edm moved from 3.423143307160e-05 to 2.980256746000e-06. ROOT edm was 2.980639536818e-06.
- nfcn moved from 150 to 168. ROOT nfcn was 162.
- Parameters moved from 0.994468340792 and 0.989120490797 to 0.998323497264 and 0.996609164995. ROOT parameters were 0.998323497642 and 0.996609165756.
- The analytical seed stored dxdi squared times the Hessian diagonal when the function had a Hessian.
- The seed stored dxdi squared times the g2 when the function had a g2 and no Hessian.
- After a step, the escape stored the g2 of the function when that g2 was finite.
- After a step, the escape stored the Hessian diagonal when the function had a Hessian and no finite g2.
- A function with no g2 and no Hessian kept the step-size estimate.
- The seed ignored a supplied g2 that held a NaN or an infinity, or that was shorter than the list of external parameters.
- The seed then used the Hessian diagonal or the step-size estimate.
- The twelve reference workloads did not supply a g2. Their values did not move.
- The machine precision followed ROOT Minuit2. `eps` became four times
  `f64::EPSILON`. `eps2` became two times the square root of `eps`.
- `MnMachinePrecision::compute_precision` set `f64::EPSILON`. It became the
  value of `new()`, four times `f64::EPSILON`, as in ROOT's
  `ComputePrecision`.
- The forward cap of a line-search vertex grew with the best step on every
  pass, up to the overall step limit of 1000.
- `rosenbrock2_migrad` and `rosenbrock2_minimize` fval moved from
  1.500342419704e-05 to 3.294567072291e-05. ROOT fval was 3.459747567500577e-05.
- edm moved from 1.527164793786e-05 to 3.423047369382e-05. ROOT edm was
  3.451850739606217e-05.
- nfcn moved from 152 to 142. ROOT nfcn was 140.
- Parameters moved from 0.996166325646 and 0.992402698612 to 0.994468418655
  and 0.989120644009. ROOT parameters were 0.994149185750 and 0.988272151209.
- The covariance of `rosenbrock2_migrad` and `rosenbrock2_minimize` moved away
  from ROOT. Elements [0][0], [0][1] and [1][1] moved from 1.018333, 2.028876
  and 4.052231 to 1.042177, 2.072850 and 4.132819. ROOT elements were 0.997644,
  1.983639 and 3.954115.
- The largest covariance difference to ROOT grew from 0.0981 to 0.1787. The
  limit is 0.5.
- The errors moved away from ROOT too. They moved from 1.009125 and 2.013015
  to 1.020871 and 2.032933. ROOT errors were 0.998821 and 1.988496.
- `rosenbrock2_migrad_strategy2` fval moved from 1.499992182613e-05 to
  3.294658278025e-05. ROOT fval was 2.975822646186212e-06.
- edm moved from 1.526804509045e-05 to 3.423143307160e-05. ROOT edm was
  2.980639536818047e-06.
- nfcn moved from 160 to 150. ROOT nfcn was 162.
- Parameters moved from 0.996166763841 and 0.992403558744 to 0.994468340792
  and 0.989120490797. ROOT parameters were 0.998323497642 and 0.996609165756.
- `quadratic2_lower_limited_migrad` fval moved from -0.6995599061026854 to
  -0.6995599061088404. ROOT fval was -0.6995599061088406.
- nfcn moved from 46 to 48. ROOT nfcn was 45.
- Parameters moved from 1.307091608844 and -2.048903410818 to 1.307091639211
  and -2.048903401202. ROOT parameters were 1.307091639211 and -2.048903401202.
- `quadratic2_no_g2_migrad` nfcn moved from 9 to 10. ROOT nfcn was 2.
- Parameters moved from 0.9999999999999992 and -2.0000000000000004 to 1 and
  -2. ROOT parameters were 1 and -2.
- `quadratic3_fixx_migrad` edm moved from 1.492603334512e-05 to
  1.492603334513e-05. fval and the parameters did not move.
- `quadratic3_fixx_hesse` edm moved from 1.49260333418e-05 to
  1.492603334305e-05. fval and the parameters did not move.
- `quadratic3_fixx_hesse` covariance elements [0][0], [0][1] and [1][1] moved
  from 0.10162601624357, -0.00406504064900 and 0.01016260162453 to
  0.10162601624459, -0.00406504065053 and 0.01016260162465.
- The errors of the two free parameters moved from 0.31878835650565 and
  0.10080972981083 to 0.31878835650724 and 0.10080972981143. The error of the
  fixed parameter did not move.
- ROOT covariance elements were 0.10162601725420, -0.00406504077661 and
  0.01016260163979. ROOT errors were 0.31878835809076 and 0.10080972988648.
  The moves are 1.2e-13 to 1.6e-12.
- Minos on parameter 0 moved its lower end from -1.0082332377850711 to
  -1.0082332377915708 and its upper end from 0.9974487592881941 to
  0.997448759281765.
- Minos on parameter 1 moved its lower end from -0.5022370303122039 to
  -0.5022370303130429 and its upper end from 0.5006039697188206 to
  0.5006039697180208.
- `quadratic2_minos_p0`, `quadratic2_minos_p1`, `quadratic2_scan_p0` and
  `quadratic2_contours_01` shared one move of fval and of the parameters.
- fval moved from -0.6995269602635593 to -0.6995269602634814. ROOT fval was
  -0.6995269602629307.
- Parameters moved from 1.312746182715 and -2.048209238576 to 1.312746182721
  and -2.048209238575. ROOT parameters were 1.312746182767 and -2.048209238569.
- edm moved from 3.306487787324e-05 to 3.306487793495e-05 on Minos and the
  contours. The scan edm moved from 3.306487786942e-05 to 3.306487794182e-05.
- `quadratic2_simplex` did not move.
- Errors and covariance moved on every workload except `quadratic2_simplex`.
- Existing fits moved within the stated limits.
- `rosenbrock2_migrad_strategy2` moved away from ROOT. The largest parameter
  difference grew from 4.2e-3 to 7.49e-3. The limit is 8e-3.
- The covariance difference to ROOT of `rosenbrock2_migrad_strategy2` grew from
  0.0491 to 0.1297. The limit is 0.5. The fval difference grew from 1.2e-5 to
  3.0e-5.
- `quadratic2_no_g2_migrad` has a relative nfcn difference to ROOT of 0.800,
  at its warn limit of 0.8.

## [0.5.5] - 2026-10-31

### Changed

- The repository, homepage and docs links now point to
  github.com/machine-made/minuit2-rs. The code does not change.
- The crate is now published to crates.io with trusted publishing instead of a
  stored token.
- The benchmark and test tooling (criterion 0.7) and the build workflows use
  newer versions. The library code and its dependencies do not change.

## [0.5.4] - 2026-10-10

### Removed

- `.github/workflows/wheels.yml` went. minuit2 has no PyPI release, and the
  Python package stays a build-from-source option with maturin.

### Fixed

- The 0.5.3 README told readers to run `pip install minuit2`, but no `minuit2`
  package exists on PyPI, so the command would install whatever someone else
  registers under that name. The README now gives only the maturin build
  steps. The code does not change.

## [0.5.3] - 2026-10-10

### Added

- `README.md` gained a Disclaimer section. It states that the software comes
  "as is", without warranty, and that it is not validated for engineering
  design, certification, or safety-critical use. This release puts it on the
  crates.io page.
- `.github/workflows/wheels.yml` built the `minuit2` Python package (abi3,
  Python 3.9 and newer) for Linux, Windows, and macOS. Its PyPI upload never
  succeeded, and 0.5.4 removed the workflow: minuit2 is not on PyPI.

### Changed

- The project entered maintenance mode. The README states that the crate is
  feature-complete. Dependabot now opens security updates only; the
  `cargo-audit` and `cargo-deny` CI gates stay.

### Security

- The optional `python` feature now uses `pyo3` 0.29 and `numpy` 0.29 (was
  0.28). `pyo3` 0.28 had RUSTSEC-2026-0176 (out-of-bounds read in `nth` and
  `nth_back` for `PyList` and `PyTuple` iterators) and RUSTSEC-2026-0177
  (missing `Sync` bound on `PyCFunction::new_closure` closures). The Python
  tests pass on a wheel built with 0.29.
- `Cargo.lock` now pins `crossbeam-epoch` 0.9.21 (was 0.9.18) for
  RUSTSEC-2026-0204, reached through the optional `rayon` feature, and `anyhow`
  1.0.104 (was 1.0.102) for RUSTSEC-2026-0190, a development dependency only.
  The `Security audit` and `cargo-deny` CI jobs pass again.

### Docs

- The README now reports measured wall-clock results against C++ Minuit2
  instead of "not benchmarked": Rust is faster on USGS EQ (0.79x) and CERN
  dimuon (0.78x), and slower on NOAA CO2 (1.6x) and NIST StRD (2.6x), with
  runs pinned to one core. `examples/comparison.png` and the CERN dimuon fits
  were regenerated against v0.5.2.
- The README verification section now matches v0.5.2: the differential suite
  gives 12 pass, 0 warn, 0 fail. Lanczos3 and MGH09 remain the out-of-the-box
  failures, shared with C++ Minuit2.
- The README banner now loads from an absolute URL, so it also shows on
  crates.io.

## [0.5.2] - 2026-06-11

### Fixed

- Migrad now performs ROOT's Hesse-verified convergence check (strategy ≥ 2, or
  strategy 1 with `dcovar > 0.05`) and continues minimization from the Hesse
  state with an enlarged call budget when Hesse raises EDM above tolerance,
  mirroring `VariableMetricBuilder.cxx`. The NIST Hahn1 dataset now reaches the
  certified solution from the standard starting values with plain strategy-1
  Migrad — previously it returned `valid=False` (`above_max_edm`). Includes
  MnHesse parity fixes (cross-term `yy` caching, sag-failure failed state), a
  NaN/negative-EDM guard in the inner loop, and correct call-limit reporting
  for fits that converge inside the extended verification budget.
  (minuit2-rs-j2j)

## [0.5.1] - 2026-06-02

### Performance

- Reduced allocation churn in the hot minimizer paths, with **bit-identical**
  numerical results (NFCN, fval, parameters, and covariance unchanged — verified
  by byte-identical regeneration of the ROOT differential outputs):
  - `Numerical2PGradientCalculator` reuses two scratch vectors across the
    2-point central-difference loop instead of cloning the parameter vector on
    every function evaluation.
  - `mn_linesearch` reuses one scratch vector (`copy_from` + `axpy`) instead of
    allocating on each of its up-to-14 step evaluations.

  These source-level changes benefit all downstream builds. Criterion on this
  crate's bench suite shows ~13–29% wall-clock improvement (Rosenbrock Migrad
  -15%, Migrad+Hesse -16%, Minos -16%, Scan serial -29%).
- Enabled `lto = true` + `codegen-units = 1` for this crate's own `release`
  builds (benches/examples). No effect on numerical output.

### Added

- NIST StRD certified-oracle tests: 8 datasets (Misra1a/b, Chwirut2, Rat42,
  Kirby2, Thurber, Gauss1, ENSO) fit to the NIST-certified parameter values.
- Property-based metamorphic tests (`proptest`): translation, scaling,
  permutation, and start-point invariances over randomized inputs.
- Expanded the `iminuit` differential drop-in harness from 20 to 27 checks.
- README: an honest "Pure Rust vs C++ Minuit2" comparison (head-to-head
  function-evaluation counts from the ROOT harness) and a Testing section.

### Fixed

- macOS CI: `cargo test --all-features` now links correctly when the pyo3
  `extension-module` + `abi3` features are enabled (added `-undefined
  dynamic_lookup` via `.cargo/config.toml`). Dev/CI only; the published library
  is unaffected.

## [0.5.0] - 2026-06-01

### Added

- Python packaging: `pyproject.toml` (maturin backend, mixed layout under
  `python/`), `cdylib` crate-type, and an abi3 wheel. The bindings are now
  installable/importable (`from minuit2 import Minuit`).
- measured iminuit-compatible Python API subset (phase 1): method chaining
  (`migrad`/`hesse`/`simplex`/`minos` return `self`); `minos` populates
  `merrors`; new `errordef`/`strategy`/`tol` (read/write) and
  `fmin`/`params`/`merrors`/`nfcn`/`npar`/`nfit`/`parameters`/`accurate`
  properties; new `FMin`/`Param`/`MError` classes.
- measured iminuit-compatible scan subset (phase 2): `mnprofile` (re-minimized 1D profile),
  `contour` (fixed-others value grid). Scan/profile/contour methods return
  numpy arrays.

### Fixed

- Migrad: a free parameter started *exactly* on a limit no longer gets stuck on
  the bound. The bound transform's zero Jacobian at the boundary produced a zero
  seed gradient (EDM=0) and a premature "converged" result; the seed now performs
  a negative-curvature escape (`escape_negative_curvature` in `src/migrad/seed.rs`)
  when `g2 <= 0`, matching ROOT/iminuit. Off-boundary fits are unchanged.
  Regression test: `tests/limit_boundary.rs`.

### Changed (breaking, Python API)

- Renamed `Minuit.contour(par_x, par_y, npoints)` -> `Minuit.mncontour(x, y, *,
  cl=None, size=100)` (Minos contour boundary points, returned as a closed
  `(size+1, 2)` numpy array). `cl` scaling is not yet implemented.
- Renamed `Minuit.scan(param, nsteps, low, high)` -> `Minuit.profile(vname, *,
  size=100, bound=2.0, subtract_min=False)` (1D profile with other params held
  fixed, returns `(x, fvals)` numpy arrays).
- `Minuit.scan(...)` now refers to iminuit's brute-force global minimizer and
  currently raises `NotImplementedError` (deferred), rather than silently doing
  a 1D profile under the old name.

### Changed

- Slimmed the published crate from 8.4 MiB (237 files) to ~97 KiB by switching
  `Cargo.toml` from an `exclude` blocklist to an anchored `include` allowlist.
  The crate now ships only the library/example/bench source, `README`,
  `CHANGELOG`, and the two license files; example datasets, verification
  reports, coverage dumps, figures, and dev tooling are no longer published.
- Fixed all generated-CSV writers (`generate_parity_report.py`,
  `compare_ref_vs_rust.py`, `generate_similarity_audit.py`) to emit `LF` line
  terminators instead of the `csv` module's default `CRLF`, which had tripped
  `git diff --check`.
- Licensed the crate as `LGPL-2.1-or-later`, matching the upstream ROOT Minuit2
  that it reimplements. (Earlier 0.5.0 development drafts had used
  `MIT OR Apache-2.0`; since this is a reimplementation of LGPL code rather than
  a clean-room work, the honest and safe license is the upstream one. No
  permissively-licensed version was ever published.) Added the `LGPL-2.1` text
  as `LICENSE` and a `NOTICE` crediting the original CERN/PH-SFT authors.
- Clarified public provenance wording: this crate is an independent Rust
  reimplementation validated against ROOT Minuit2 as a numerical reference; ROOT
  Minuit2 is not a source dependency. Removed loose "clean-room" phrasing.
- Updated the release workflow gate to require release tags to be reachable
  from `origin/main` and to run `cargo package --locked` before
  `cargo publish --locked`.
- Retargeted the scheduled scientific-demo artifact commit gate from
  `master` to `main`.

## [0.4.2] - 2026-04-21

Security hardening and supply-chain release. No runtime API or behavior
changes; the Rust library surface is identical to 0.4.1.

### Security

- Upgraded `pyo3` 0.28.0 → 0.28.3. The 0.28.0 lockfile entry was yanked on
  crates.io and carried RUSTSEC-2026-0013 (type confusion in
  `#[pyclass(extends=<native>)]` with the `abi3` feature on Python 3.12+).
  The Python bindings in this crate do not use `extends=`, so the library
  was not materially exposed, but the yanked version was replaced anyway.
- Hardened all four GitHub Actions workflows against shell injection:
  every `${{ inputs.* }}` and `${{ github.event.inputs.* }}` used in a
  `run:` block is now routed through a step-level `env:` map so user
  input lands in bash as data, never as program source. Numeric inputs
  are validated at the top of the step.
- Hardened the release workflow against awk code injection through a
  maliciously-named tag. `VERSION` is now passed to awk via `-v` (data
  binding) and matched with `index()` (literal prefix), not spliced into
  the awk program.
- Pinned every third-party GitHub Action by commit SHA across `ci.yml`,
  `release.yml`, `scientific-demo.yml`, and `scientific-demo-scheduled.yml`
  (8 distinct actions). Floating refs like `@v4` and `@stable` can be
  force-moved upstream; SHA pins close that supply-chain vector.
- Added a `cargo-deny` CI gate alongside the existing `cargo-audit`
  check. `deny.toml` enforces a permissive license allowlist, denies
  unknown registries and git sources, and gates the advisory database
  against the entire dependency graph.
- Split scientific-demo workflows so the benchmark job runs with
  `permissions: contents: read`. Only the follow-up `commit-refresh` /
  `commit-artifacts` job carries `contents: write`, and only when gated
  on `master` (scheduled) or the `auto_commit` input (manual).
- Release job now rejects any tag whose commit is not an ancestor of
  `origin/master`: a release cannot ship from an unreviewed branch.

### Data

- Added `examples/data/SHA256SUMS` pinning the SHA-256 of every scientific
  reference data file to the blob tracked in git.
- Rewrote `scripts/fetch_scientific_demo_data.sh` to verify content
  against the manifest, stage downloads to a sibling tempfile, and only
  move into place after a successful checksum match. Restored TLS
  verification on every download (no more `curl -k`); `--proto '=https'`
  refuses any non-HTTPS scheme including on redirects. Any mismatch or
  missing manifest entry aborts with a clear diagnostic.
- Refreshed the NOAA CO₂ reference snapshot to the current upstream
  (additive February 2026 row and small retroactive uncertainty
  revisions, consistent with NOAA post-analysis practice).

### CI/CD

- Removed Miri CI job (nightly toolchain incompatibility with `quote`
  crate build script).

## [0.4.1] - 2026-02-12

### Added

- ROOT v6-36-08 verification harness with parity baseline for cross-validation
- Executed-surface gate with expanded verification workloads
- Scientific demonstration examples (Rosenbrock, Himmelbayas, Rastrigin, sphere functions)
- Fair solver benchmark workflow with automated artifact refresh

### Changed

- Internal refactoring: tightened executed-surface mapping and reduced false positives
- Enhanced `MnMinimize` internal logic and state handling
- Improved `MnScan`, `MnContours`, and `MnMinos` internal implementations
- Expanded test coverage (93 tests, up from 92)

### Fixed

- Security: upgraded PyO3 0.23 → 0.28 to resolve RUSTSEC-2025-0020
- Resolved all clippy warnings for clean CI lint gate
- Corrected two API errors in README examples
- Applied cargo fmt across entire codebase

### CI/CD

- Added scheduled scientific benchmark workflow with auto-commit to master
- Removed main branch alias, targeting master only
- Finalized examples cleanup for cleaner CI workflow

### Docs

- Expanded README with comprehensive user guide integrating DOC.md content
- Added verification snapshot and developer guide sections

## [0.4.0] - 2026-02-09

### Added

- Analytical gradient support via `AnalyticalGradientCalculator` trait — users can provide exact gradients for faster convergence
- Parallel parameter scans with rayon (opt-in via `parallel` feature) — simultaneous evaluation across scan ranges
- Python bindings via PyO3 (opt-in via `python` feature) — call minuit2 minimizers from Python
- `MnMinimize` combined minimizer — automatic Simplex → Migrad fallback strategy
- Criterion benchmarks for performance regression tracking
- Robustness test suite — NaN/Inf resilience, boundary edge cases, high-dimensional stress tests

### Changed

- Enhanced panic safety in core minimizers with comprehensive edge-case handling
- CI/CD improvements: security audit, code coverage, benchmark jobs

### Docs

- Streamlined README and DOC.md for v0.4.0 features
- Added error analysis examples
- CHANGELOG.md now used for GitHub release notes

## [0.3.0] - 2026-02-09

### Added

- `MnHesse` — Full Hessian computation for accurate parameter errors and covariance
- `MnMinos` — Profile-likelihood asymmetric confidence intervals
- `MnScan` — 1D parameter scans with minimum tracking
- `MnContours` — 2D confidence contour computation
- Global correlation coefficients and covariance squeeze utilities
- 76 tests (48 unit + 8 Migrad + 6 Simplex + 4 Hesse + 3 Minos + 3 Scan + 2 Contours + 2 doctests)

## [0.2.0] - 2026-02-09

### Added

- `MnMigrad` — quasi-Newton minimizer with DFP/BFGS inverse Hessian update
- Quadratic interpolation (`MnParabola`, `MnParabolaPoint`)
- Parabolic line search (`mn_linesearch`)
- Positive-definite matrix correction (`make_pos_def`)
- Two-point central difference gradient (`Numerical2PGradientCalculator`)
- Full Migrad implementation: seed, builder, minimizer, public API
- 58 tests (42 unit + 8 Migrad + 6 Simplex + 2 doctests)

## [0.1.0] - 2026-02-09

### Added

- `FCN` trait with blanket impl for closures
- Parameter management: free, bounded, fixed, constant
- Sin/sqrt parameter transforms for bounded optimization
- `MnSimplex` minimizer (Minuit variant of Nelder-Mead)
- `FunctionMinimum` result type with Display impl
- `MnStrategy` (low/medium/high optimization presets)
- `nalgebra` for all linear algebra

[Unreleased]: https://github.com/machine-made/minuit2-rs/compare/v0.5.5...HEAD
[0.5.5]: https://github.com/machine-made/minuit2-rs/compare/v0.5.4...v0.5.5
[0.5.4]: https://github.com/machine-made/minuit2-rs/compare/v0.5.3...v0.5.4
[0.5.3]: https://github.com/machine-made/minuit2-rs/compare/v0.5.2...v0.5.3
[0.5.2]: https://github.com/machine-made/minuit2-rs/compare/v0.5.1...v0.5.2
[0.5.1]: https://github.com/machine-made/minuit2-rs/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/machine-made/minuit2-rs/compare/v0.4.2...v0.5.0
[0.4.2]: https://github.com/machine-made/minuit2-rs/compare/v0.4.1...v0.4.2
[0.4.1]: https://github.com/machine-made/minuit2-rs/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/machine-made/minuit2-rs/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/machine-made/minuit2-rs/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/machine-made/minuit2-rs/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/machine-made/minuit2-rs/commits/v0.1.0
