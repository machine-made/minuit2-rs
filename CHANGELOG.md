# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

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

[Unreleased]: https://github.com/ricardofrantz/minuit2-rs/compare/v0.5.4...HEAD
[0.5.4]: https://github.com/ricardofrantz/minuit2-rs/compare/v0.5.3...v0.5.4
[0.5.3]: https://github.com/ricardofrantz/minuit2-rs/compare/v0.5.2...v0.5.3
[0.5.2]: https://github.com/ricardofrantz/minuit2-rs/compare/v0.5.1...v0.5.2
[0.5.1]: https://github.com/ricardofrantz/minuit2-rs/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/ricardofrantz/minuit2-rs/compare/v0.4.2...v0.5.0
[0.4.2]: https://github.com/ricardofrantz/minuit2-rs/compare/v0.4.1...v0.4.2
[0.4.1]: https://github.com/ricardofrantz/minuit2-rs/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/ricardofrantz/minuit2-rs/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/ricardofrantz/minuit2-rs/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/ricardofrantz/minuit2-rs/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/ricardofrantz/minuit2-rs/commits/v0.1.0
