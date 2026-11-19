# Contributing

Keep changes focused and preserve the public API unless an API change has been agreed.
Explain changes to numerical results and keep existing tolerances and waivers.

Run these checks before submitting a pull request:

```sh
cargo fmt --all -- --check
cargo clippy --locked --all-targets --all-features -- -D warnings
cargo test --locked --all-features
cargo test --locked --no-default-features
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --all-features
cargo deny check
cargo test --locked --features parallel --test parallel_determinism
python3 -I scripts/check_rust_baseline.py
```

Run the selected house checks from outside the repository, using a local copy of the standards audit:

```sh
python3 -I path/to/standards/scripts/audit.py path/to/minuit2-rs --only changelog,files,deny,rust-toolchain,actions-pinned,dependabot
```

Every numerical change needs a differential run against the ROOT reference runner:

```sh
python3 -I scripts/compare_ref_vs_rust.py --report-dir "$PWD/target/verification"
```

Build the reference runner with `scripts/build_root_reference_runner.sh` before running that comparison.
Review every changed field, including function call counts.
Keep the committed workload results unchanged when checking lint or documentation edits.

The criterion benchmarks no longer run on GitHub Actions. To compare your branch with main, run this on a quiet machine:
`bash scripts/bench_x3.sh origin/main HEAD`
