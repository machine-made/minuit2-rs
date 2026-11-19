// Hahn1 phase counts using the parity test's exact objective and NIST loader.
// From the repository root, build the library, then this file (the included test
// helpers read CARGO_MANIFEST_DIR). Cargo keeps the dependency libraries in
// target/debug/build/<crate>/<hash>/out, so the loop in the command substitution
// adds each of those directories to the search path.
// The run prints one JSON document, which is the stored file:
//   cargo build --locked --lib
//   mkdir -p target/nfcn-accounting-4
//   CARGO_MANIFEST_DIR=. rustc --edition=2024 -g tools/root_evidence/hahn1_phase_counts.rs \
//     --extern minuit2=target/debug/libminuit2.rlib -L dependency=target/debug/deps \
//     $(for d in target/debug/build/*/*/out; do printf ' -L dependency=%s' $d; done) \
//     -o target/nfcn-accounting-4/hahn1_phase_counts_rust
//   target/nfcn-accounting-4/hahn1_phase_counts_rust --trace \
//     > reports/verification/raw/trace/hahn1_phase_counts_hahn1_start2.rust.json
// The optional dataset argument defaults to examples/data/nist/Hahn1.dat.
// Other options are rejected with a usage message. The document holds "provenance"
// (the build and run commands above) and "result". The program exits with 1 if
// reported_nfcn and user_calls differ.
// Phase classification is in phase_counts.rs; M2_TRACE_STACKS prints stacks.
#[allow(dead_code)]
mod phase_counts;
#[allow(dead_code)]
mod hahn1 {
    include!("../../tests/hahn1_core_parity.rs");

    use super::phase_counts::{begin_record, end_record, field, json_number};

    // Commands that produce this program (see the header); the run command is added to it.
    fn provenance(command: &str) -> String {
        const HEAD: &str = concat!(
            r#"{
    "program": "tools/root_evidence/hahn1_phase_counts.rs",
    "root_tag": null,
    "library_build_command": "cargo build --locked --lib",
    "build_command": "CARGO_MANIFEST_DIR=. rustc --edition=2024 -g "#,
            r#"tools/root_evidence/hahn1_phase_counts.rs "#,
            r#"--extern minuit2=target/debug/libminuit2.rlib -L dependency=target/debug/deps "#,
            r#"$(for d in target/debug/build/*/*/out; do printf ' -L dependency=%s' $d; done) "#,
            r#"-o target/nfcn-accounting-4/hahn1_phase_counts_rust",
    "command": ""#
        );
        format!("{HEAD}{command}\"\n  }}")
    }

    pub fn counts() -> Result<(), String> {
        let mut path = None;
        let mut trace = false;
        let mut command = String::from("target/nfcn-accounting-4/hahn1_phase_counts_rust");
        for arg in std::env::args().skip(1) {
            command = format!("{command} {arg}");
            if arg == "--trace" {
                trace = true;
            } else if arg.starts_with("--") || path.is_some() {
                return Err(format!(
                    "rejected argument {arg}\nusage: hahn1_phase_counts_rust [--trace] [DATASET]"
                ));
            } else {
                path = Some(PathBuf::from(arg));
            }
        }
        let path = path.unwrap_or_else(|| PathBuf::from("examples/data/nist/Hahn1.dat"));
        let dataset = load_hahn1(&path);
        let fcn = Hahn1Fcn {
            x: dataset.x,
            y: dataset.y,
        };
        let counts = std::cell::RefCell::new(super::phase_counts::PhaseCounts::new(trace));
        let mut migrad = MnMigrad::new();
        for (i, start) in dataset.start2.into_iter().enumerate() {
            migrad = migrad.add(format!("b{}", i + 1), start, 0.1);
        }
        let minimum = migrad
            .with_strategy(1)
            .tolerance(0.1)
            .max_fcn(1_000_000)
            .minimize(&|p: &[f64]| {
                counts.borrow_mut().record(p);
                fcn.value(p)
            });
        begin_record(&provenance(&command), "hahn1_start2");
        field("reported_nfcn", minimum.nfcn());
        field("valid", minimum.is_valid());
        field("fval", json_number(minimum.fval()));
        counts.borrow().print();
        end_record();
        let user_calls = counts.borrow().calls;
        if minimum.nfcn() != user_calls {
            return Err(format!(
                "reported_nfcn {} differs from user_calls {user_calls}",
                minimum.nfcn()
            ));
        }
        Ok(())
    }
}
fn main() -> std::process::ExitCode {
    match hahn1::counts() {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("{message}");
            std::process::ExitCode::FAILURE
        }
    }
}
