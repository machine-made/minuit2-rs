// Rust counterpart of nfcn_accounting.cpp; reuses the runner's objectives.
// From the repository root, build the library, then this file. Cargo keeps the
// dependency libraries in target/debug/build/<crate>/<hash>/out, so the loop in
// the command substitution adds each of those directories to the search path.
// Each run prints one JSON document, which is the stored file:
//   cargo build --locked --lib
//   mkdir -p target/nfcn-accounting-4
//   rustc --edition=2024 -g tools/root_evidence/nfcn_accounting.rs \
//     --extern minuit2=target/debug/libminuit2.rlib -L dependency=target/debug/deps \
//     $(for d in target/debug/build/*/*/out; do printf ' -L dependency=%s' $d; done) \
//     -o target/nfcn-accounting-4/nfcn_accounting_rust
//   target/nfcn-accounting-4/nfcn_accounting_rust --trace quadratic2_no_g2_migrad \
//     > reports/verification/raw/trace/nfcn_accounting_quadratic2_no_g2_migrad.rust.json
//   target/nfcn-accounting-4/nfcn_accounting_rust --trace rosenbrock2_maxfcn140 \
//     > reports/verification/raw/trace/nfcn_accounting_rosenbrock2_maxfcn140.rust.json
//   target/nfcn-accounting-4/nfcn_accounting_rust --trace quadratic2_minos_before_minos \
//     > reports/verification/raw/trace/nfcn_accounting_quadratic2_minos_before_minos.rust.json
// The document holds "provenance" (the build and run commands above) and "result".
// The program exits with 1 if reported_nfcn and user_calls differ.
// Pass --trace to include every evaluation and phase; M2_TRACE_STACKS adds stacks.
#[allow(dead_code)]
mod phase_counts;

#[allow(dead_code)]
mod runner {
    include!("../../src/bin/ref_compare_runner.rs");

    use super::phase_counts::{PhaseCounts, begin_record, end_record, field, provenance};
    use std::cell::{Cell, RefCell};

    // Commands that build this program; the run command is added to the provenance.
    const PROGRAM: &str = "tools/root_evidence/nfcn_accounting.rs";
    const BUILD_COMMAND: &str = concat!(
        "mkdir -p target/nfcn-accounting-4 && rustc --edition=2024 -g tools/root_evidence/nfcn_accounting.rs ",
        "--extern minuit2=target/debug/libminuit2.rlib -L dependency=target/debug/deps ",
        "$(for d in target/debug/build/*/*/out; do printf ' -L dependency=%s' $d; done) ",
        "-o target/nfcn-accounting-4/nfcn_accounting_rust"
    );

    struct Counted<T> {
        objective: T,
        counts: RefCell<PhaseCounts>,
        hessian_calls: Cell<usize>,
    }

    impl<T: FCN> FCN for Counted<T> {
        fn value(&self, p: &[f64]) -> f64 {
            self.counts.borrow_mut().record(p);
            self.objective.value(p)
        }

        fn has_g2(&self) -> bool {
            self.objective.has_g2()
        }

        fn g2(&self, p: &[f64]) -> Vec<f64> {
            self.objective.g2(p)
        }

        fn has_hessian(&self) -> bool {
            self.objective.has_hessian()
        }

        fn hessian(&self, p: &[f64]) -> Vec<f64> {
            self.hessian_calls.set(self.hessian_calls.get() + 1);
            self.objective.hessian(p)
        }
    }

    impl<T: FCNGradient> FCNGradient for Counted<T> {
        fn gradient(&self, p: &[f64]) -> Vec<f64> {
            self.objective.gradient(p)
        }
    }

    pub fn counts() -> Result<(), String> {
        let usage = "usage: nfcn_accounting_rust [--trace] \
            quadratic2_no_g2_migrad|rosenbrock2_maxfcn140|quadratic2_minos_before_minos";
        let mut trace = false;
        let mut workload = None;
        let mut command = String::from("target/nfcn-accounting-4/nfcn_accounting_rust");
        for arg in std::env::args().skip(1) {
            command = format!("{command} {arg}");
            if arg == "--trace" {
                trace = true;
            } else if arg.starts_with("--") || workload.is_some() {
                return Err(format!("rejected argument {arg}\n{usage}"));
            } else {
                workload = Some(arg);
            }
        }
        match workload.as_deref() {
            Some("quadratic2_no_g2_migrad") => quadratic2_no_g2_migrad(trace, &command),
            Some("rosenbrock2_maxfcn140") => rosenbrock2_maxfcn140(trace, &command),
            Some("quadratic2_minos_before_minos") => quadratic2_minos_before_minos(trace, &command),
            _ => Err(usage.to_owned()),
        }
    }

    fn quadratic2_no_g2_migrad(trace: bool, command: &str) -> Result<(), String> {
        let quadratic = Counted {
            objective: QuadraticNoG2,
            counts: RefCell::new(PhaseCounts::new(trace)),
            hessian_calls: Cell::new(0),
        };
        let minimum = MnMigrad::new()
            .add("x", 0.4, 0.1)
            .add("y", -1.0, 0.1)
            .tolerance(0.1)
            .minimize_grad(&quadratic);
        begin_record(
            &provenance(PROGRAM, BUILD_COMMAND, command),
            "quadratic2_no_g2_migrad",
        );
        field("reported_nfcn", minimum.nfcn());
        field("hessian_calls", quadratic.hessian_calls.get());
        print!(",\n    \"state_nfcn\": [");
        for (i, state) in minimum.states().iter().enumerate() {
            print!("{}\n      {}", if i > 0 { "," } else { "" }, state.nfcn());
        }
        print!("\n    ]");
        quadratic.counts.borrow().print();
        end_record();
        quadratic.counts.borrow().check_reported(minimum.nfcn())
    }

    fn rosenbrock2_maxfcn140(trace: bool, command: &str) -> Result<(), String> {
        let rosenbrock = Counted {
            objective: Rosenbrock2,
            counts: RefCell::new(PhaseCounts::new(trace)),
            hessian_calls: Cell::new(0),
        };
        let minimum = MnMigrad::new()
            .add("x", 0.0, 0.1)
            .add("y", 0.0, 0.1)
            .max_fcn(140)
            .tolerance(0.1)
            .minimize(&rosenbrock);
        begin_record(
            &provenance(PROGRAM, BUILD_COMMAND, command),
            "rosenbrock2_maxfcn140",
        );
        field("reported_nfcn", minimum.nfcn());
        field("valid", minimum.is_valid());
        field("call_limit", minimum.reached_call_limit());
        field("above_edm", minimum.is_above_max_edm());
        field("edm", minimum.edm());
        rosenbrock.counts.borrow().print();
        end_record();
        rosenbrock.counts.borrow().check_reported(minimum.nfcn())
    }

    fn quadratic2_minos_before_minos(trace: bool, command: &str) -> Result<(), String> {
        let minos = Counted {
            objective: Quadratic2,
            counts: RefCell::new(PhaseCounts::new(trace)),
            hessian_calls: Cell::new(0),
        };
        let minimum = MnMigrad::new()
            .add("x", 0.4, 0.1)
            .add("y", -1.0, 0.1)
            .tolerance(0.1)
            .minimize(&minos);
        let migrad_calls = minimum.nfcn();
        minos.counts.borrow_mut().begin_standalone();
        let minimum = MnHesse::new().calculate(&minos, &minimum);
        begin_record(
            &provenance(PROGRAM, BUILD_COMMAND, command),
            "quadratic2_minos_before_minos",
        );
        field("migrad_nfcn", migrad_calls);
        field("reported_nfcn", minimum.nfcn());
        field("hessian_calls", minos.hessian_calls.get());
        minos.counts.borrow().print();
        end_record();
        minos.counts.borrow().check_reported(minimum.nfcn())
    }
}

fn main() -> std::process::ExitCode {
    match runner::counts() {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("{message}");
            std::process::ExitCode::FAILURE
        }
    }
}
