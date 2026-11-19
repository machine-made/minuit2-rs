// Rust counterpart of nfcn_accounting.cpp; reuses the runner's objectives.
// Build the library with cargo build --locked --lib --message-format=json.
// Compile this file with rustc --edition=2024, --extern minuit2=<rlib from
// that build's compiler-artifact>, and -L dependency=<each dependency directory>.
// Run the resulting executable to print reported and user-side call counts.
// Pass --trace for every Rosenbrock evaluation; M2_TRACE_STACKS adds stacks.
#[allow(dead_code)]
mod runner {
    include!("../../src/bin/ref_compare_runner.rs");

    use std::cell::Cell;

    struct Counted<T> {
        objective: T,
        calls: Cell<usize>,
        hessian_calls: Cell<usize>,
        trace: bool,
    }

    impl<T: FCN> FCN for Counted<T> {
        fn value(&self, p: &[f64]) -> f64 {
            self.calls.set(self.calls.get() + 1);
            if self.trace {
                println!("TRACE {} {:.16e} {:.16e}", self.calls.get(), p[0], p[1]);
                if std::env::var_os("M2_TRACE_STACKS").is_some() {
                    eprintln!(
                        "CALL {}\n{}",
                        self.calls.get(),
                        std::backtrace::Backtrace::force_capture()
                    );
                }
            }
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

    pub fn counts() {
        let quadratic = Counted {
            objective: QuadraticNoG2,
            calls: Cell::new(0),
            hessian_calls: Cell::new(0),
            trace: false,
        };
        let minimum = MnMigrad::new()
            .add("x", 0.4, 0.1)
            .add("y", -1.0, 0.1)
            .tolerance(0.1)
            .minimize_grad(&quadratic);
        println!(
            "{{\"workload\":\"quadratic2_no_g2_migrad\",\"reported_nfcn\":{},\"user_calls\":{},\"hessian_calls\":{},\"state_nfcn\":{:?}}}",
            minimum.nfcn(),
            quadratic.calls.get(),
            quadratic.hessian_calls.get(),
            minimum
                .states()
                .iter()
                .map(|state| state.nfcn())
                .collect::<Vec<_>>()
        );

        let rosenbrock = Counted {
            objective: Rosenbrock2,
            calls: Cell::new(0),
            hessian_calls: Cell::new(0),
            trace: std::env::args().any(|arg| arg == "--trace"),
        };
        let minimum = MnMigrad::new()
            .add("x", 0.0, 0.1)
            .add("y", 0.0, 0.1)
            .max_fcn(140)
            .tolerance(0.1)
            .minimize(&rosenbrock);
        println!(
            "{{\"workload\":\"rosenbrock2_maxfcn140\",\"reported_nfcn\":{},\"user_calls\":{},\"valid\":{},\"call_limit\":{},\"above_edm\":{},\"edm\":{:.17e}}}",
            minimum.nfcn(),
            rosenbrock.calls.get(),
            minimum.is_valid(),
            minimum.reached_call_limit(),
            minimum.is_above_max_edm(),
            minimum.edm()
        );
    }
}

fn main() {
    runner::counts();
}
