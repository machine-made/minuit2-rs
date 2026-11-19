// Rust counterpart of nfcn_accounting.cpp; reuses the runner's objectives.
// Build the library with cargo build --locked --lib --message-format=json.
// Compile this file with rustc --edition=2024, --extern minuit2=<rlib from
// that build's compiler-artifact>, and -L dependency=<each dependency directory>.
// Run the resulting executable to print reported and user-side call counts.
// Pass --trace to include every evaluation and phase; M2_TRACE_STACKS adds stacks.
#[allow(dead_code)]
mod phase_counts;

#[allow(dead_code)]
mod runner {
    include!("../../src/bin/ref_compare_runner.rs");

    use super::phase_counts::PhaseCounts;
    use std::cell::{Cell, RefCell};

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

    pub fn counts() {
        let trace = std::env::args().any(|arg| arg == "--trace");
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
        print!(
            "{{\"workload\":\"quadratic2_no_g2_migrad\",\"reported_nfcn\":{},\"hessian_calls\":{},\"state_nfcn\":{:?}",
            minimum.nfcn(),
            quadratic.hessian_calls.get(),
            minimum
                .states()
                .iter()
                .map(|state| state.nfcn())
                .collect::<Vec<_>>()
        );

        quadratic.counts.borrow().print();
        println!("}}");

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
        print!(
            "{{\"workload\":\"rosenbrock2_maxfcn140\",\"reported_nfcn\":{},\"valid\":{},\"call_limit\":{},\"above_edm\":{},\"edm\":{:.17e}",
            minimum.nfcn(),
            minimum.is_valid(),
            minimum.reached_call_limit(),
            minimum.is_above_max_edm(),
            minimum.edm()
        );
        rosenbrock.counts.borrow().print();
        println!("}}");

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
        print!(
            "{{\"workload\":\"quadratic2_minos_before_minos\",\"migrad_nfcn\":{},\"reported_nfcn\":{},\"hessian_calls\":{}",
            migrad_calls,
            minimum.nfcn(),
            minos.hessian_calls.get()
        );
        minos.counts.borrow().print();
        println!("}}");
    }
}

fn main() {
    runner::counts();
}
