// Hahn1 phase counts using the parity test's exact objective and NIST loader.
// Build as described in nfcn_accounting.rs, also setting CARGO_MANIFEST_DIR to
// the repository root for the included test helpers. Run from the repository:
//   target/nfcn-accounting-4/hahn1_phase_counts_rust --trace
// The optional dataset argument defaults to examples/data/nist/Hahn1.dat.
// Phase classification is in phase_counts.rs; M2_TRACE_STACKS prints stacks.
#[allow(dead_code)]
mod phase_counts;
#[allow(dead_code)]
mod hahn1 {
    include!("../../tests/hahn1_core_parity.rs");

    pub fn counts() {
        let mut path = PathBuf::from("examples/data/nist/Hahn1.dat");
        let mut trace = false;
        for arg in std::env::args().skip(1) {
            if arg == "--trace" {
                trace = true;
            } else {
                path = PathBuf::from(arg);
            }
        }
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
        print!(
            "{{\"workload\":\"hahn1_start2\",\"reported_nfcn\":{},\"valid\":{},\"fval\":{:.17e}",
            minimum.nfcn(),
            minimum.is_valid(),
            minimum.fval()
        );
        counts.borrow().print();
        println!("}}");
    }
}
fn main() {
    hahn1::counts();
}
