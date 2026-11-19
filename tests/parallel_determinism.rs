#![cfg(feature = "parallel")]

use minuit2::{MnMigrad, MnScan};
use rayon::ThreadPoolBuilder;

fn quadratic2(p: &[f64]) -> f64 {
    let x = p[0] - 1.0;
    let y = p[1] + 2.0;
    x * x + 4.0 * y * y + 0.3 * p[0] * p[1]
}

#[test]
fn scan_bits_equal_for_1_and_4_threads() {
    let minimum = MnMigrad::new()
        .add("x", 0.0, 0.1)
        .add("y", 0.0, 0.1)
        .minimize(&quadratic2);
    assert!(minimum.is_valid());
    let scan = MnScan::new(&quadratic2, &minimum);
    let serial = scan.scan_serial(0, 100, -3.0, 3.0);
    assert_eq!(serial.len(), 101);
    let serial_bits: Vec<_> = serial
        .iter()
        .map(|(x, f)| (x.to_bits(), f.to_bits()))
        .collect();

    for threads in [1, 4] {
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("thread pool must start");
        let parallel = pool.install(|| scan.scan_parallel(0, 100, -3.0, 3.0));
        let parallel_bits: Vec<_> = parallel
            .iter()
            .map(|(x, f)| (x.to_bits(), f.to_bits()))
            .collect();
        assert_eq!(serial_bits, parallel_bits, "thread count {threads}");
    }
}
