# Benchmark Baseline

- Command (default): `cargo bench --bench benchmarks -- --noplot`
- Command (parallel): `cargo bench --features parallel --bench benchmarks -- --noplot`

## Median Time By Benchmark

| Benchmark | Default Median | Parallel Median | Parallel vs Default |
|---|---:|---:|---:|
| Correlated 2D: MnContours 12 points | 27.1 µs | 26.51 µs | 0.98x |
| Gaussian fit (chi-square, 3 params): Migrad + Hesse | 10.88 µs | 10.88 µs | 1.00x |
| Kernel: line search 1D quadratic | 80.71 ns | 80.75 ns | 1.00x |
| Kernel: make_pos_def 20x20 | 15.64 µs | 14.75 µs | 0.94x |
| Quadratic 2D: MnMigrad + MnHesse | 5.099 µs | 5.173 µs | 1.01x |
| Quadratic 2D: MnMinos error(0) | 10.58 µs | 10.41 µs | 0.98x |
| Quadratic 2D: MnScan parallel (101 points) | - | 28.89 µs | - |
| Quadratic 2D: MnScan serial (101 points) | 2.441 µs | 2.138 µs | 0.88x |
| Quadratic 4D: MnMigrad minimize | 4.287 µs | 4.303 µs | 1.00x |
| Rosenbrock 2D: MnMigrad minimize | 17.87 µs | 17.95 µs | 1.00x |
| Rosenbrock 2D: MnMinimize hybrid | 18.01 µs | 19.25 µs | 1.07x |
| Rosenbrock 2D: MnSimplex minimize | 3.448 µs | 3.384 µs | 0.98x |

## Scan Comparison (`parallel` run)

- Serial median: **2.138 µs**
- Parallel median: **28.89 µs**
- Speedup (`serial / parallel`): **0.07x**
