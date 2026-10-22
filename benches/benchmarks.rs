use criterion::{Criterion, black_box, criterion_group, criterion_main};
use minuit2::linesearch::mn_linesearch;
use minuit2::minimum::parameters::MinimumParameters;
use minuit2::mn_fcn::MnFcn;
use minuit2::posdef::make_pos_def;
use minuit2::{
    FCN, MinuitParameter, MnContours, MnHesse, MnMachinePrecision, MnMigrad, MnMinimize, MnMinos,
    MnScan, MnSimplex, MnUserTransformation,
};
use nalgebra::{DMatrix, DVector};

struct GaussianChi2 {
    x: Vec<f64>,
    y: Vec<f64>,
}

impl FCN for GaussianChi2 {
    fn value(&self, p: &[f64]) -> f64 {
        let (a, mu, sigma) = (p[0], p[1], p[2]);
        self.x
            .iter()
            .zip(self.y.iter())
            .map(|(&xi, &yi)| {
                let model = a * (-0.5 * ((xi - mu) / sigma).powi(2)).exp();
                (yi - model).powi(2)
            })
            .sum()
    }
}

struct Quadratic1D;

impl FCN for Quadratic1D {
    fn value(&self, p: &[f64]) -> f64 {
        p[0] * p[0]
    }
}

fn bench_rosenbrock_migrad(c: &mut Criterion) {
    let rosenbrock = |p: &[f64]| (1.0 - p[0]).powi(2) + 100.0 * (p[1] - p[0] * p[0]).powi(2);

    c.bench_function("Rosenbrock 2D: MnMigrad minimize", |b| {
        b.iter(|| {
            let result = MnMigrad::new()
                .add("x", 0.0, 0.1)
                .add("y", 0.0, 0.1)
                .minimize(&rosenbrock);
            black_box(result);
        })
    });
}

fn bench_rosenbrock_minimize(c: &mut Criterion) {
    let rosenbrock = |p: &[f64]| (1.0 - p[0]).powi(2) + 100.0 * (p[1] - p[0] * p[0]).powi(2);

    c.bench_function("Rosenbrock 2D: MnMinimize hybrid", |b| {
        b.iter(|| {
            let result = MnMinimize::new()
                .add("x", 0.0, 0.1)
                .add("y", 0.0, 0.1)
                .minimize(&rosenbrock);
            black_box(result);
        })
    });
}

fn bench_rosenbrock_simplex(c: &mut Criterion) {
    let rosenbrock = |p: &[f64]| (1.0 - p[0]).powi(2) + 100.0 * (p[1] - p[0] * p[0]).powi(2);

    c.bench_function("Rosenbrock 2D: MnSimplex minimize", |b| {
        b.iter(|| {
            let result = MnSimplex::new()
                .add("x", 0.0, 0.1)
                .add("y", 0.0, 0.1)
                .minimize(&rosenbrock);
            black_box(result);
        })
    });
}

fn bench_quadratic_4d_migrad(c: &mut Criterion) {
    let quadratic = |p: &[f64]| p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + p[3] * p[3];

    c.bench_function("Quadratic 4D: MnMigrad minimize", |b| {
        b.iter(|| {
            let result = MnMigrad::new()
                .add("x1", 1.0, 0.1)
                .add("x2", 2.0, 0.1)
                .add("x3", 3.0, 0.1)
                .add("x4", 4.0, 0.1)
                .minimize(&quadratic);
            black_box(result);
        })
    });
}

fn bench_quadratic_2d_migrad_hesse(c: &mut Criterion) {
    let quadratic = |p: &[f64]| p[0] * p[0] + p[1] * p[1];

    c.bench_function("Quadratic 2D: MnMigrad + MnHesse", |b| {
        b.iter(|| {
            let migrad_result = MnMigrad::new()
                .add("x", 1.0, 0.1)
                .add("y", 2.0, 0.1)
                .minimize(&quadratic);
            let hesse_result = MnHesse::new().calculate(&quadratic, &migrad_result);
            black_box((migrad_result, hesse_result));
        })
    });
}

fn bench_gaussian_fit_migrad_hesse(c: &mut Criterion) {
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = vec![0.1, 2.3, 7.8, 6.1, 0.2];
    let fcn = GaussianChi2 { x, y };

    c.bench_function("Gaussian fit (chi-square, 3 params): Migrad + Hesse", |b| {
        b.iter(|| {
            let migrad_result = MnMigrad::new()
                .add("A", 8.0, 1.0)
                .add("mu", 4.0, 0.5)
                .add_lower_limited("sigma", 2.0, 0.5, 0.01)
                .minimize(&fcn);
            let hesse_result = MnHesse::new().calculate(&fcn, &migrad_result);
            black_box((migrad_result, hesse_result));
        })
    });
}

fn bench_minos_error(c: &mut Criterion) {
    let quadratic = |p: &[f64]| 2.0 * p[0] * p[0] + 8.0 * p[1] * p[1];
    let minimum = MnMigrad::new()
        .add("x", 5.0, 1.0)
        .add("y", -3.0, 1.0)
        .minimize(&quadratic);
    let hesse = MnHesse::new().calculate(&quadratic, &minimum);

    c.bench_function("Quadratic 2D: MnMinos error(0)", |b| {
        b.iter(|| {
            let me = MnMinos::new(&quadratic, &hesse).minos_error(0);
            black_box(me);
        })
    });
}

fn bench_contours(c: &mut Criterion) {
    let correlated = |p: &[f64]| p[0] * p[0] + p[1] * p[1] + p[0] * p[1];
    let minimum = MnMigrad::new()
        .add("x", 4.0, 1.0)
        .add("y", -2.0, 1.0)
        .minimize(&correlated);
    let hesse = MnHesse::new().calculate(&correlated, &minimum);
    let contours = MnContours::new(&correlated, &hesse);

    c.bench_function("Correlated 2D: MnContours 12 points", |b| {
        b.iter(|| {
            let pts = contours.points(0, 1, 12);
            black_box(pts);
        })
    });
}

fn bench_scan_serial(c: &mut Criterion) {
    let quadratic = |p: &[f64]| p[0] * p[0] + p[1] * p[1];
    let minimum = MnMigrad::new()
        .add("x", 2.0, 0.5)
        .add("y", -1.0, 0.5)
        .minimize(&quadratic);
    let scan = MnScan::new(&quadratic, &minimum);

    c.bench_function("Quadratic 2D: MnScan serial (101 points)", |b| {
        b.iter(|| {
            let pts = scan.scan_serial(0, 100, -3.0, 3.0);
            black_box(pts);
        })
    });
}

#[cfg(feature = "parallel")]
fn bench_scan_parallel(c: &mut Criterion) {
    let quadratic = |p: &[f64]| p[0] * p[0] + p[1] * p[1];
    let minimum = MnMigrad::new()
        .add("x", 2.0, 0.5)
        .add("y", -1.0, 0.5)
        .minimize(&quadratic);
    let scan = MnScan::new(&quadratic, &minimum);

    c.bench_function("Quadratic 2D: MnScan parallel (101 points)", |b| {
        b.iter(|| {
            let pts = scan.scan_parallel(0, 100, -3.0, 3.0);
            black_box(pts);
        })
    });
}

fn bench_kernel_posdef(c: &mut Criterion) {
    let n = 20;
    let mut matrix = DMatrix::<f64>::identity(n, n);
    for i in 0..n {
        matrix[(i, i)] = if i % 3 == 0 {
            -0.1
        } else {
            1.0 + i as f64 * 1e-3
        };
    }
    for i in 0..n {
        for j in 0..n {
            if i != j {
                matrix[(i, j)] = 0.02 * ((i + j) as f64).sin();
            }
        }
    }
    let prec = MnMachinePrecision::new();

    c.bench_function("Kernel: make_pos_def 20x20", |b| {
        b.iter(|| {
            let out = make_pos_def(&matrix, &prec);
            black_box(out);
        })
    });
}

fn bench_kernel_linesearch(c: &mut Criterion) {
    let model = Quadratic1D;
    let trafo = MnUserTransformation::new(vec![MinuitParameter::new(0, "x", 2.0, 0.1)]);
    let fcn = MnFcn::new(&model, &trafo);
    let start = MinimumParameters::new(DVector::from_vec(vec![2.0]), 4.0);
    let step = DVector::from_vec(vec![-1.0]);
    let gdel = step.dot(&DVector::from_vec(vec![4.0]));
    let prec = MnMachinePrecision::new();

    c.bench_function("Kernel: line search 1D quadratic", |b| {
        b.iter(|| {
            let result = mn_linesearch(&fcn, &start, &step, gdel, &prec);
            black_box(result);
        })
    });
}

#[cfg(feature = "parallel")]
criterion_group!(
    benches,
    bench_rosenbrock_migrad,
    bench_rosenbrock_minimize,
    bench_rosenbrock_simplex,
    bench_quadratic_4d_migrad,
    bench_quadratic_2d_migrad_hesse,
    bench_gaussian_fit_migrad_hesse,
    bench_minos_error,
    bench_contours,
    bench_scan_serial,
    bench_scan_parallel,
    bench_kernel_posdef,
    bench_kernel_linesearch
);

#[cfg(not(feature = "parallel"))]
criterion_group!(
    benches,
    bench_rosenbrock_migrad,
    bench_rosenbrock_minimize,
    bench_rosenbrock_simplex,
    bench_quadratic_4d_migrad,
    bench_quadratic_2d_migrad_hesse,
    bench_gaussian_fit_migrad_hesse,
    bench_minos_error,
    bench_contours,
    bench_scan_serial,
    bench_kernel_posdef,
    bench_kernel_linesearch
);
criterion_main!(benches);
