use minuit2::MnUserParameters;
use nalgebra::DMatrix;

#[derive(Clone, Copy)]
enum BoundType {
    Unbounded,
    Upper,
    Lower,
    Double,
}

/// Reference behavior checked against ROOT Minuit2 `math/minuit2/test/testCovariance.cxx`:
/// transforming covariance internal->external should preserve a known
/// external covariance for representative bound configurations.
fn run_covariance_roundtrip(bound_type: BoundType) {
    let mut upar = MnUserParameters::new();
    upar.add("x", 1.0, 0.1);
    upar.add("y", 1.0, 0.1);
    upar.add("z", 1.0, 0.1);
    upar.add("x0", 2.0, 0.1);
    upar.add("y0", 2.0, 0.1);
    upar.add("z0", 2.0, 0.1);

    match bound_type {
        BoundType::Upper => {
            upar.set_upper_limit(0, 5.0);
            upar.set_upper_limit(4, 5.0);
        }
        BoundType::Lower => {
            upar.set_lower_limit(0, -5.0);
            upar.set_lower_limit(4, -5.0);
        }
        BoundType::Double => {
            upar.set_limits(0, -5.0, 5.0);
            upar.set_limits(4, -5.0, 5.0);
        }
        BoundType::Unbounded => {}
    }

    let trafo = upar.trafo();
    let n = trafo.variable_parameters();
    assert_eq!(n, 6);

    // Reference fixture uses diag=2 and first off-diagonal=-1.
    let mut ext_cov = DMatrix::<f64>::zeros(n, n);
    for i in 0..n {
        ext_cov[(i, i)] = 2.0;
        if i + 1 < n {
            ext_cov[(i, i + 1)] = -1.0;
            ext_cov[(i + 1, i)] = -1.0;
        }
    }

    // Build internal covariance that should map back exactly to ext_cov.
    let internal = trafo.initial_internal_values();
    let jac: Vec<f64> = (0..n)
        .map(|i| {
            let ext = trafo.ext_of_int(i);
            trafo.dint2ext(ext, internal[i])
        })
        .collect();

    let mut int_cov = DMatrix::<f64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            int_cov[(i, j)] = ext_cov[(i, j)] / (jac[i] * jac[j]);
        }
    }

    let ext_back = trafo.int2ext_covariance(&internal, &int_cov);
    for i in 0..n {
        for j in i..n {
            let got = ext_back.get(i, j);
            let want = ext_cov[(i, j)];
            assert!(
                (got - want).abs() <= 1e-10,
                "bound={:?} cov({}, {}) expected {}, got {}",
                bound_name(bound_type),
                i,
                j,
                want,
                got
            );
        }
    }
}

fn bound_name(bound_type: BoundType) -> &'static str {
    match bound_type {
        BoundType::Unbounded => "unbounded",
        BoundType::Upper => "upper",
        BoundType::Lower => "lower",
        BoundType::Double => "double",
    }
}

#[test]
fn root_covariance_unbounded() {
    run_covariance_roundtrip(BoundType::Unbounded);
}

#[test]
fn root_covariance_upper() {
    run_covariance_roundtrip(BoundType::Upper);
}

#[test]
fn root_covariance_lower() {
    run_covariance_roundtrip(BoundType::Lower);
}

#[test]
fn root_covariance_double() {
    run_covariance_roundtrip(BoundType::Double);
}
