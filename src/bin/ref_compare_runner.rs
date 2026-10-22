use minuit2::{
    FCN, FCNGradient, MnContours, MnHesse, MnMigrad, MnMinimize, MnMinos, MnScan, MnSimplex,
};

#[derive(Debug, Clone)]
struct RunResult {
    workload: String,
    algorithm: String,
    valid: bool,
    fval: f64,
    edm: f64,
    nfcn: usize,
    params: Vec<f64>,
    errors: Vec<f64>,
    covariance: Option<Vec<Vec<f64>>>,
    minos: Option<MinosResult>,
}

#[derive(Debug, Clone)]
struct MinosResult {
    valid: bool,
    parameter: usize,
    lower: f64,
    upper: f64,
}

fn parse_workload_arg() -> Option<String> {
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == "--workload" {
            return args.next();
        }
        if let Some(rest) = arg.strip_prefix("--workload=") {
            return Some(rest.to_string());
        }
    }
    None
}

fn vec_json(v: &[f64]) -> String {
    let items = v
        .iter()
        .map(|x| format!("{x:.17}"))
        .collect::<Vec<_>>()
        .join(",");
    format!("[{items}]")
}

fn mat_json(m: &[Vec<f64>]) -> String {
    let rows = m.iter().map(|r| vec_json(r)).collect::<Vec<_>>().join(",");
    format!("[{rows}]")
}

fn to_json(result: &RunResult) -> String {
    let cov_json = result
        .covariance
        .as_ref()
        .map_or_else(|| "null".to_string(), |m| mat_json(m));
    let minos_json = result.minos.as_ref().map_or_else(
        || "null".to_string(),
        |m| {
            format!(
                "{{\"valid\":{},\"parameter\":{},\"lower\":{:.17},\"upper\":{:.17}}}",
                if m.valid { "true" } else { "false" },
                m.parameter,
                m.lower,
                m.upper
            )
        },
    );

    format!(
        concat!(
            "{{",
            "\"runner\":\"minuit2-rs\",",
            "\"workload\":\"{}\",",
            "\"algorithm\":\"{}\",",
            "\"valid\":{},",
            "\"fval\":{:.17},",
            "\"edm\":{:.17},",
            "\"nfcn\":{},",
            "\"params\":{},",
            "\"errors\":{},",
            "\"has_covariance\":{},",
            "\"covariance\":{},",
            "\"has_minos\":{},",
            "\"minos\":{}",
            "}}"
        ),
        result.workload,
        result.algorithm,
        if result.valid { "true" } else { "false" },
        result.fval,
        result.edm,
        result.nfcn,
        vec_json(&result.params),
        vec_json(&result.errors),
        if result.covariance.is_some() {
            "true"
        } else {
            "false"
        },
        cov_json,
        if result.minos.is_some() {
            "true"
        } else {
            "false"
        },
        minos_json
    )
}

fn covariance_dense(min: &minuit2::FunctionMinimum) -> Option<Vec<Vec<f64>>> {
    let cov = min.user_state().covariance()?;
    let n = cov.nrow();
    let mut out = vec![vec![0.0; n]; n];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            *value = cov.get(i, j);
        }
    }
    Some(out)
}

fn common_result(workload: &str, algorithm: &str, min: &minuit2::FunctionMinimum) -> RunResult {
    RunResult {
        workload: workload.to_string(),
        algorithm: algorithm.to_string(),
        valid: min.is_valid(),
        fval: min.fval(),
        edm: min.edm(),
        nfcn: min.nfcn(),
        params: min.params(),
        errors: min.user_state().errors(),
        covariance: covariance_dense(min),
        minos: None,
    }
}

struct Quadratic3;

impl FCN for Quadratic3 {
    fn value(&self, p: &[f64]) -> f64 {
        let (x, y, z) = (p[0], p[1], p[2]);
        x * x + 10.0 * y * y + 100.0 * z * z + 2.0 * x * y + 4.0 * x * z + 8.0 * y * z
    }
}

impl FCNGradient for Quadratic3 {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        let (x, y, z) = (p[0], p[1], p[2]);
        vec![
            2.0 * x + 2.0 * y + 4.0 * z,
            2.0 * x + 20.0 * y + 8.0 * z,
            4.0 * x + 8.0 * y + 200.0 * z,
        ]
    }
}

struct Rosenbrock2;

impl FCN for Rosenbrock2 {
    fn value(&self, p: &[f64]) -> f64 {
        let (x, y) = (p[0], p[1]);
        let t1 = y - x * x;
        let t2 = 1.0 - x;
        100.0 * t1 * t1 + t2 * t2
    }
}

struct Quadratic2;

impl FCN for Quadratic2 {
    fn value(&self, p: &[f64]) -> f64 {
        let (x, y) = (p[0], p[1]);
        let dx = x - 1.0;
        let dy = y + 2.0;
        dx * dx + 4.0 * dy * dy + 0.3 * x * y
    }
}

struct QuadraticNoG2;

impl FCN for QuadraticNoG2 {
    fn value(&self, p: &[f64]) -> f64 {
        let (x, y) = (p[0], p[1]);
        let dx = x - 1.0;
        let dy = y + 2.0;
        dx * dx + dy * dy
    }

    fn has_hessian(&self) -> bool {
        true
    }

    fn hessian(&self, _p: &[f64]) -> Vec<f64> {
        vec![2.0, 0.0, 2.0]
    }

    fn has_g2(&self) -> bool {
        false
    }
}

impl FCNGradient for QuadraticNoG2 {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        vec![2.0 * (p[0] - 1.0), 2.0 * (p[1] + 2.0)]
    }
}

fn run_quadratic3_fixx_migrad() -> RunResult {
    let f = Quadratic3;
    let min = MnMigrad::new()
        .add("x", 1.0, 0.1)
        .add("y", 2.0, 0.1)
        .add("z", 3.0, 0.1)
        .fix(0)
        .tolerance(0.1)
        .minimize_grad(&f);
    common_result("quadratic3_fixx_migrad", "migrad", &min)
}

fn run_quadratic3_fixx_hesse() -> RunResult {
    let f = Quadratic3;
    let min = MnMigrad::new()
        .add("x", 1.0, 0.1)
        .add("y", 2.0, 0.1)
        .add("z", 3.0, 0.1)
        .fix(0)
        .tolerance(0.1)
        .minimize_grad(&f);
    let min_h = MnHesse::new().calculate(&f, &min);
    common_result("quadratic3_fixx_hesse", "migrad+hesse", &min_h)
}

fn run_rosenbrock2_migrad() -> RunResult {
    let f = Rosenbrock2;
    let min = MnMigrad::new()
        .add("x", 0.0, 0.1)
        .add("y", 0.0, 0.1)
        .tolerance(0.1)
        .minimize(&f);
    common_result("rosenbrock2_migrad", "migrad", &min)
}

fn run_quadratic2_minos_p0() -> RunResult {
    let f = Quadratic2;
    let min = MnMigrad::new()
        .add("x", 0.4, 0.1)
        .add("y", -1.0, 0.1)
        .tolerance(0.1)
        .minimize(&f);
    let min_h = MnHesse::new().calculate(&f, &min);
    let me = MnMinos::new(&f, &min_h).minos_error(0);

    let mut result = common_result("quadratic2_minos_p0", "migrad+hesse+minos", &min_h);
    result.minos = Some(MinosResult {
        valid: me.is_valid(),
        parameter: me.parameter(),
        lower: me.lower_error(),
        upper: me.upper_error(),
    });
    result
}

fn run_quadratic2_minos_p1() -> RunResult {
    let f = Quadratic2;
    let min = MnMigrad::new()
        .add("x", 0.4, 0.1)
        .add("y", -1.0, 0.1)
        .tolerance(0.1)
        .minimize(&f);
    let min_h = MnHesse::new().calculate(&f, &min);
    let me = MnMinos::new(&f, &min_h).minos_error(1);

    let mut result = common_result("quadratic2_minos_p1", "migrad+hesse+minos", &min_h);
    result.minos = Some(MinosResult {
        valid: me.is_valid(),
        parameter: me.parameter(),
        lower: me.lower_error(),
        upper: me.upper_error(),
    });
    result
}

fn run_quadratic2_simplex() -> RunResult {
    let f = Quadratic2;
    let min = MnSimplex::new()
        .add("x", 0.4, 0.1)
        .add("y", -1.0, 0.1)
        .tolerance(0.1)
        .minimize(&f);
    common_result("quadratic2_simplex", "simplex", &min)
}

fn run_rosenbrock2_minimize() -> RunResult {
    let f = Rosenbrock2;
    let min = MnMinimize::new()
        .add("x", 0.0, 0.1)
        .add("y", 0.0, 0.1)
        .tolerance(0.1)
        .minimize(&f);
    common_result("rosenbrock2_minimize", "minimize", &min)
}

fn run_quadratic2_limited_migrad() -> RunResult {
    let f = Quadratic2;
    let min = MnMigrad::new()
        .add_limited("x", 0.4, 0.1, 0.0, 2.0)
        .add_limited("y", -1.0, 0.1, -3.0, -1.0)
        .tolerance(0.1)
        .minimize(&f);
    common_result("quadratic2_limited_migrad", "migrad", &min)
}

fn run_quadratic2_lower_limited_migrad() -> RunResult {
    let f = Quadratic2;
    let min = MnMigrad::new()
        .add_lower_limited("x", 0.4, 0.1, 0.0)
        .add_lower_limited("y", -1.0, 0.1, -2.5)
        .tolerance(0.1)
        .minimize(&f);
    common_result("quadratic2_lower_limited_migrad", "migrad", &min)
}

fn run_quadratic2_upper_limited_migrad() -> RunResult {
    let f = Quadratic2;
    let min = MnMigrad::new()
        .add_upper_limited("x", 0.4, 0.1, 1.8)
        .add_upper_limited("y", -1.0, 0.1, -1.5)
        .tolerance(0.1)
        .minimize(&f);
    common_result("quadratic2_upper_limited_migrad", "migrad", &min)
}

fn run_rosenbrock2_migrad_strategy2() -> RunResult {
    let f = Rosenbrock2;
    let min = MnMigrad::new()
        .add("x", 0.0, 0.1)
        .add("y", 0.0, 0.1)
        .with_strategy(2)
        .tolerance(0.1)
        .minimize(&f);
    common_result("rosenbrock2_migrad_strategy2", "migrad_s2", &min)
}

fn run_quadratic2_scan_p0() -> RunResult {
    let f = Quadratic2;
    let min = MnMigrad::new()
        .add("x", 0.4, 0.1)
        .add("y", -1.0, 0.1)
        .tolerance(0.1)
        .minimize(&f);
    let _points = MnScan::new(&f, &min).scan(0, 61, 0.0, 0.0);
    common_result("quadratic2_scan_p0", "migrad+scan", &min)
}

fn run_quadratic2_scan_p1_limited() -> RunResult {
    let f = Quadratic2;
    let min = MnMigrad::new()
        .add_limited("x", 0.4, 0.1, 0.0, 2.0)
        .add_limited("y", -1.0, 0.1, -3.0, -1.0)
        .tolerance(0.1)
        .minimize(&f);
    let _points = MnScan::new(&f, &min).scan(1, 61, 0.0, 0.0);
    common_result("quadratic2_scan_p1_limited", "migrad+scan", &min)
}

fn run_quadratic2_contours_01() -> RunResult {
    let f = Quadratic2;
    let min = MnMigrad::new()
        .add("x", 0.4, 0.1)
        .add("y", -1.0, 0.1)
        .tolerance(0.1)
        .minimize(&f);
    let min_h = MnHesse::new().calculate(&f, &min);
    let _points = MnContours::new(&f, &min_h).points(0, 1, 12);
    common_result("quadratic2_contours_01", "migrad+hesse+contours", &min_h)
}

fn run_quadratic2_no_g2_migrad() -> RunResult {
    let f = QuadraticNoG2;
    let min = MnMigrad::new()
        .add("x", 0.4, 0.1)
        .add("y", -1.0, 0.1)
        .tolerance(0.1)
        .minimize_grad(&f);
    common_result("quadratic2_no_g2_migrad", "migrad_no_g2", &min)
}

fn main() {
    let Some(workload) = parse_workload_arg() else {
        eprintln!("usage: cargo run --bin ref_compare_runner -- --workload <id>");
        std::process::exit(2);
    };

    let result = match workload.as_str() {
        "quadratic3_fixx_migrad" => run_quadratic3_fixx_migrad(),
        "quadratic3_fixx_hesse" => run_quadratic3_fixx_hesse(),
        "rosenbrock2_migrad" => run_rosenbrock2_migrad(),
        "quadratic2_minos_p0" => run_quadratic2_minos_p0(),
        "quadratic2_minos_p1" => run_quadratic2_minos_p1(),
        "quadratic2_simplex" => run_quadratic2_simplex(),
        "rosenbrock2_minimize" => run_rosenbrock2_minimize(),
        "quadratic2_limited_migrad" => run_quadratic2_limited_migrad(),
        "quadratic2_lower_limited_migrad" => run_quadratic2_lower_limited_migrad(),
        "quadratic2_upper_limited_migrad" => run_quadratic2_upper_limited_migrad(),
        "rosenbrock2_migrad_strategy2" => run_rosenbrock2_migrad_strategy2(),
        "quadratic2_scan_p0" => run_quadratic2_scan_p0(),
        "quadratic2_scan_p1_limited" => run_quadratic2_scan_p1_limited(),
        "quadratic2_contours_01" => run_quadratic2_contours_01(),
        "quadratic2_no_g2_migrad" => run_quadratic2_no_g2_migrad(),
        _ => {
            eprintln!("unknown workload: {workload}");
            std::process::exit(3);
        }
    };

    println!("{}", to_json(&result));
}
