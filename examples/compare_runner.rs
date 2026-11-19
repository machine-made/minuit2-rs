//! Count calls inside the objective for a JSON comparison case.

use minuit2::{FCN, FCNGradient, FunctionMinimum, MnHesse, MnMigrad};
use serde_json::{Value, json};
use std::cell::Cell;
use std::io::{self, Read};

struct Objective {
    case: Value,
    values: Cell<usize>,
    gradients: Cell<usize>,
}

fn number(value: &Value) -> f64 {
    value.as_f64().expect("expected a finite number")
}

impl FCN for Objective {
    #[expect(
        clippy::many_single_char_names,
        reason = "keeps the existing mathematical variable names"
    )]
    fn value(&self, p: &[f64]) -> f64 {
        self.values.set(self.values.get() + 1);
        match self.case["function"].as_str().expect("missing function") {
            "dense" => {
                let center = self.case["center"].as_array().unwrap();
                let h = self.case["hessian"].as_array().unwrap();
                let mut value = 0.0;
                for i in 0..p.len() {
                    for j in 0..p.len() {
                        value += 0.5
                            * (p[i] - number(&center[i]))
                            * number(&h[i][j])
                            * (p[j] - number(&center[j]));
                    }
                }
                value
            }
            "quadratic2" => (p[0] - 1.0).powi(2) + 4.0 * (p[1] + 2.0).powi(2) + 0.3 * p[0] * p[1],
            "rosenbrock" => 100.0 * (p[1] - p[0] * p[0]).powi(2) + (1.0 - p[0]).powi(2),
            "bound" => (p[0] - number(&self.case["target"])).powi(2),
            "gaussian" | "poisson" => {
                let poisson = self.case["function"] == "poisson";
                self.case["data"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|row| {
                        let x = number(&row[0]);
                        let y = number(&row[1]);
                        let z = (x - p[1]) / p[2];
                        let peak = p[0] * (-0.5 * z * z).exp();
                        if poisson {
                            let mu = peak + p[3] * (-p[4] * (x + 5.0)).exp();
                            if mu <= 0.0 {
                                return f64::INFINITY;
                            }
                            // The offset keeps the minimum near zero.
                            let d = (mu - y) / y;
                            y * (d - d.ln_1p())
                        } else {
                            ((peak - y) / number(&row[2])).powi(2)
                        }
                    })
                    .sum()
            }
            other => panic!("unknown function: {other}"),
        }
    }

    fn error_def(&self) -> f64 {
        self.case["errordef"].as_f64().unwrap_or(1.0)
    }
}

impl FCNGradient for Objective {
    fn gradient(&self, p: &[f64]) -> Vec<f64> {
        self.gradients.set(self.gradients.get() + 1);
        match self.case["function"].as_str().unwrap() {
            "quadratic2" => vec![
                2.0 * (p[0] - 1.0) + 0.3 * p[1],
                8.0 * (p[1] + 2.0) + 0.3 * p[0],
            ],
            "rosenbrock" => vec![
                -400.0 * p[0] * (p[1] - p[0] * p[0]) + 2.0 * (p[0] - 1.0),
                200.0 * (p[1] - p[0] * p[0]),
            ],
            other => panic!("no user gradient for {other}"),
        }
    }
}

fn snapshot(minimum: &FunctionMinimum, objective: &Objective) -> Value {
    let covariance = minimum.user_state().covariance().map(|cov| {
        (0..cov.nrow())
            .map(|i| (0..cov.nrow()).map(|j| cov.get(i, j)).collect::<Vec<_>>())
            .collect::<Vec<_>>()
    });
    json!({
        "valid": minimum.is_valid(), "fval": minimum.fval(), "params": minimum.params(),
        "edm": minimum.edm(), "covariance": covariance,
        "value_calls": objective.values.get(), "gradient_calls": objective.gradients.get(),
    })
}

#[expect(
    clippy::cast_possible_truncation,
    reason = "keeps the existing JSON integer conversions"
)]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let case: Value = serde_json::from_str(&input)?;
    let strategy = case["strategy"].as_u64().unwrap_or(1) as u32;
    if strategy > 2 {
        return Err("Rust strategy is not implemented".into());
    }
    let mut migrad = MnMigrad::new()
        .with_strategy(strategy)
        .tolerance(case["tolerance"].as_f64().unwrap_or(0.1));
    let start = case["start"].as_array().ok_or("missing start")?;
    let steps = case["steps"].as_array().ok_or("missing steps")?;
    for (i, value) in start.iter().enumerate() {
        let name = format!("p{i}");
        let value = number(value);
        let step = number(&steps[i]);
        migrad = match (case["limits"][i][0].as_f64(), case["limits"][i][1].as_f64()) {
            (Some(lo), Some(hi)) => migrad.add_limited(name, value, step, lo, hi),
            (Some(lo), None) => migrad.add_lower_limited(name, value, step, lo),
            (None, Some(hi)) => migrad.add_upper_limited(name, value, step, hi),
            (None, None) => migrad.add(name, value, step),
        };
    }
    if let Some(fixed) = case["fixed"].as_array() {
        for i in fixed {
            migrad = migrad.fix(i.as_u64().ok_or("invalid fixed index")? as usize);
        }
    }
    let objective = Objective {
        case,
        values: Cell::new(0),
        gradients: Cell::new(0),
    };
    let minimum = if objective.case["gradient"].as_bool().unwrap_or(false) {
        migrad.minimize_grad(&objective)
    } else {
        migrad.minimize(&objective)
    };
    let before = snapshot(&minimum, &objective);
    let values = objective.values.get();
    let gradients = objective.gradients.get();
    let hesse = MnHesse::new()
        .with_strategy(strategy)
        .calculate(&objective, &minimum);
    println!(
        "{}",
        json!({"migrad": before, "hesse": snapshot(&hesse, &objective),
        "hesse_value_calls": objective.values.get() - values,
        "hesse_gradient_calls": objective.gradients.get() - gradients})
    );
    Ok(())
}
