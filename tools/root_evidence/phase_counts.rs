// Evidence-only stack classification. Compile the library and harness with
// debug information and without optimisation so Minuit2 frames remain visible.
// Centre, diagonal and mixed calls use the same coordinate rule as phase_counts.h.
use std::collections::BTreeMap;

/// JSON has no NaN or infinity: write null for a non-finite number. Otherwise write
/// the shortest digits that read back as the same double.
pub fn json_number(value: f64) -> String {
    if value.is_finite() {
        format!("{value:?}")
    } else {
        "null".to_owned()
    }
}

/// A record is the whole output of one run, in the layout of the stored files:
/// two-space indent, one key per line, a trailing newline. The caller prints the
/// result fields, then `PhaseCounts::print`, then `end_record`.
pub fn begin_record(provenance: &str, workload: &str) {
    print!(
        "{{\n  \"provenance\": {provenance},\n  \"result\": {{\n    \"workload\": \"{workload}\""
    );
}

pub fn field(name: &str, value: impl std::fmt::Display) {
    print!(",\n    \"{name}\": {value}");
}

pub fn end_record() {
    println!("\n  }}\n}}");
}

#[derive(Default)]
pub struct PhaseCounts {
    pub calls: usize,
    pub passes: usize,
    pub before_hesse: usize,
    in_hesse: bool,
    standalone: bool,
    pub trace: bool,
    centre: Vec<f64>,
    phases: BTreeMap<String, usize>,
    evaluations: Vec<String>,
}

impl PhaseCounts {
    pub fn new(trace: bool) -> Self {
        Self {
            trace,
            ..Default::default()
        }
    }
    pub fn begin_standalone(&mut self) {
        self.standalone = true;
        self.in_hesse = false;
    }

    pub fn record(&mut self, p: &[f64]) {
        let stack = std::backtrace::Backtrace::force_capture().to_string();
        let hesse = stack.contains("minuit2::hesse::calculator::calculate");
        let phase = if hesse {
            let suffix = if !self.in_hesse {
                if self.passes == 0 {
                    self.before_hesse = self.calls;
                }
                self.passes += 1;
                self.centre = p.to_vec();
                "centre"
            } else if stack.contains("minuit2::hesse::gradient::") {
                "refinement"
            } else {
                match p.iter().zip(&self.centre).filter(|(x, c)| x != c).count() {
                    1 => "diagonal",
                    2 => "mixed",
                    changed => panic!(
                        "unclassified Hesse call {}: {changed} changed coordinates at p = {p:?}",
                        self.calls + 1
                    ),
                }
            };
            format!(
                "{}hesse_{suffix}",
                if self.standalone { "standalone_" } else { "" }
            )
        } else if self.passes > 0 {
            "continuation".to_owned()
        } else if stack.contains("minuit2::migrad::seed::") {
            "seed".to_owned()
        } else {
            "migrad_iteration".to_owned()
        };
        self.in_hesse = hesse;
        self.calls += 1;
        *self.phases.entry(phase.clone()).or_default() += 1;
        if std::env::var_os("M2_TRACE_STACKS").is_some() {
            eprintln!("CALL {}\n{stack}", self.calls);
        }
        if self.trace {
            let coordinates: Vec<String> = p.iter().map(|&x| json_number(x)).collect();
            self.evaluations.push(format!(
                "{{\"call\":{},\"phase\":\"{phase}\",\"p\":[{}]}}",
                self.calls,
                coordinates.join(",")
            ));
        }
    }

    pub fn print(&self) {
        print!(
            ",\n    \"user_calls\": {},\n    \"hesse_passes\": {},\n    \"before_hesse\": {},\n    \"phases\": {{",
            self.calls, self.passes, self.before_hesse
        );
        for (i, (phase, count)) in self.phases.iter().enumerate() {
            print!(
                "{}\n      \"{phase}\": {count}",
                if i > 0 { "," } else { "" }
            );
        }
        print!("\n    }}");
        if self.trace {
            print!(",\n    \"trace\": [");
            for (i, entry) in self.evaluations.iter().enumerate() {
                print!("{}\n      {entry}", if i > 0 { "," } else { "" });
            }
            print!("\n    ]");
        }
    }
}
