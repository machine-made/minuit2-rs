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

/// JSON string literal for any text: quote, backslash and control characters escaped.
pub fn json_string(text: &str) -> String {
    let mut out = String::from("\"");
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04X}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// Provenance block of a stored file: the program, how the library was built, and
/// the commands that build and run the program.
pub fn provenance(program: &str, build_command: &str, run_command: &str) -> String {
    format!(
        "{{\n    \"program\": {},\n    \"root_tag\": null,\n    \"library_build_command\": \"cargo build --locked --lib\",\n    \"build_command\": {},\n    \"command\": {}\n  }}",
        json_string(program),
        json_string(build_command),
        json_string(run_command)
    )
}

/// A record is the whole output of one run, in the layout of the stored files:
/// two-space indent, one key per line, a trailing newline. The caller prints the
/// result fields, then `PhaseCounts::print`, then `end_record`.
pub fn begin_record(provenance: &str, workload: &str) {
    print!(
        "{{\n  \"provenance\": {provenance},\n  \"result\": {{\n    \"workload\": {}",
        json_string(workload)
    );
}

/// Values of a field: integers, bool and f64 only, and every f64 goes through `json_number`.
pub trait JsonValue {
    fn json(&self) -> String;
}

impl JsonValue for usize {
    fn json(&self) -> String {
        self.to_string()
    }
}

impl JsonValue for bool {
    fn json(&self) -> String {
        self.to_string()
    }
}

impl JsonValue for f64 {
    fn json(&self) -> String {
        json_number(*self)
    }
}

pub fn field(name: &str, value: impl JsonValue) {
    print!(",\n    {}: {}", json_string(name), value.json());
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
                "{{\"call\":{},\"phase\":{},\"p\":[{}]}}",
                self.calls,
                json_string(&phase),
                coordinates.join(",")
            ));
        }
    }

    /// The invariant of the call-count change: Minuit2 reports exactly the calls the user function saw.
    pub fn check_reported(&self, reported_nfcn: usize) -> Result<(), String> {
        if reported_nfcn == self.calls {
            Ok(())
        } else {
            Err(format!(
                "reported_nfcn {reported_nfcn} differs from user_calls {}",
                self.calls
            ))
        }
    }

    pub fn print(&self) {
        print!(
            ",\n    \"user_calls\": {},\n    \"hesse_passes\": {},\n    \"before_hesse\": {},\n    \"phases\": {{",
            self.calls, self.passes, self.before_hesse
        );
        for (i, (phase, count)) in self.phases.iter().enumerate() {
            print!(
                "{}\n      {}: {count}",
                if i > 0 { "," } else { "" },
                json_string(phase)
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
