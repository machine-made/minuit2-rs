// Evidence-only stack classification. Compile the library and harness with
// debug information and without optimisation so Minuit2 frames remain visible.
// Centre, diagonal and mixed calls use the same coordinate rule as phase_counts.h.
use std::collections::BTreeMap;

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
                    _ => panic!("unclassified Hesse call"),
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
            self.evaluations.push(format!(
                "{{\"call\":{},\"phase\":\"{phase}\",\"p\":{p:?}}}",
                self.calls
            ));
        }
    }

    pub fn print(&self) {
        print!(
            ",\"user_calls\":{},\"hesse_passes\":{},\"before_hesse\":{},\"phases\":{{",
            self.calls, self.passes, self.before_hesse
        );
        for (i, (phase, count)) in self.phases.iter().enumerate() {
            if i > 0 {
                print!(",");
            }
            print!("\"{phase}\":{count}");
        }
        print!("}}");
        if self.trace {
            print!(",\"trace\":[{}]", self.evaluations.join(","));
        }
    }
}
