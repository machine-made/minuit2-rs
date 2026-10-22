//! Simplex vertex storage with best/worst tracking.
//!
//! Stores `N + 1` vertices, each as a parameter vector plus function value,
//! and tracks the best and worst function values.

#[derive(Debug, Clone)]
pub struct SimplexParameters {
    /// (function_value, parameter_vector) pairs for each vertex.
    params: Vec<(f64, Vec<f64>)>,
    /// Index of vertex with lowest function value.
    jlow: usize,
    /// Index of vertex with highest function value.
    jhigh: usize,
}

impl SimplexParameters {
    pub fn new(params: Vec<(f64, Vec<f64>)>) -> Self {
        let (jlow, jhigh) = Self::find_extremes(&params);
        Self {
            params,
            jlow,
            jhigh,
        }
    }

    fn find_extremes(params: &[(f64, Vec<f64>)]) -> (usize, usize) {
        let mut jlow = 0;
        let mut jhigh = 0;
        for (i, (fval, _)) in params.iter().enumerate() {
            if *fval < params[jlow].0 {
                jlow = i;
            }
            if *fval > params[jhigh].0 {
                jhigh = i;
            }
        }
        (jlow, jhigh)
    }

    /// Update a vertex and recompute extremes.
    pub fn update(&mut self, index: usize, fval: f64, vec: Vec<f64>) {
        self.params[index] = (fval, vec);
        let (jlow, jhigh) = Self::find_extremes(&self.params);
        self.jlow = jlow;
        self.jhigh = jhigh;
    }

    pub fn params(&self) -> &[(f64, Vec<f64>)] {
        &self.params
    }

    /// Index of best (lowest) vertex.
    pub fn jlow(&self) -> usize {
        self.jlow
    }

    pub fn jl(&self) -> usize {
        self.jlow()
    }

    /// Index of worst (highest) vertex.
    pub fn jhigh(&self) -> usize {
        self.jhigh
    }

    pub fn jh(&self) -> usize {
        self.jhigh()
    }

    /// Function value at best vertex.
    pub fn fval_best(&self) -> f64 {
        self.params[self.jlow].0
    }

    /// Function value at worst vertex.
    pub fn fval_worst(&self) -> f64 {
        self.params[self.jhigh].0
    }

    /// Parameter vector at best vertex.
    pub fn best(&self) -> &[f64] {
        &self.params[self.jlow].1
    }

    /// Compatibility helper for the current direction-like span in simplex space.
    pub fn dirin(&self) -> Vec<f64> {
        if self.params.is_empty() {
            return Vec::new();
        }
        let best = &self.params[self.jlow].1;
        let worst = &self.params[self.jhigh].1;
        best.iter().zip(worst).map(|(b, w)| w - b).collect()
    }

    /// Number of vertices (N+1).
    pub fn len(&self) -> usize {
        self.params.len()
    }

    /// Whether the simplex has no vertices.
    pub fn is_empty(&self) -> bool {
        self.params.is_empty()
    }

    /// EDM estimate: difference between worst and best function values.
    pub fn edm(&self) -> f64 {
        self.fval_worst() - self.fval_best()
    }

    /// Index of second-worst vertex.
    pub fn jsecond_high(&self) -> usize {
        let mut jsec = if self.jhigh == 0 { 1 } else { 0 };
        for (i, (fval, _)) in self.params.iter().enumerate() {
            if i != self.jhigh && *fval > self.params[jsec].0 {
                jsec = i;
            }
        }
        jsec
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extremes() {
        let params = vec![
            (3.0, vec![1.0, 0.0]),
            (1.0, vec![0.0, 1.0]),
            (5.0, vec![1.0, 1.0]),
        ];
        let sp = SimplexParameters::new(params);
        assert_eq!(sp.jlow(), 1);
        assert_eq!(sp.jhigh(), 2);
        assert!((sp.fval_best() - 1.0).abs() < 1e-15);
        assert!((sp.fval_worst() - 5.0).abs() < 1e-15);
    }

    #[test]
    fn update_recomputes() {
        let params = vec![(3.0, vec![1.0]), (1.0, vec![2.0]), (5.0, vec![3.0])];
        let mut sp = SimplexParameters::new(params);
        sp.update(2, 0.5, vec![4.0]); // vertex 2 is now best
        assert_eq!(sp.jlow(), 2);
        assert!((sp.fval_best() - 0.5).abs() < 1e-15);
    }
}
