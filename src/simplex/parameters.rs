//! Simplex vertex storage with best/worst tracking.
//!
//! Stores `N + 1` vertices, each as a parameter vector plus function value,
//! and tracks the best and worst function values.

/// Simplex vertices with their objective values and best and worst indices.
/// # Examples
///
/// ```
/// use minuit2::simplex::parameters::SimplexParameters;
/// let s = SimplexParameters::new(vec![(4.0, vec![2.0]), (1.0, vec![1.0])]);
/// assert_eq!(s.jlow(), 1);
/// assert!((s.edm() - 3.0).abs() < 1e-12);
/// ```
#[derive(Debug, Clone)]
pub struct SimplexParameters {
    /// (`function_value`, `parameter_vector`) pairs for each vertex.
    params: Vec<(f64, Vec<f64>)>,
    /// Index of vertex with lowest function value.
    jlow: usize,
    /// Index of vertex with highest function value.
    jhigh: usize,
}

impl SimplexParameters {
    /// Store simplex vertices as objective values paired with internal parameter vectors.
    #[must_use]
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
    ///
    /// # Panics
    /// Panics if the vertex index is out of range.
    pub fn update(&mut self, index: usize, fval: f64, vec: Vec<f64>) {
        self.params[index] = (fval, vec);
        let (jlow, jhigh) = Self::find_extremes(&self.params);
        self.jlow = jlow;
        self.jhigh = jhigh;
    }

    /// Return the stored objective values and internal parameter vectors.
    #[must_use]
    pub fn params(&self) -> &[(f64, Vec<f64>)] {
        &self.params
    }

    /// Index of best (lowest) vertex.
    #[must_use]
    pub fn jlow(&self) -> usize {
        self.jlow
    }

    /// Return the index of the best vertex.
    #[must_use]
    pub fn jl(&self) -> usize {
        self.jlow()
    }

    /// Index of worst (highest) vertex.
    #[must_use]
    pub fn jhigh(&self) -> usize {
        self.jhigh
    }

    /// Return the index of the worst vertex.
    #[must_use]
    pub fn jh(&self) -> usize {
        self.jhigh()
    }

    /// Function value at best vertex.
    ///
    /// # Panics
    /// Panics if the simplex contains no vertices.
    #[must_use]
    pub fn fval_best(&self) -> f64 {
        self.params[self.jlow].0
    }

    /// Function value at worst vertex.
    ///
    /// # Panics
    /// Panics if the simplex contains no vertices.
    #[must_use]
    pub fn fval_worst(&self) -> f64 {
        self.params[self.jhigh].0
    }

    /// Parameter vector at best vertex.
    ///
    /// # Panics
    /// Panics if the simplex contains no vertices.
    #[must_use]
    pub fn best(&self) -> &[f64] {
        &self.params[self.jlow].1
    }

    /// Compatibility helper for the current direction-like span in simplex space.
    #[must_use]
    pub fn dirin(&self) -> Vec<f64> {
        if self.params.is_empty() {
            return Vec::new();
        }
        let best = &self.params[self.jlow].1;
        let worst = &self.params[self.jhigh].1;
        best.iter().zip(worst).map(|(b, w)| w - b).collect()
    }

    /// Number of vertices (N+1).
    #[must_use]
    pub fn len(&self) -> usize {
        self.params.len()
    }

    /// Whether the simplex has no vertices.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.params.is_empty()
    }

    /// EDM estimate: difference between worst and best function values.
    ///
    /// # Panics
    /// Panics if the simplex contains no vertices.
    #[must_use]
    pub fn edm(&self) -> f64 {
        self.fval_worst() - self.fval_best()
    }

    /// Index of second-worst vertex.
    ///
    /// With fewer than two vertices, returns one, which is outside the vertex range.
    #[must_use]
    pub fn jsecond_high(&self) -> usize {
        let mut jsec = usize::from(self.jhigh == 0);
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
