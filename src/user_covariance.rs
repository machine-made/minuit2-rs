/// User-level covariance matrix stored as upper triangle.
///
/// Stores a symmetric matrix as `n * (n + 1) / 2` entries, by columns of its upper triangle.
/// Each entry has the product of the corresponding external parameter units.
/// # Examples
///
/// ```
/// use minuit2::MnUserCovariance;
/// let mut cov = MnUserCovariance::new(2);
/// cov.set(0, 1, 0.25);
/// assert!((cov.get(1, 0) - 0.25).abs() < 1e-12);
/// assert_eq!(cov.size(), 3);
/// ```
#[derive(Debug, Clone)]
pub struct MnUserCovariance {
    data: Vec<f64>,
    nrow: usize,
}

impl MnUserCovariance {
    /// Create a zero covariance matrix for `n` parameters.
    #[must_use]
    pub fn new(n: usize) -> Self {
        Self {
            data: vec![0.0; n * (n + 1) / 2],
            nrow: n,
        }
    }

    /// Create from raw upper-triangle data.
    ///
    /// # Panics
    /// Panics if the data length differs from `n * (n + 1) / 2`.
    #[must_use]
    pub fn from_vec(data: Vec<f64>, n: usize) -> Self {
        assert_eq!(data.len(), n * (n + 1) / 2, "data size mismatch");
        Self { data, nrow: n }
    }

    /// Get the number of rows (parameters).
    #[must_use]
    pub fn nrow(&self) -> usize {
        self.nrow
    }

    /// Access element (row, col). Symmetric: (i,j) == (j,i).
    ///
    /// # Panics
    /// Panics if the computed packed index is outside the stored covariance data.
    #[must_use]
    pub fn get(&self, row: usize, col: usize) -> f64 {
        let (r, c) = if row <= col { (row, col) } else { (col, row) };
        self.data[r + c * (c + 1) / 2]
    }

    /// Set element (row, col). Also sets (col, row).
    ///
    /// # Panics
    /// Panics if the computed packed index is outside the stored covariance data.
    pub fn set(&mut self, row: usize, col: usize, val: f64) {
        let (r, c) = if row <= col { (row, col) } else { (col, row) };
        self.data[r + c * (c + 1) / 2] = val;
    }

    /// Get the raw upper-triangle data.
    #[must_use]
    pub fn data(&self) -> &[f64] {
        &self.data
    }

    /// Return mutable packed covariance data in column order of the upper triangle.
    pub fn data_mut(&mut self) -> &mut [f64] {
        &mut self.data
    }

    /// Multiply every covariance entry by the supplied dimensionless factor.
    pub fn scale(&mut self, value: f64) {
        for x in &mut self.data {
            *x *= value;
        }
    }

    /// Return the number of stored upper triangle entries.
    #[must_use]
    pub fn size(&self) -> usize {
        self.data.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symmetric_access() {
        let mut cov = MnUserCovariance::new(3);
        cov.set(0, 1, 0.5);
        assert!((cov.get(0, 1) - 0.5).abs() < 1e-15);
        assert!((cov.get(1, 0) - 0.5).abs() < 1e-15);
    }

    #[test]
    fn diagonal() {
        let mut cov = MnUserCovariance::new(2);
        cov.set(0, 0, 1.0);
        cov.set(1, 1, 4.0);
        assert!((cov.get(0, 0) - 1.0).abs() < 1e-15);
        assert!((cov.get(1, 1) - 4.0).abs() < 1e-15);
    }

    #[test]
    fn data_length() {
        let cov = MnUserCovariance::new(4);
        assert_eq!(cov.data().len(), 10); // 4*5/2
    }
}
