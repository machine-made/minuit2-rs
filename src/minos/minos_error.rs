//! `MinosError`: paired upper/lower crossing results for one parameter.
//!
//! Wraps two crossing results (upper and lower) and computes the actual
//! parameter error in external space.

use super::cross::MnCross;

/// Asymmetric MINOS errors for a single parameter.
/// # Examples
///
/// ```
/// use minuit2::MnMigrad;
/// let f = |p: &[f64]| (p[0] - 1.0).powi(2);
/// let minimum = MnMigrad::new().add("x", 0.0, 0.1).minimize(&f);
/// let errors = minuit2::MnMinos::new(&f, &minimum).minos(0);
/// assert!(errors.is_valid());
/// assert!(errors.lower_error() < 0.0);
/// assert!(errors.upper_error() > 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct MinosError {
    /// External parameter index.
    parameter: usize,
    /// Parameter value at the original minimum.
    min: f64,
    /// Hesse (parabolic) error for this parameter.
    hesse_error: f64,
    /// Lower crossing result.
    lower: MnCross,
    /// Upper crossing result.
    upper: MnCross,
}

impl MinosError {
    /// Store both crossings, the external parameter index, the minimum value, and its parabolic error.
    #[must_use]
    pub fn new(
        parameter: usize,
        min: f64,
        hesse_error: f64,
        lower: MnCross,
        upper: MnCross,
    ) -> Self {
        Self {
            parameter,
            min,
            hesse_error,
            lower,
            upper,
        }
    }

    /// The lower (negative) MINOS error.
    ///
    /// Returns: -err * (1 + lower.value) if valid, else -`hesse_error`.
    #[must_use]
    pub fn lower_error(&self) -> f64 {
        if self.lower.is_valid() {
            -self.hesse_error * (1.0 + self.lower.value())
        } else {
            -self.hesse_error
        }
    }

    /// The upper (positive) MINOS error.
    ///
    /// Returns: err * (1 + upper.value) if valid, else `hesse_error`.
    #[must_use]
    pub fn upper_error(&self) -> f64 {
        if self.upper.is_valid() {
            self.hesse_error * (1.0 + self.upper.value())
        } else {
            self.hesse_error
        }
    }

    /// Return the external parameter index.
    #[must_use]
    pub fn parameter(&self) -> usize {
        self.parameter
    }

    /// Return the lower crossing result.
    #[must_use]
    pub fn lower(&self) -> &MnCross {
        &self.lower
    }

    /// Return the upper crossing result.
    #[must_use]
    pub fn upper(&self) -> &MnCross {
        &self.upper
    }

    /// Check if both upper and lower errors are valid.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.lower.is_valid() && self.upper.is_valid()
    }

    /// Return whether the lower crossing is valid.
    #[must_use]
    pub fn lower_valid(&self) -> bool {
        self.lower.is_valid()
    }

    /// Return whether the upper crossing is valid.
    #[must_use]
    pub fn upper_valid(&self) -> bool {
        self.upper.is_valid()
    }

    /// Check if the lower crossing is at a limit.
    #[must_use]
    pub fn at_lower_limit(&self) -> bool {
        self.lower.at_limit()
    }

    /// Check if the upper crossing is at a limit.
    #[must_use]
    pub fn at_upper_limit(&self) -> bool {
        self.upper.at_limit()
    }

    /// Return whether the lower crossing reached its call limit.
    #[must_use]
    pub fn at_lower_max_fcn(&self) -> bool {
        self.lower.at_max_fcn()
    }

    /// Return whether the upper crossing reached its call limit.
    #[must_use]
    pub fn at_upper_max_fcn(&self) -> bool {
        self.upper.at_max_fcn()
    }

    /// Return whether the lower search found a new minimum.
    #[must_use]
    pub fn lower_new_min(&self) -> bool {
        self.lower.new_minimum()
    }

    /// Return whether the upper search found a new minimum.
    #[must_use]
    pub fn upper_new_min(&self) -> bool {
        self.upper.new_minimum()
    }

    /// Return the number of objective evaluations.
    #[must_use]
    pub fn nfcn(&self) -> usize {
        self.lower.nfcn() + self.upper.nfcn()
    }

    /// Parameter value at the minimum.
    #[must_use]
    pub fn min(&self) -> f64 {
        self.min
    }
}
