//! `ContoursError`: result type for 2D contour computation.
//!
//! Contains the contour points and the MINOS errors for both axes.

use crate::minos::MinosError;

/// Result of a 2D contour computation.
/// # Examples
///
/// ```
/// use minuit2::{MnMigrad, MnContours};
/// let f = |p: &[f64]| p[0] * p[0] + p[1] * p[1];
/// let minimum = MnMigrad::new().add("x", 0.1, 0.1).add("y", 0.1, 0.1).minimize(&f);
/// let result = MnContours::new(&f, &minimum).contour(0, 1, 4);
/// assert_eq!(result.points.len(), 4);
/// assert_eq!(result.xpar(), 0);
/// ```
#[derive(Debug, Clone)]
pub struct ContoursError {
    /// External index of x parameter.
    pub par_x: usize,
    /// External index of y parameter.
    pub par_y: usize,
    /// Contour points as (x, y) pairs in external space.
    pub points: Vec<(f64, f64)>,
    /// MINOS errors for the x parameter.
    pub x_minos: MinosError,
    /// MINOS errors for the y parameter.
    pub y_minos: MinosError,
    /// Stored call count; the current contour builder sets this to zero.
    pub nfcn: usize,
}

impl ContoursError {
    /// Return the external index of the horizontal parameter.
    #[must_use]
    pub fn xpar(&self) -> usize {
        self.par_x
    }

    /// Return the external index of the vertical parameter.
    #[must_use]
    pub fn ypar(&self) -> usize {
        self.par_y
    }

    /// Return the number of objective evaluations.
    #[must_use]
    pub fn nfcn(&self) -> usize {
        self.nfcn
    }

    /// Return the horizontal parameter value at the original minimum.
    #[must_use]
    pub fn x_min(&self) -> f64 {
        self.x_minos.min()
    }

    /// Return the vertical parameter value at the original minimum.
    #[must_use]
    pub fn y_min(&self) -> f64 {
        self.y_minos.min()
    }
}
