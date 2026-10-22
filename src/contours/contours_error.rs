//! ContoursError: result type for 2D contour computation.
//!
//! Contains the contour points and the MINOS errors for both axes.

use crate::minos::MinosError;

/// Result of a 2D contour computation.
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
    /// Total function calls.
    pub nfcn: usize,
}

impl ContoursError {
    pub fn xpar(&self) -> usize {
        self.par_x
    }

    pub fn ypar(&self) -> usize {
        self.par_y
    }

    pub fn nfcn(&self) -> usize {
        self.nfcn
    }

    pub fn x_min(&self) -> f64 {
        self.x_minos.min()
    }

    pub fn y_min(&self) -> f64 {
        self.y_minos.min()
    }
}
