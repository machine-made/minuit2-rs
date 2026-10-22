//! MnContours: 2D confidence contours.
//!
//! Computes points tracing a contour where `F(x, y) = Fmin + Up` in the
//! `(par_x, par_y)` plane.

pub mod contours_error;

pub use contours_error::ContoursError;

use crate::fcn::FCN;
use crate::minimum::FunctionMinimum;
use crate::minos::MnMinos;
use crate::minos::minos_error::MinosError;
use crate::strategy::MnStrategy;

/// Compute 2D confidence contours.
pub struct MnContours<'a> {
    fcn: &'a dyn FCN,
    minimum: &'a FunctionMinimum,
    strategy: MnStrategy,
}

impl<'a> MnContours<'a> {
    /// Create a new contour calculator.
    pub fn new(fcn: &'a dyn FCN, minimum: &'a FunctionMinimum) -> Self {
        Self {
            fcn,
            minimum,
            strategy: MnStrategy::default(),
        }
    }

    /// Set strategy level.
    pub fn with_strategy(mut self, level: u32) -> Self {
        self.strategy = MnStrategy::new(level);
        self
    }

    /// Compute contour points for parameters `par_x` and `par_y`.
    ///
    /// Returns `npoints` points tracing the F = Fmin + Up contour.
    /// Minimum 4 points (the MINOS cardinal points).
    pub fn points(&self, par_x: usize, par_y: usize, npoints: usize) -> Vec<(f64, f64)> {
        let npoints = npoints.max(4);
        let nvar = self.minimum.n_variable_params();
        let _maxcalls = 100 * (npoints + 5) * (nvar + 1);

        let up = self.minimum.up();
        let user_state = self.minimum.user_state();

        let (x_minos, y_minos) = self.minos_errors(par_x, par_y);

        if !x_minos.is_valid() || !y_minos.is_valid() {
            return Vec::new();
        }

        // 4 cardinal points from MINOS
        let x_val = user_state.parameter(par_x).value();
        let y_val = user_state.parameter(par_y).value();

        // x_upper: fix x at upper, find y
        // x_lower: fix x at lower, find y
        // y_upper: fix y at upper, find x
        // y_lower: fix y at lower, find x
        let x_up = x_val + x_minos.upper_error();
        let x_lo = x_val + x_minos.lower_error(); // lower_error is negative
        let y_up = y_val + y_minos.upper_error();
        let y_lo = y_val + y_minos.lower_error();

        // The 4 cardinal points (approximate — on the contour boundary)
        let mut pts = vec![
            (x_up, y_val), // right
            (x_val, y_up), // top
            (x_lo, y_val), // left
            (x_val, y_lo), // bottom
        ];

        if npoints <= 4 {
            return pts;
        }

        // Scale factors for distance computation
        let scalx = if (x_up - x_lo).abs() > 1e-15 {
            1.0 / (x_up - x_lo)
        } else {
            1.0
        };
        let scaly = if (y_up - y_lo).abs() > 1e-15 {
            1.0 / (y_up - y_lo)
        } else {
            1.0
        };

        // Add more points by bisecting largest gaps
        let remaining = npoints - 4;
        for _ in 0..remaining {
            if pts.len() < 2 {
                break;
            }

            // Find largest gap (in scaled distance)
            let mut max_dist = 0.0_f64;
            let mut max_idx = 0;
            for i in 0..pts.len() {
                let j = (i + 1) % pts.len();
                let dx = (pts[j].0 - pts[i].0) * scalx;
                let dy = (pts[j].1 - pts[i].1) * scaly;
                let dist = (dx * dx + dy * dy).sqrt();
                if dist > max_dist {
                    max_dist = dist;
                    max_idx = i;
                }
            }

            let j = (max_idx + 1) % pts.len();
            let mid_x = 0.5 * (pts[max_idx].0 + pts[j].0);
            let mid_y = 0.5 * (pts[max_idx].1 + pts[j].1);

            // Perpendicular direction from center
            let dx = pts[j].0 - pts[max_idx].0;
            let dy = pts[j].1 - pts[max_idx].1;
            // Normal to the gap segment, pointing outward from minimum
            let nx = -dy * scalx;
            let ny = dx * scaly;
            let norm = (nx * nx + ny * ny).sqrt();
            if norm < 1e-15 {
                continue;
            }

            // Try to find a contour point near the midpoint.
            // Use MnFunctionCross along the perpendicular from the minimum.
            // For simplicity, we do a scan along the perpendicular direction
            // from the center (x_val, y_val) through the midpoint.
            let dir_x = mid_x - x_val;
            let dir_y = mid_y - y_val;

            // Evaluate function at midpoint to see if we're close to contour
            let nparams = user_state.len();
            let mut pars: Vec<f64> = (0..nparams)
                .map(|i| user_state.parameter(i).value())
                .collect();
            pars[par_x] = mid_x;
            pars[par_y] = mid_y;
            let f_mid = self.fcn.value(&pars);

            // Adjust: scale to hit the contour F = fmin + up
            let fmin = self.minimum.fval();
            let target = fmin + up;
            let ratio = if (f_mid - fmin).abs() > 1e-15 {
                (target / (f_mid - fmin)).sqrt()
            } else {
                1.0
            };

            let new_x = x_val + dir_x * ratio;
            let new_y = y_val + dir_y * ratio;

            // Check it's a reasonable point
            let seg_dist =
                ((new_x - pts[max_idx].0).powi(2) + (new_y - pts[max_idx].1).powi(2)).sqrt();
            if seg_dist < 1e-10 {
                continue;
            }

            // Insert after max_idx
            pts.insert(max_idx + 1, (new_x, new_y));
        }

        pts
    }

    /// Compute full contour with MINOS errors for both axes.
    pub fn contour(&self, par_x: usize, par_y: usize, npoints: usize) -> ContoursError {
        let (x_minos, y_minos) = self.minos_errors(par_x, par_y);

        let pts = self.points(par_x, par_y, npoints);

        ContoursError {
            par_x,
            par_y,
            points: pts,
            x_minos,
            y_minos,
            nfcn: 0,
        }
    }

    fn minos_errors(&self, par_x: usize, par_y: usize) -> (MinosError, MinosError) {
        let minos = MnMinos::new(self.fcn, self.minimum).with_strategy(self.strategy.strategy());

        (minos.minos_error(par_x), minos.minos_error(par_y))
    }
}
