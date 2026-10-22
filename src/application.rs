//! Shared minimization logic: default maxfcn calculation.
//!
//! The default maximum number of function calls is `200 + 100*n + 5*n^2`,
//! where `n` is the number of variable parameters.

/// Compute default maximum function calls for `n` variable parameters.
pub fn default_max_fcn(n: usize) -> usize {
    200 + 100 * n + 5 * n * n
}

/// Default tolerance.
pub const DEFAULT_TOLERANCE: f64 = 0.1;
