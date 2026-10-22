/// Strategy presets controlling gradient/Hessian calculation effort.
///
/// Three Minuit-style effort levels: low (0), medium (1), high (2).
/// Medium is the default.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MnStrategy {
    strategy: u32,
    grad_ncycles: u32,
    hess_ncycles: u32,
    hess_grad_ncycles: u32,
    grad_step_tol: u32, // stored as x10 to avoid float equality issues
    grad_tol: u32,      // stored as x100
    hess_step_tol: u32, // stored as x10
    hess_g2_tol: u32,   // stored as x100
    hess_cfd_g2: u32,
    hess_force_pos_def: u32,
    store_level: u32,
}

impl MnStrategy {
    /// Create a new strategy with the given level (0=low, 1=medium, 2=high).
    pub fn new(level: u32) -> Self {
        let mut s = Self {
            strategy: level,
            grad_ncycles: 0,
            hess_ncycles: 0,
            hess_grad_ncycles: 0,
            grad_step_tol: 0,
            grad_tol: 0,
            hess_step_tol: 0,
            hess_g2_tol: 0,
            hess_cfd_g2: 0,
            hess_force_pos_def: 1,
            store_level: 1,
        };
        match level {
            0 => s.set_low_strategy(),
            2 => s.set_high_strategy(),
            _ => s.set_medium_strategy(),
        }
        s
    }

    fn set_low_strategy(&mut self) {
        self.strategy = 0;
        self.grad_ncycles = 2;
        self.grad_step_tol = 5; // 0.5
        self.grad_tol = 10; // 0.1
        self.hess_ncycles = 3;
        self.hess_step_tol = 5; // 0.5
        self.hess_g2_tol = 10; // 0.1
        self.hess_grad_ncycles = 1;
        self.hess_cfd_g2 = 0;
        self.hess_force_pos_def = 1;
        self.store_level = 1;
    }

    fn set_medium_strategy(&mut self) {
        self.strategy = 1;
        self.grad_ncycles = 3;
        self.grad_step_tol = 3; // 0.3
        self.grad_tol = 5; // 0.05
        self.hess_ncycles = 5;
        self.hess_step_tol = 3; // 0.3
        self.hess_g2_tol = 5; // 0.05
        self.hess_grad_ncycles = 2;
        self.hess_cfd_g2 = 0;
        self.hess_force_pos_def = 1;
        self.store_level = 1;
    }

    fn set_high_strategy(&mut self) {
        self.strategy = 2;
        self.grad_ncycles = 5;
        self.grad_step_tol = 1; // 0.1
        self.grad_tol = 2; // 0.02
        self.hess_ncycles = 7;
        self.hess_step_tol = 1; // 0.1
        self.hess_g2_tol = 2; // 0.02
        self.hess_grad_ncycles = 6;
        self.hess_cfd_g2 = 0;
        self.hess_force_pos_def = 1;
        self.store_level = 1;
    }

    /// Get the strategy level.
    pub fn strategy(&self) -> u32 {
        self.strategy
    }

    /// Get the number of gradient calculation cycles.
    pub fn grad_ncycles(&self) -> u32 {
        self.grad_ncycles
    }

    pub fn gradient_ncycles(&self) -> u32 {
        self.grad_ncycles()
    }

    /// Get the gradient step tolerance.
    pub fn grad_step_tol(&self) -> f64 {
        self.grad_step_tol as f64 / 10.0
    }

    pub fn gradient_step_tolerance(&self) -> f64 {
        self.grad_step_tol()
    }

    /// Get the gradient tolerance.
    pub fn grad_tol(&self) -> f64 {
        self.grad_tol as f64 / 100.0
    }

    pub fn gradient_tolerance(&self) -> f64 {
        self.grad_tol()
    }

    /// Get the number of Hessian calculation cycles.
    pub fn hess_ncycles(&self) -> u32 {
        self.hess_ncycles
    }

    pub fn hessian_ncycles(&self) -> u32 {
        self.hess_ncycles()
    }

    /// Get the Hessian step tolerance.
    pub fn hess_step_tol(&self) -> f64 {
        self.hess_step_tol as f64 / 10.0
    }

    pub fn hessian_step_tolerance(&self) -> f64 {
        self.hess_step_tol()
    }

    /// Get the Hessian g2 tolerance.
    pub fn hess_g2_tol(&self) -> f64 {
        self.hess_g2_tol as f64 / 100.0
    }

    pub fn hessian_g2_tolerance(&self) -> f64 {
        self.hess_g2_tol()
    }

    /// Get the number of Hessian gradient calculation cycles.
    pub fn hess_grad_ncycles(&self) -> u32 {
        self.hess_grad_ncycles
    }

    pub fn hessian_gradient_ncycles(&self) -> u32 {
        self.hess_grad_ncycles()
    }

    pub fn storage_level(&self) -> u32 {
        self.store_level
    }

    pub fn hessian_central_fd_mixed_derivatives(&self) -> u32 {
        self.hess_cfd_g2
    }

    pub fn hessian_force_pos_def(&self) -> u32 {
        self.hess_force_pos_def
    }

    pub fn set_gradient_ncycles(&mut self, ncycles: u32) {
        self.grad_ncycles = ncycles;
    }

    pub fn set_gradient_step_tolerance(&mut self, tol: f64) {
        self.grad_step_tol = (tol.max(0.0) * 10.0).round() as u32;
    }

    pub fn set_gradient_tolerance(&mut self, tol: f64) {
        self.grad_tol = (tol.max(0.0) * 100.0).round() as u32;
    }

    pub fn set_hessian_ncycles(&mut self, ncycles: u32) {
        self.hess_ncycles = ncycles;
    }

    pub fn set_hessian_step_tolerance(&mut self, tol: f64) {
        self.hess_step_tol = (tol.max(0.0) * 10.0).round() as u32;
    }

    pub fn set_hessian_g2_tolerance(&mut self, tol: f64) {
        self.hess_g2_tol = (tol.max(0.0) * 100.0).round() as u32;
    }

    pub fn set_hessian_gradient_ncycles(&mut self, ncycles: u32) {
        self.hess_grad_ncycles = ncycles;
    }

    pub fn set_hessian_central_fd_mixed_derivatives(&mut self, flag: u32) {
        self.hess_cfd_g2 = flag;
    }

    pub fn set_hessian_force_pos_def(&mut self, flag: u32) {
        self.hess_force_pos_def = flag;
    }

    pub fn set_storage_level(&mut self, level: u32) {
        self.store_level = level;
    }

    /// Check if this is a low strategy.
    pub fn is_low(&self) -> bool {
        self.strategy == 0
    }

    /// Check if this is a medium strategy.
    pub fn is_medium(&self) -> bool {
        self.strategy == 1
    }

    /// Check if this is a high strategy.
    pub fn is_high(&self) -> bool {
        self.strategy >= 2
    }

    pub fn is_very_high(&self) -> bool {
        self.strategy >= 3
    }

    pub fn set_very_high_strategy(&mut self) {
        self.set_high_strategy();
        self.strategy = 3;
    }
}

impl Default for MnStrategy {
    fn default() -> Self {
        Self::new(1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn low_strategy() {
        let s = MnStrategy::new(0);
        assert!(s.is_low());
        assert_eq!(s.grad_ncycles(), 2);
        assert!((s.grad_step_tol() - 0.5).abs() < 1e-15);
        assert!((s.grad_tol() - 0.1).abs() < 1e-15);
        assert_eq!(s.hess_ncycles(), 3);
        assert_eq!(s.hess_grad_ncycles(), 1);
    }

    #[test]
    fn medium_strategy() {
        let s = MnStrategy::new(1);
        assert!(s.is_medium());
        assert_eq!(s.grad_ncycles(), 3);
        assert!((s.grad_step_tol() - 0.3).abs() < 1e-15);
        assert!((s.grad_tol() - 0.05).abs() < 1e-15);
        assert_eq!(s.hess_ncycles(), 5);
        assert_eq!(s.hess_grad_ncycles(), 2);
    }

    #[test]
    fn high_strategy() {
        let s = MnStrategy::new(2);
        assert!(s.is_high());
        assert_eq!(s.grad_ncycles(), 5);
        assert!((s.grad_step_tol() - 0.1).abs() < 1e-15);
        assert!((s.grad_tol() - 0.02).abs() < 1e-15);
        assert_eq!(s.hess_ncycles(), 7);
        assert_eq!(s.hess_grad_ncycles(), 6);
    }

    #[test]
    fn default_is_medium() {
        let s = MnStrategy::default();
        assert!(s.is_medium());
    }

    #[test]
    fn strategy_alias_getters_match_core_getters() {
        let s = MnStrategy::new(1);
        assert_eq!(s.gradient_ncycles(), s.grad_ncycles());
        assert!((s.gradient_step_tolerance() - s.grad_step_tol()).abs() < 1e-15);
        assert!((s.gradient_tolerance() - s.grad_tol()).abs() < 1e-15);
        assert_eq!(s.hessian_ncycles(), s.hess_ncycles());
        assert!((s.hessian_step_tolerance() - s.hess_step_tol()).abs() < 1e-15);
        assert!((s.hessian_g2_tolerance() - s.hess_g2_tol()).abs() < 1e-15);
        assert_eq!(s.hessian_gradient_ncycles(), s.hess_grad_ncycles());
        assert_eq!(s.storage_level(), 1);
    }

    #[test]
    fn strategy_setters_override_values() {
        let mut s = MnStrategy::new(1);
        s.set_gradient_ncycles(9);
        s.set_gradient_step_tolerance(0.7);
        s.set_gradient_tolerance(0.09);
        s.set_hessian_ncycles(11);
        s.set_hessian_step_tolerance(0.8);
        s.set_hessian_g2_tolerance(0.11);
        s.set_hessian_gradient_ncycles(12);
        s.set_hessian_central_fd_mixed_derivatives(1);
        s.set_hessian_force_pos_def(0);
        s.set_storage_level(3);

        assert_eq!(s.gradient_ncycles(), 9);
        assert!((s.gradient_step_tolerance() - 0.7).abs() < 1e-15);
        assert!((s.gradient_tolerance() - 0.09).abs() < 1e-15);
        assert_eq!(s.hessian_ncycles(), 11);
        assert!((s.hessian_step_tolerance() - 0.8).abs() < 1e-15);
        assert!((s.hessian_g2_tolerance() - 0.11).abs() < 1e-15);
        assert_eq!(s.hessian_gradient_ncycles(), 12);
        assert_eq!(s.hessian_central_fd_mixed_derivatives(), 1);
        assert_eq!(s.hessian_force_pos_def(), 0);
        assert_eq!(s.storage_level(), 3);
    }

    #[test]
    fn very_high_strategy() {
        let mut s = MnStrategy::new(1);
        assert!(!s.is_very_high());
        s.set_very_high_strategy();
        assert!(s.is_very_high());
        assert!(s.is_high());
    }
}
