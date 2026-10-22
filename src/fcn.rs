/// Function to be minimized.
///
/// The core trait that users implement. `error_def()` returns the change in FCN
/// that corresponds to one standard deviation (1.0 for chi-square, 0.5 for
/// log-likelihood).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GradientParameterSpace {
    Internal,
    External,
}

pub trait FCN {
    /// Evaluate the function at the given parameter values.
    fn value(&self, par: &[f64]) -> f64;

    /// FCN change corresponding to 1-sigma. Default = 1.0 (chi-square).
    fn error_def(&self) -> f64 {
        1.0
    }

    /// Alias for `error_def()`.
    fn up(&self) -> f64 {
        self.error_def()
    }

    /// Optional mutable setter for FCN error definition.
    ///
    /// Default no-op keeps backwards compatibility for immutable FCNs.
    fn set_error_def(&mut self, _up: f64) {}

    /// Whether this FCN provides analytical gradients.
    fn has_gradient(&self) -> bool {
        false
    }

    /// Optional gradient with access to previous derivative state.
    fn gradient_with_prev_result(
        &self,
        _par: &[f64],
        _previous_grad: Option<&[f64]>,
        _previous_g2: Option<&[f64]>,
        _previous_gstep: Option<&[f64]>,
    ) -> Vec<f64> {
        Vec::new()
    }

    /// Defines whether supplied derivatives live in internal or external space.
    fn grad_parameter_space(&self) -> GradientParameterSpace {
        GradientParameterSpace::External
    }

    /// Optional diagonal second-derivative contract.
    fn g2(&self, _par: &[f64]) -> Vec<f64> {
        Vec::new()
    }

    /// Optional Hessian contract (packed upper-triangle format).
    fn hessian(&self, _par: &[f64]) -> Vec<f64> {
        Vec::new()
    }

    /// Whether `hessian()` is implemented.
    fn has_hessian(&self) -> bool {
        false
    }

    /// Whether `g2()` is implemented.
    fn has_g2(&self) -> bool {
        false
    }
}

/// Blanket impl: any `Fn(&[f64]) -> f64` is a valid FCN with error_def = 1.0.
impl<F> FCN for F
where
    F: Fn(&[f64]) -> f64,
{
    fn value(&self, par: &[f64]) -> f64 {
        self(par)
    }
}

/// FCN that also provides analytical gradients.
pub trait FCNGradient: FCN {
    /// Compute the gradient vector at the given parameter values.
    fn gradient(&self, par: &[f64]) -> Vec<f64>;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn closure_as_fcn() {
        let f = |x: &[f64]| x[0] * x[0] + x[1] * x[1];
        assert!((FCN::value(&f, &[3.0, 4.0]) - 25.0).abs() < 1e-15);
        assert!((f.error_def() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn struct_fcn() {
        struct Rosenbrock;
        impl FCN for Rosenbrock {
            fn value(&self, p: &[f64]) -> f64 {
                (1.0 - p[0]).powi(2) + 100.0 * (p[1] - p[0] * p[0]).powi(2)
            }
            fn error_def(&self) -> f64 {
                0.5
            }
        }
        let f = Rosenbrock;
        assert!((f.value(&[1.0, 1.0])).abs() < 1e-15);
        assert!((f.error_def() - 0.5).abs() < 1e-15);
        assert!((f.up() - 0.5).abs() < 1e-15);
        assert!(!f.has_hessian());
        assert!(!f.has_g2());
        assert_eq!(f.grad_parameter_space(), GradientParameterSpace::External);
    }

    #[test]
    fn gradient_contract_defaults() {
        struct QuadGrad;
        impl FCN for QuadGrad {
            fn value(&self, p: &[f64]) -> f64 {
                p[0] * p[0]
            }
        }
        impl FCNGradient for QuadGrad {
            fn gradient(&self, p: &[f64]) -> Vec<f64> {
                vec![2.0 * p[0]]
            }
        }

        let f = QuadGrad;
        assert_eq!(FCNGradient::gradient(&f, &[3.0]), vec![6.0]);
        assert!(!FCN::has_gradient(&f));
        assert!(FCN::gradient_with_prev_result(&f, &[3.0], None, None, None).is_empty());
    }
}
