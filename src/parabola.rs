//! Quadratic parabola utilities for line search.
//!
//! A simple quadratic `f(x) = a*x^2 + b*x + c` used during parabolic line
//! search to interpolate function values along a search direction.

/// A point (x, y) on the parabola.
#[derive(Debug, Clone, Copy)]
pub struct MnParabolaPoint {
    pub x: f64,
    pub y: f64,
}

impl MnParabolaPoint {
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    pub fn x(&self) -> f64 {
        self.x
    }

    pub fn y(&self) -> f64 {
        self.y
    }
}

/// Quadratic f(x) = a*x² + b*x + c.
#[derive(Debug, Clone, Copy)]
pub struct MnParabola {
    a: f64,
    b: f64,
    c: f64,
}

impl MnParabola {
    pub fn new(a: f64, b: f64, c: f64) -> Self {
        Self { a, b, c }
    }

    /// x coordinate of the minimum: -b / (2a).
    pub fn min(&self) -> f64 {
        -self.b / (2.0 * self.a)
    }

    /// y value at the minimum.
    pub fn y_min(&self) -> f64 {
        let xmin = self.min();
        self.a * xmin * xmin + self.b * xmin + self.c
    }

    /// Evaluate at a given x.
    pub fn y(&self, x: f64) -> f64 {
        self.a * x * x + self.b * x + self.c
    }

    /// Larger x root for a target y value.
    pub fn x_pos(&self, y: f64) -> f64 {
        (y / self.a + self.min() * self.min() - self.c / self.a).sqrt() + self.min()
    }

    /// Smaller x root for a target y value.
    pub fn x_neg(&self, y: f64) -> f64 {
        -(y / self.a + self.min() * self.min() - self.c / self.a).sqrt() + self.min()
    }

    /// Get the quadratic coefficient.
    pub fn a(&self) -> f64 {
        self.a
    }

    /// Get the linear coefficient.
    pub fn b(&self) -> f64 {
        self.b
    }

    /// Get the constant coefficient.
    pub fn c(&self) -> f64 {
        self.c
    }
}

/// Fit a parabola through two points with a known derivative at the first.
///
/// Given (x1, y1) with dfdx at x1, and (x2, y2), solve for a, b, c.
pub fn from_2_points_gradient(
    p1: MnParabolaPoint,
    p2: MnParabolaPoint,
    dfdx_at_p1: f64,
) -> MnParabola {
    // f(x) = a*x² + b*x + c
    // f'(x1) = 2*a*x1 + b = dfdx_at_p1
    // f(x1) = a*x1² + b*x1 + c = y1
    // f(x2) = a*x2² + b*x2 + c = y2
    //
    // From the first two:  b = dfdx_at_p1 - 2*a*x1
    // Substituting into the third:
    //   a*x2² + (dfdx_at_p1 - 2*a*x1)*x2 + c = y2
    //   a*(x2² - 2*x1*x2) + dfdx_at_p1*x2 + c = y2
    // From the second: c = y1 - a*x1² - b*x1 = y1 - a*x1² - (dfdx_at_p1 - 2*a*x1)*x1
    //   c = y1 - a*x1² - dfdx_at_p1*x1 + 2*a*x1² = y1 + a*x1² - dfdx_at_p1*x1
    // Substituting:
    //   a*(x2² - 2*x1*x2) + dfdx_at_p1*x2 + y1 + a*x1² - dfdx_at_p1*x1 = y2
    //   a*(x2² - 2*x1*x2 + x1²) = y2 - y1 - dfdx_at_p1*(x2 - x1)
    //   a*(x2 - x1)² = y2 - y1 - dfdx_at_p1*(x2 - x1)
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let a = (dy - dfdx_at_p1 * dx) / (dx * dx);
    let b = dfdx_at_p1 - 2.0 * a * p1.x;
    let c = p1.y - a * p1.x * p1.x - b * p1.x;
    MnParabola::new(a, b, c)
}

/// Fit a parabola through three points using the C++ MnParabola::operator() formula.
pub fn from_3_points(p1: MnParabolaPoint, p2: MnParabolaPoint, p3: MnParabolaPoint) -> MnParabola {
    // Standard Lagrange interpolation for 3 points, rearranged to ax²+bx+c
    let x1 = p1.x;
    let x2 = p2.x;
    let x3 = p3.x;
    let y1 = p1.y;
    let y2 = p2.y;
    let y3 = p3.y;

    let d12 = x1 - x2;
    let d13 = x1 - x3;
    let d23 = x2 - x3;

    // Lagrange interpolation coefficients:
    // L1 = y1/((x1-x2)(x1-x3)), L2 = y2/((x2-x1)(x2-x3)), L3 = y3/((x3-x1)(x3-x2))
    let l1 = y1 / (d12 * d13);
    let l2 = y2 / ((-d12) * d23);
    let l3 = y3 / ((-d13) * (-d23));
    let a = l1 + l2 + l3;
    let b = -(x2 + x3) * l1 - (x1 + x3) * l2 - (x1 + x2) * l3;
    let c = x2 * x3 * l1 + x1 * x3 * l2 + x1 * x2 * l3;

    MnParabola::new(a, b, c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parabola_minimum() {
        // f(x) = (x-3)^2 + 1 = x² - 6x + 10
        let p = MnParabola::new(1.0, -6.0, 10.0);
        assert!((p.min() - 3.0).abs() < 1e-14);
        assert!((p.y_min() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn from_2_points_with_gradient() {
        // f(x) = x² → f(0)=0, f'(0)=0, f(1)=1
        let p = from_2_points_gradient(
            MnParabolaPoint::new(0.0, 0.0),
            MnParabolaPoint::new(1.0, 1.0),
            0.0,
        );
        assert!((p.a() - 1.0).abs() < 1e-14);
        assert!(p.b().abs() < 1e-14);
        assert!(p.c().abs() < 1e-14);
        assert!(p.min().abs() < 1e-14);
    }

    #[test]
    fn from_3_points_recovery() {
        // f(x) = 2(x-1)^2 + 3 = 2x² - 4x + 5
        // f(0)=5, f(1)=3, f(2)=5
        let p = from_3_points(
            MnParabolaPoint::new(0.0, 5.0),
            MnParabolaPoint::new(1.0, 3.0),
            MnParabolaPoint::new(2.0, 5.0),
        );
        assert!((p.a() - 2.0).abs() < 1e-12);
        assert!((p.b() - (-4.0)).abs() < 1e-12);
        assert!((p.c() - 5.0).abs() < 1e-12);
        assert!((p.min() - 1.0).abs() < 1e-12);
        assert!((p.y_min() - 3.0).abs() < 1e-12);
    }
}
