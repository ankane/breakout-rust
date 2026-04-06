#[cfg(feature = "no_std")]
pub use libm::{ceil, log as ln, pow as powf};

#[cfg(not(feature = "no_std"))]
pub fn ceil(x: f64) -> f64 {
    x.ceil()
}

#[cfg(not(feature = "no_std"))]
pub fn ln(x: f64) -> f64 {
    x.ln()
}

pub fn pow2(x: f64) -> f64 {
    x * x
}

#[cfg(not(feature = "no_std"))]
pub fn powf(x: f64, n: f64) -> f64 {
    x.powf(n)
}
