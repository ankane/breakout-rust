//! Breakout detection for Rust
//!
//! [View the docs](https://github.com/ankane/breakout-rust)

mod amoc;
mod edm_multi;
mod edm_tail;
mod edmx;
mod error;
mod multi;
mod multiset;

pub use amoc::{amoc, AmocParams};
pub use error::Error;
pub use multi::{multi, MultiParams};
