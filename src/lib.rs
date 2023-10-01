#![doc = include_str!("../README.md")]

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
