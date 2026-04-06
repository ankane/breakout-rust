#![doc = include_str!("../README.md")]
#![cfg_attr(feature = "no_std", no_std)]
#![allow(clippy::comparison_chain)]
#![allow(clippy::needless_range_loop)]

extern crate alloc;

mod amoc;
mod edm_multi;
mod edm_tail;
mod edmx;
mod error;
mod math;
mod multi;
mod multiset;

pub use amoc::{amoc, AmocParams};
pub use error::Error;
pub use multi::{multi, MultiParams};
