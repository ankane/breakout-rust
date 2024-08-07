use std::error;
use std::fmt;

/// An error.
#[derive(Debug, Eq, PartialEq)]
pub enum Error {
    Parameter(String),
}

impl error::Error for Error {}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Error::Parameter(ref err) => write!(f, "{}", err.as_str()),
        }
    }
}
