// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

use pyo3::{
    exceptions::{PyIOError, PyRuntimeError, PyValueError},
    prelude::*,
};
use std::io::{self, ErrorKind};

pub fn trim_newline(line: &mut String) {
    if line.ends_with('\n') {
        line.pop();
        if line.ends_with('\r') {
            line.pop();
        }
    }
}

pub fn to_py_error(err: io::Error) -> PyErr {
    match err.kind() {
        ErrorKind::InvalidData | ErrorKind::InvalidInput => PyValueError::new_err(err.to_string()),
        ErrorKind::NotFound => PyIOError::new_err(err.to_string()),
        _ => PyRuntimeError::new_err(err.to_string()),
    }
}

pub fn from_py_error(err: PyErr) -> io::Error {
    io::Error::other(err.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // from_py_error
    // -----------------------------------------------------------------------

    #[test]
    fn test_from_py_error() {
        let py_err = PyValueError::new_err("value error occurred");
        let io_err = from_py_error(py_err);
        assert_eq!(io_err.to_string(), "ValueError: value error occurred");
    }

    // -----------------------------------------------------------------------
    // to_py_error — error kind mapping
    // -----------------------------------------------------------------------

    #[test]
    fn test_to_py_error_invalid_data_maps_to_value_error() {
        let err = io::Error::new(ErrorKind::InvalidData, "bad data");
        let py_err = to_py_error(err);

        Python::initialize();
        let type_name = Python::attach(|py| py_err.get_type(py).name().unwrap().to_string());

        assert_eq!(type_name, "ValueError");
    }

    #[test]
    fn test_to_py_error_not_found_maps_to_io_error() {
        let err = io::Error::new(ErrorKind::NotFound, "not found");
        let py_err = to_py_error(err);

        Python::initialize();
        let type_name = Python::attach(|py| py_err.get_type(py).name().unwrap().to_string());

        assert_eq!(type_name, "OSError");
    }

    #[test]
    fn test_to_py_error_kind_routing() {
        // Verify some expected error kinds are mapped without panicking.
        for kind in [
            ErrorKind::InvalidData,
            ErrorKind::InvalidInput,
            ErrorKind::NotFound,
            ErrorKind::Other,
            ErrorKind::PermissionDenied,
        ] {
            let _ = to_py_error(io::Error::new(kind, "test"));
        }
    }

    // -----------------------------------------------------------------------
    // trim_newline
    // -----------------------------------------------------------------------

    #[test]
    fn test_trim_newline_lf() {
        let mut line = "hello\n".to_string();
        trim_newline(&mut line);
        assert_eq!(line, "hello");
    }

    #[test]
    fn test_trim_newline_crlf() {
        let mut line = "hello\r\n".to_string();
        trim_newline(&mut line);
        assert_eq!(line, "hello");
    }

    #[test]
    fn test_trim_newline_no_newline() {
        let mut line = "hello".to_string();
        trim_newline(&mut line);
        assert_eq!(line, "hello");
    }

    #[test]
    fn test_trim_newline_empty() {
        let mut line = String::new();
        trim_newline(&mut line);
        assert_eq!(line, "");
    }
}
