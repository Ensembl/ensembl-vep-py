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

use std::io::{self, ErrorKind};

use noodles_core::Region;
use pyo3::{exceptions::PyValueError, PyResult};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ReaderMode {
    Auto,
    Sequential,
    Query,
}

impl ReaderMode {
    pub(crate) fn parse(value: Option<&str>) -> PyResult<Self> {
        match value.unwrap_or("auto") {
            "auto" => Ok(Self::Auto),
            "sequential" => Ok(Self::Sequential),
            "query" => Ok(Self::Query),
            value => Err(PyValueError::new_err(format!(
                "unsupported reader mode {value:?}: expected 'auto', 'query' or 'sequential'"
            ))),
        }
    }

    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Sequential => "sequential",
            Self::Query => "query",
        }
    }
}

pub(crate) fn parse_region(region: &str) -> io::Result<Region> {
    region.parse().map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidInput,
            format!("invalid tabix region {region:?}: {e}"),
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_region_valid() {
        assert!(parse_region("chr1:1-1000").is_ok());
    }

    #[test]
    fn test_parse_region_invalid() {
        let err = parse_region("not_a_region::??").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput);
        assert!(err.to_string().contains("not_a_region::??"));
    }

    #[test]
    fn test_parse_region_empty() {
        let err = parse_region("").unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput);
        assert!(err.to_string().contains("invalid tabix region"));
    }

    #[test]
    fn test_reader_mode_parse_defaults_to_auto() {
        assert_eq!(ReaderMode::parse(None).unwrap(), ReaderMode::Auto);
    }

    #[test]
    fn test_reader_mode_parse_accepts_aliases() {
        assert_eq!(ReaderMode::parse(Some("query")).unwrap(), ReaderMode::Query);
        assert_eq!(ReaderMode::parse(Some("sequential")).unwrap(), ReaderMode::Sequential);
        assert_eq!(ReaderMode::parse(Some("auto")).unwrap(), ReaderMode::Auto);
    }

    #[test]
    fn test_reader_mode_parse_rejects_invalid_mode() {
        let err = ReaderMode::parse(Some("invalid-mode")).unwrap_err();
        assert!(
            err.to_string()
                .contains("unsupported reader mode \"invalid-mode\": expected 'auto', 'query' or 'sequential'")
        );
    }

    #[test]
    fn test_reader_mode_as_str() {
        assert_eq!(ReaderMode::Sequential.as_str(), "sequential");
        assert_eq!(ReaderMode::Query.as_str(), "query");
        assert_eq!(ReaderMode::Auto.as_str(), "auto");
    }
}