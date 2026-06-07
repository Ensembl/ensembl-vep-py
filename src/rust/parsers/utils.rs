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

use noodles_core::Region;
use std::io::{self, ErrorKind};

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
}
