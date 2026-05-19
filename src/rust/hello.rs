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

use log::info;
use pyo3::prelude::*;

#[pyfunction]
pub fn sum_as_string(a: &str, b: &str) -> String {
    info!("Adding the numbers - {} and {}", a, b);
    a.to_owned() + b
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sum_as_string() {
        assert_eq!(sum_as_string("1", "2"), "12");
        assert_eq!(sum_as_string("foo", "bar"), "foobar");
        assert_eq!(sum_as_string("", "x"), "x");
    }
}
