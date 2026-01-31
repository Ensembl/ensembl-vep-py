// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

use pyo3::prelude::*;

// plain Rust helper, not exposed to Python: useful for threading and unit test
fn sum_as_string_impl(a: &str, b: &str) -> String {
    a.to_owned() + b
}

#[pymodule]
mod _ensembl_vep_core {
    use super::sum_as_string_impl;
    use pyo3::prelude::*;

    /// Can we do this more elegantly such as PyModule.add_function?
    #[pyfunction]
    fn sum_as_string(a: &str, b: &str) -> PyResult<String> {
        Ok(sum_as_string_impl(a, b))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sum_as_string_impl() {
        assert_eq!(sum_as_string_impl("1", "2"), "12");
        assert_eq!(sum_as_string_impl("foo", "bar"), "foobar");
        assert_eq!(sum_as_string_impl("", "x"), "x");
    }
}
