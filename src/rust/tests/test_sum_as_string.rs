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

#[pymodule]
mod _ensembl_vep_core {
    use pyo3::prelude::*;

    /// Concatenates two strings and returns the result.
    #[pyfunction]
    fn sum_as_string(a: &str, b: &str) -> PyResult<String> {
        Ok((a.to_owned() + b).to_string())
    }
}

#[cfg(test)]
mod py_integration_tests {
    use super::*;
    use pyo3::prelude::*;

    #[test]
    fn test_sum_as_string_python() {
        // Register this module for the embedded Python interpreter
        pyo3::append_to_inittab!(_ensembl_vep_core);

        // Prepare Python for use from multiple threads (required in tests)
        pyo3::prepare_freethreaded_python();

        Python::with_gil(|py| {
            // Import the module as Python would
            let module = py
                .import_bound("_ensembl_vep_core")
                .expect("Failed to import _ensembl_vep_core");

            // Get the Python-exposed function
            let func = module
                .getattr("sum_as_string")
                .expect("Failed to get sum_as_string");

            // Call it with arguments
            let result: String = func
                .call1(("1", "2"))
                .expect("Python call failed")
                .extract()
                .expect("Failed to extract String");

            assert_eq!(result, "12");
        });
    }
}
