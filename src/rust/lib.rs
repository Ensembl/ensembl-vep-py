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

// Keep this crate root as PyO3 module registration glue only. Rust logic should
// live in submodules where it can be tested directly; coverage jobs ignore this
// file with cargo-llvm-cov --ignore-filename-regex so stable Rust is sufficient.
pub mod hello;
pub mod parsers;
pub mod test_utils;
pub mod utils;

#[pymodule]
mod ensembl_vep_core {
    use pyo3::prelude::*;

    #[pymodule_init]
    fn init(_m: &Bound<'_, PyModule>) -> PyResult<()> {
        pyo3_log::init();
        Ok(())
    }

    #[pymodule_export]
    use crate::hello::sum_as_string;

    #[pymodule_export]
    use crate::parsers::tabix::TabixParser;

    #[pymodule_export]
    use crate::parsers::gff3::{Gff3Parser, Gff3Record};
}
