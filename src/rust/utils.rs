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

use pyo3::{prelude::*, types::PyDict};

#[repr(i8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Strand {
    Forward = 1,
    Reverse = -1,
    Unknown = 0,
}


impl<'py> IntoPyObject<'py> for Strand {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        let module = py.import("ensembl_vep_core")?;
        let strand_cls = module.getattr("Strand")?;

        strand_cls.call1((self as i8,))
    }
}

#[repr(i8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExonPhase {
    Zero = 0,
    One = 1,
    Two = 2,
    Unknown = -1,
}

impl<'py> IntoPyObject<'py> for ExonPhase {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        let module = py.import("ensembl_vep_core")?;
        let exon_phase_cls = module.getattr("ExonPhase")?;

        exon_phase_cls.call1((self as i8,))
    }
}

pub(crate) fn register_enums(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    let enum_mod = py.import("enum")?;

    // Register Strand enum.
    let members = PyDict::new(py);
    members.set_item("FORWARD", Strand::Forward as i8)?;
    members.set_item("REVERSE", Strand::Reverse as i8)?;
    members.set_item("UNKNOWN", Strand::Unknown as i8)?;
    let strand = enum_mod.getattr("IntEnum")?.call1(("Strand", members))?;
    m.add("Strand", strand)?;

    // Register ExonPhase enum.
    let members = PyDict::new(py);
    members.set_item("ZERO", ExonPhase::Zero as i8)?;
    members.set_item("ONE", ExonPhase::One as i8)?;
    members.set_item("TWO", ExonPhase::Two as i8)?;
    members.set_item("UNKNOWN", ExonPhase::Unknown as i8)?;
    let exon_phase = enum_mod.getattr("IntEnum")?.call1(("ExonPhase", members))?;
    m.add("ExonPhase", exon_phase)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strand_enum() {
        assert_eq!(Strand::Forward as i8, 1);
        assert_eq!(Strand::Reverse as i8, -1);
        assert_eq!(Strand::Unknown as i8, 0);
    }

    #[test]
    fn test_exon_phase_enum() {
        assert_eq!(ExonPhase::Zero as i8, 0);
        assert_eq!(ExonPhase::One as i8, 1);
        assert_eq!(ExonPhase::Two as i8, 2);
        assert_eq!(ExonPhase::Unknown as i8, -1);
    }

    #[test]
    fn test_register_enums_strand() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "ensembl_vep_core").unwrap();
            register_enums(&module).unwrap();

            // Test Strand enum registration.
            let actual_strand = module.getattr("Strand").unwrap();
            assert_eq!(
                actual_strand
                    .getattr("FORWARD")
                    .unwrap()
                    .extract::<i8>()
                    .unwrap(),
                Strand::Forward as i8
            );
            assert_eq!(
                actual_strand
                    .getattr("REVERSE")
                    .unwrap()
                    .extract::<i8>()
                    .unwrap(),
                Strand::Reverse as i8
            );
            assert_eq!(
                actual_strand
                    .getattr("UNKNOWN")
                    .unwrap()
                    .extract::<i8>()
                    .unwrap(),
                Strand::Unknown as i8
            );
        });
    }

    #[test]
    fn test_register_enums_exon_phase() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "ensembl_vep_core").unwrap();
            register_enums(&module).unwrap();

            // Test ExonPhase enum registration.
            let actual_exon_phase = module.getattr("ExonPhase").unwrap();
            assert_eq!(
                actual_exon_phase
                    .getattr("ZERO")
                    .unwrap()
                    .extract::<i8>()
                    .unwrap(),
                ExonPhase::Zero as i8
            );
            assert_eq!(
                actual_exon_phase
                    .getattr("ONE")
                    .unwrap()
                    .extract::<i8>()
                    .unwrap(),
                ExonPhase::One as i8
            );
            assert_eq!(
                actual_exon_phase
                    .getattr("TWO")
                    .unwrap()
                    .extract::<i8>()
                    .unwrap(),
                ExonPhase::Two as i8
            );
            assert_eq!(
                actual_exon_phase
                    .getattr("UNKNOWN")
                    .unwrap()
                    .extract::<i8>()
                    .unwrap(),
                ExonPhase::Unknown as i8
            );
        });
    }
}
