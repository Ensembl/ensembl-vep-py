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

//! NCBI genetic-code translation kernel.

use std::{error::Error, fmt};

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

const GENETIC_CODE_TABLES: &[(u8, &[u8; 64])] = &[
    (
        1,
        b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        2,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    ),
    (
        3,
        b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        4,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        5,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    ),
    (
        6,
        b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        9,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    ),
    (
        10,
        b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        11,
        b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        12,
        b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        13,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    ),
    (
        14,
        b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    ),
    (
        15,
        b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        16,
        b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        21,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    ),
    (
        22,
        b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        23,
        b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        24,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
    ),
    (
        25,
        b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        26,
        b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        27,
        b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        28,
        b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        29,
        b"FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        30,
        b"FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        31,
        b"FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        32,
        b"FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        33,
        b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
    ),
];
const DEFAULT_GENETIC_CODE_TABLE_ID: u8 = 1;
const UNKNOWN: u8 = b'X';

/// Error returned for an unsupported NCBI genetic-code table identifier.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct UnsupportedGeneticCodeError {
    /// Unsupported table identifier.
    pub table_id: u8,
}

impl fmt::Display for UnsupportedGeneticCodeError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "unsupported NCBI genetic-code table ID: {}",
            self.table_id
        )
    }
}

impl Error for UnsupportedGeneticCodeError {}

/// Return an NCBI genetic-code table in NCBI codon order.
///
/// The amino-acid rows match NCBI `ncbieaa`. Alternative initiation codons are represented by
/// NCBI's separate `sncbieaa` rows and do not change translation within a coding sequence.
pub fn build_codon_table(table_id: u8) -> Result<&'static [u8; 64], UnsupportedGeneticCodeError> {
    GENETIC_CODE_TABLES
        .iter()
        .find_map(|(candidate_id, amino_acids)| (*candidate_id == table_id).then_some(*amino_acids))
        .ok_or(UnsupportedGeneticCodeError { table_id })
}

/// Return the table index for one uppercase standard DNA base.
#[inline]
fn base_index(base: u8) -> Option<usize> {
    match base {
        b'T' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// Translate one uppercase DNA triplet, or return the unknown residue.
#[inline]
fn translate_codon(codon_table: &[u8; 64], first: u8, second: u8, third: u8) -> u8 {
    match (base_index(first), base_index(second), base_index(third)) {
        (Some(first), Some(second), Some(third)) => codon_table[first * 16 + second * 4 + third],
        _ => UNKNOWN,
    }
}

/// Translate coding-orientation DNA with the standard nuclear codon table.
///
/// Only uppercase `T`, `C`, `A`, and `G` are standard bases. Ambiguous, lowercase, and incomplete
/// codons each yield `X`; translation continues after stop codons. Non-ASCII input has no standard
/// codons and is counted by Unicode scalar values, matching Python string slicing.
pub fn translate_dna(coding_sequence: &str) -> String {
    translate_dna_with_table(coding_sequence, DEFAULT_GENETIC_CODE_TABLE_ID)
        .expect("the default NCBI genetic-code table is embedded")
}

/// Translate coding-orientation DNA with an NCBI genetic-code table.
///
/// Only uppercase `T`, `C`, `A`, and `G` are standard bases. Ambiguous, lowercase, and incomplete
/// codons each yield `X`; translation continues after stop codons. Non-ASCII input has no standard
/// codons and is counted by Unicode scalar values, matching Python string slicing.
pub fn translate_dna_with_table(
    coding_sequence: &str,
    table_id: u8,
) -> Result<String, UnsupportedGeneticCodeError> {
    let codon_table = build_codon_table(table_id)?;

    if !coding_sequence.is_ascii() {
        let characters: Vec<char> = coding_sequence.chars().collect();
        let mut peptide = String::with_capacity(characters.len().div_ceil(3));
        for codon in characters.chunks(3) {
            if codon.len() == 3 && codon.iter().all(|base| base.is_ascii()) {
                peptide.push(translate_codon(
                    codon_table,
                    codon[0] as u8,
                    codon[1] as u8,
                    codon[2] as u8,
                ) as char);
            } else {
                peptide.push(UNKNOWN as char);
            }
        }
        return Ok(peptide);
    }

    let sequence = coding_sequence.as_bytes();
    let mut peptide = Vec::with_capacity(sequence.len().div_ceil(3));
    let mut position = 0;
    while position + 2 < sequence.len() {
        peptide.push(translate_codon(
            codon_table,
            sequence[position],
            sequence[position + 1],
            sequence[position + 2],
        ));
        position += 3;
    }
    if position < sequence.len() {
        peptide.push(UNKNOWN);
    }
    // Every byte is an ASCII amino-acid symbol inserted above.
    Ok(String::from_utf8(peptide).expect("DNA translation emits ASCII residues"))
}

/// Expose DNA translation to Python with table 1 as the default.
#[pyfunction(name = "translate_dna", signature = (coding_sequence, table_id = 1))]
pub fn translate_dna_python(coding_sequence: &str, table_id: i64) -> PyResult<String> {
    let table_id = u8::try_from(table_id).map_err(|_| {
        PyValueError::new_err(format!(
            "Unsupported NCBI genetic-code table ID: {table_id}"
        ))
    })?;
    translate_dna_with_table(coding_sequence, table_id).map_err(|error| {
        PyValueError::new_err(format!(
            "Unsupported NCBI genetic-code table ID: {}",
            error.table_id
        ))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_current_ncbi_genetic_code_builds() {
        let expected_ids = [
            1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29,
            30, 31, 32, 33,
        ];
        let actual_ids: Vec<u8> = GENETIC_CODE_TABLES
            .iter()
            .map(|(table_id, _)| *table_id)
            .collect();

        assert_eq!(actual_ids, expected_ids);
        for table_id in expected_ids {
            assert!(build_codon_table(table_id).is_ok());
        }
    }

    #[test]
    fn standard_codons_translate_in_table_order() {
        assert_eq!(translate_dna("ATGGCTTAA"), "MA*");
    }

    #[test]
    fn selected_genetic_code_changes_translation() {
        assert_eq!(
            translate_dna_with_table("TGAATAAGAAGG", 2),
            Ok("WM**".into())
        );
    }

    #[test]
    fn python_translation_accepts_and_validates_table_ids() {
        Python::initialize();

        assert_eq!(translate_dna_python("TGA", 2).unwrap(), "W");
        for table_id in [-1, 7] {
            let error = translate_dna_python("ATG", table_id).unwrap_err();
            Python::attach(|python| assert!(error.is_instance_of::<PyValueError>(python)));
        }
    }

    #[test]
    fn unsupported_genetic_code_returns_error() {
        let error = UnsupportedGeneticCodeError { table_id: 7 };

        assert_eq!(build_codon_table(7), Err(error));
        assert_eq!(translate_dna_with_table("ATG", 7), Err(error));
        assert_eq!(
            error.to_string(),
            "unsupported NCBI genetic-code table ID: 7"
        );
    }

    #[test]
    fn invalid_and_incomplete_codons_translate_to_unknown() {
        assert_eq!(translate_dna("ATGANNatgA"), "MXXX");
    }

    #[test]
    fn non_ascii_input_uses_python_character_groups() {
        assert_eq!(translate_dna("ATG😀"), "MX");
    }

    #[test]
    fn translation_keeps_codon_after_stop() {
        assert_eq!(translate_dna("TAAGCTATG"), "*AM");
    }
}
