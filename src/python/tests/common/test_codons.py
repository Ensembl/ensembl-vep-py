# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Unit tests for codon utilities."""

import pytest

from ensembl.vep.common.codons import CODON_TABLE, STOP_CODONS, translate_dna


class TestCodonTable:
    """Tests for the standard nuclear codon table."""

    def test_contains_every_standard_codon(self) -> None:
        """Test the table contains all 64 unique DNA codons."""
        expected_codons = {
            first_base + second_base + third_base
            for first_base in "ACGT"
            for second_base in "ACGT"
            for third_base in "ACGT"
        }

        assert set(CODON_TABLE) == expected_codons

    @pytest.mark.parametrize(
        ("codon", "amino_acid"),
        [
            ("ATG", "M"),
            ("TGG", "W"),
            ("TTT", "F"),
            ("GCT", "A"),
        ],
    )
    def test_maps_codon_to_amino_acid(self, codon: str, amino_acid: str) -> None:
        """Test representative codons map to their expected amino acids."""
        assert CODON_TABLE[codon] == amino_acid

    def test_stop_codons_match_codon_table(self) -> None:
        """Test the stop-codon collection contains exactly the standard stops."""
        assert STOP_CODONS == {"TAA", "TAG", "TGA"}
        assert all(CODON_TABLE[codon] == "*" for codon in STOP_CODONS)


class TestTranslateDNA:
    """Tests for translating coding-orientation DNA sequences."""

    def test_translates_complete_codons(self) -> None:
        """Test a sequence of complete codons translates in order."""
        assert translate_dna("ATGGCTTAA") == "MA*"

    def test_empty_sequence_returns_empty_peptide(self) -> None:
        """Test an empty DNA sequence produces an empty peptide."""
        assert translate_dna("") == ""

    @pytest.mark.parametrize("partial_codon", ["A", "AT"])
    def test_incomplete_final_codon_translates_to_unknown(self, partial_codon: str) -> None:
        """Test a trailing partial codon translates to an unknown residue."""
        assert translate_dna(f"ATG{partial_codon}") == "MX"

    @pytest.mark.parametrize("unknown_codon", ["ANN", "atg", "---"])
    def test_unknown_codon_translates_to_unknown(self, unknown_codon: str) -> None:
        """Test non-standard codons translate to an unknown residue."""
        assert translate_dna(unknown_codon) == "X"

    def test_translation_continues_after_stop_codon(self) -> None:
        """Test translation returns all codons following an internal stop."""
        assert translate_dna("TAAGCTATG") == "*AM"
