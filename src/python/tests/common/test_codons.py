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

from ensembl.vep.common import codons as codons_module
from ensembl.vep.common.codons import (
    CODON_TABLE,
    DEFAULT_GENETIC_CODE_TABLE_ID,
    GENETIC_CODE_TABLE_IDS,
    STOP_CODONS,
    get_stop_codons,
    translate_dna,
)


class TestCodonTable:
    """NCBI genetic-code table construction."""

    @pytest.mark.parametrize(
        ("table_id", "amino_acids"),
        [
            (1, "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (2, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"),
            (3, "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (4, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (5, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"),
            (6, "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (9, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"),
            (10, "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (11, "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (12, "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (13, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"),
            (14, "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"),
            (15, "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (16, "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (21, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"),
            (22, "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (23, "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (24, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"),
            (25, "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (26, "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (27, "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (28, "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (29, "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (30, "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (31, "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (32, "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            (33, "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"),
        ],
    )
    def test_builds_ncbi_genetic_code(self, table_id: int, amino_acids: str) -> None:
        """Build each supported NCBI amino-acid row in NCBI codon order."""
        codon_table = codons_module._build_codon_table(table_id)  # pylint: disable=protected-access
        codons = (
            first_base + second_base + third_base
            for first_base in "TCAG"
            for second_base in "TCAG"
            for third_base in "TCAG"
        )

        assert "".join(codon_table[codon] for codon in codons) == amino_acids

    def test_supported_ids_match_current_ncbi_tables(self) -> None:
        """Expose exactly the current NCBI genetic-code table identifiers."""
        assert GENETIC_CODE_TABLE_IDS == {
            1,
            2,
            3,
            4,
            5,
            6,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            21,
            22,
            23,
            24,
            25,
            26,
            27,
            28,
            29,
            30,
            31,
            32,
            33,
        }

    def test_standard_table_is_the_default(self) -> None:
        """Keep the standard nuclear table as the module default."""
        assert DEFAULT_GENETIC_CODE_TABLE_ID == 1

    def test_unsupported_id_does_not_build_table(self) -> None:
        """Reject an identifier that NCBI does not assign to a table."""
        with pytest.raises(ValueError, match="Unsupported NCBI genetic-code table ID: 7"):
            codons_module._build_codon_table(7)  # pylint: disable=protected-access

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

    @pytest.mark.parametrize(
        ("table_id", "expected"),
        [
            pytest.param(1, {"TAA", "TAG", "TGA"}, id="standard"),
            pytest.param(2, {"TAA", "TAG", "AGA", "AGG"}, id="vertebrate-mitochondrial"),
            pytest.param(6, {"TGA"}, id="ciliate-nuclear"),
            pytest.param(27, set(), id="karyorelict-nuclear"),
        ],
    )
    def test_stop_codons_follow_genetic_code(self, table_id: int, expected: set[str]) -> None:
        """Derive stop codons from the selected NCBI table."""
        assert get_stop_codons(table_id) == expected

    def test_stop_codons_reject_unsupported_table(self) -> None:
        """Reject stop-codon lookup for an unsupported table."""
        with pytest.raises(ValueError, match="Unsupported NCBI genetic-code table ID: 7"):
            get_stop_codons(7)


class TestTranslateDNA:
    """Tests for translating coding-orientation DNA sequences."""

    def test_translates_complete_codons(self) -> None:
        """Test a sequence of complete codons translates in order."""
        assert translate_dna("ATGGCTTAA") == "MA*"

    def test_translates_with_selected_genetic_code(self) -> None:
        """Use the selected NCBI table for reassigned codons."""
        assert translate_dna("TGAATAAGAAGG", table_id=2) == "WM**"

    def test_unsupported_genetic_code_is_rejected(self) -> None:
        """Reject translation with an unsupported NCBI table identifier."""
        with pytest.raises(ValueError, match="Unsupported NCBI genetic-code table ID: 7"):
            translate_dna("ATG", table_id=7)

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

    def test_non_ascii_unknown_codon_translates_to_unknown(self) -> None:
        """Treat non-ASCII characters as unknown DNA bases."""
        assert translate_dna("ATG😀") == "MX"

    def test_translation_continues_after_stop_codon(self) -> None:
        """Test translation returns all codons following an internal stop."""
        assert translate_dna("TAAGCTATG") == "*AM"

    def test_python_fallback_preserves_native_semantics(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Use the table implementation when the optional native extension is unavailable."""
        monkeypatch.setattr(codons_module, "_translate_dna_native", None)

        assert codons_module.translate_dna("ATGANNatgA") == "MXXX"

    def test_stale_native_extension_uses_python_for_selected_table(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Use Python for nonstandard tables when the installed native binding only accepts DNA."""
        monkeypatch.setattr(codons_module, "_translate_dna_native", lambda sequence: sequence)
        monkeypatch.setattr(codons_module, "_NATIVE_SUPPORTS_TABLE_ID", False)

        assert codons_module.translate_dna("TGA", table_id=2) == "W"

    def test_native_extension_receives_selected_table(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Pass a nonstandard table identifier to a compatible native binding."""
        received: list[tuple[str, int]] = []

        def translate_native(sequence: str, table_id: int) -> str:
            received.append((sequence, table_id))
            return "native"

        monkeypatch.setattr(codons_module, "_translate_dna_native", translate_native)
        monkeypatch.setattr(codons_module, "_NATIVE_SUPPORTS_TABLE_ID", True)

        assert codons_module.translate_dna("TGA", table_id=2) == "native"
        assert received == [("TGA", 2)]
