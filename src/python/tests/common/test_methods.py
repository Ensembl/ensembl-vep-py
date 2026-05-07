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
"""Unit tests for common methods."""

import pytest

from ensembl.vep.common import trim_sequences, shift_sequence, Strand, ClassSOTerm, infer_variant_so_term


class TestTrimSequences:
    """Tests for the trim_sequences function."""

    def test_trim_sequences_basic_left_right(self) -> None:
        """Trims matching bases from both ends."""
        ref, alt, start, end, changed = trim_sequences("ACGT", "ACCT", start=10)
        assert (ref, alt, start, end, changed) == ("G", "C", 12, 12, True)

    def test_trim_sequences_end_first(self) -> None:
        """Trims from right first, then left."""
        ref, alt, start, end, changed = trim_sequences("ACGT", "ACCT", start=10, end_first=True)
        assert (ref, alt, start, end, changed) == ("G", "C", 12, 12, True)

    def test_trim_sequences_empty_to_dash(self) -> None:
        """Empty alleles replaced with '-'."""
        ref, alt, start, end, changed = trim_sequences("A", "A", start=5, empty_to_dash=True)
        assert (ref, alt, start, end, changed) == ("-", "-", 6, 5, True)

    def test_trim_sequences_reverse_strand(self) -> None:
        """Trimming with strand = -1."""
        ref, alt, start, end, changed = trim_sequences("ACGT", "ACCT", start=10, strand=Strand.REVERSE)
        assert (ref, alt, start, end, changed) == ("G", "C", 11, 11, True)

    def test_trim_sequences_none_input(self) -> None:
        """Should raise ValueError if ref or alt is None."""
        with pytest.raises(ValueError):
            trim_sequences(None, "A")  # type: ignore[arg-type]
        with pytest.raises(ValueError):
            trim_sequences("A", None)  # type: ignore[arg-type]

    def test_trim_sequences_no_trimming(self) -> None:
        """No common bases, no trimming."""
        ref, alt, start, end, changed = trim_sequences("A", "T", start=1)
        assert (ref, alt, start, end, changed) == ("A", "T", 1, 1, False)


class TestShiftSequence:
    """Tests for the shift_sequence function."""

    def test_shift_sequence_3prime(self) -> None:
        """Shift to the right (3')."""
        seq = "A"
        post_seq = "AA"
        pre_seq = "TT"
        start, end = 5, 5
        shift, shifted_seq, new_start, new_end, _ = shift_sequence(
            seq, post_seq=post_seq, pre_seq=pre_seq, start=start, end=end
        )
        assert (shift, shifted_seq, new_start, new_end) == (2, "A", 7, 7)

    def test_shift_sequence_5prime(self) -> None:
        """Shift to the left (5')."""
        seq = "A"
        post_seq = "TT"
        pre_seq = "AA"
        start, end = 5, 5
        shift, shifted_seq, new_start, new_end, _ = shift_sequence(
            seq, post_seq=post_seq, pre_seq=pre_seq, start=start, end=end, reverse=True
        )

        assert (shift, shifted_seq, new_start, new_end) == (2, "A", 3, 3)

    def test_shift_sequence_no_shift(self) -> None:
        """No shift possible."""
        seq = "A"
        post_seq = "T"
        pre_seq = "G"
        start, end = 5, 5
        shift, shifted_seq, new_start, new_end, _ = shift_sequence(
            seq, post_seq=post_seq, pre_seq=pre_seq, start=start, end=end
        )

        assert (shift, shifted_seq, new_start, new_end) == (0, "A", 5, 5)

    def test_shift_sequence_full_shift(self) -> None:
        """Shift as much as possible."""
        seq = "A"
        post_seq = "AAAA"
        pre_seq = "AAAA"
        start, end = 5, 5
        shift, shifted_seq, new_start, new_end, _ = shift_sequence(
            seq, post_seq=post_seq, pre_seq=pre_seq, start=start, end=end
        )

        assert (shift, shifted_seq, new_start, new_end) == (4, "A", 9, 9)

        # 5' direction
        shift, shifted_seq, new_start, new_end, _ = shift_sequence(
            seq, post_seq=post_seq, pre_seq=pre_seq, start=start, end=end, reverse=True
        )

        assert (shift, shifted_seq, new_start, new_end) == (4, "A", 1, 1)

    def test_shift_sequence_empty_seq(self) -> None:
        """Edge case: empty sequence."""
        seq = ""
        post_seq = "A"
        pre_seq = "A"
        start, end = 1, 1
        shift, shifted_seq, new_start, new_end, _ = shift_sequence(
            seq, post_seq=post_seq, pre_seq=pre_seq, start=start, end=end
        )
        assert (shift, shifted_seq, new_start, new_end) == (0, "", 1, 1)

        seq = "A"
        post_seq = "A"
        pre_seq = ""
        start, end = 1, 1
        shift, shifted_seq, new_start, new_end, _ = shift_sequence(
            seq, post_seq=post_seq, pre_seq=pre_seq, start=start, end=end
        )

        assert (shift, shifted_seq, new_start, new_end) == (1, "A", 2, 2)

        seq = "A"
        post_seq = ""
        pre_seq = "A"
        start, end = 1, 1
        shift, shifted_seq, new_start, new_end, _ = shift_sequence(
            seq, post_seq=post_seq, pre_seq=pre_seq, start=start, end=end, reverse=True
        )

        assert (shift, shifted_seq, new_start, new_end) == (1, "A", 0, 0)

    def test_shift_sequence_tandem_repeat_3prime_no_seq_change(self) -> None:
        """Shift to the right (3') on a tandem repeat.
        Result should have changed position but no sequence change."""
        seq = "AT"
        post_seq = "ATATTT"
        pre_seq = "TTAT"
        start, end = 5, 6
        shift, shifted_seq, new_start, new_end, _ = shift_sequence(
            seq, post_seq=post_seq, pre_seq=pre_seq, start=start, end=end
        )
        assert (shift, shifted_seq, new_start, new_end) == (4, "AT", 9, 10)

    def test_shift_sequence_tandem_repeat_3prime_w_seq_change(self) -> None:
        """Shift to the right (3') on a tandem repeat.
        Result should have changed position and sequence change (shifted till end of sequence)."""
        seq = "TAGC"
        post_seq = "TAGCTAGCTATGCGA"
        pre_seq = "TAGC"
        start, end = 5, 4
        shift, shifted_seq, new_start, new_end, _ = shift_sequence(
            seq, post_seq=post_seq, pre_seq=pre_seq, start=start, end=end
        )
        assert (shift, shifted_seq, new_start, new_end) == (10, "GCTA", 15, 14)

    def test_shift_sequence_incomplete_dt_loop_limit(self) -> None:
        """Shift to the right (3') on a tandem repeat.
        Result should have changed position and sequence change but not to complete right.
        The post_seq is short and thus reach loop limit before reaching end of sequence."""
        seq = "TAGC"
        post_seq = "TAGCTAGCTA"
        pre_seq = "TAGC"
        start, end = 5, 4
        shift, shifted_seq, new_start, new_end, _ = shift_sequence(
            seq, post_seq=post_seq, pre_seq=pre_seq, start=start, end=end
        )
        assert (shift, shifted_seq, new_start, new_end) == (7, "CTAG", 12, 11)

    def test_shift_sequence_tandem_repeat_5prime_w_seq_change(self) -> None:
        """Shift to the left (5') on a tandem repeat."""
        seq = "AT"
        post_seq = "ATATTT"
        pre_seq = "TTAT"
        start, end = 5, 6
        shift, shifted_seq, new_start, new_end, _ = shift_sequence(
            seq, post_seq=post_seq, pre_seq=pre_seq, start=start, end=end, reverse=True
        )
        assert (shift, shifted_seq, new_start, new_end) == (3, "TA", 2, 3)

    def test_shift_sequence_none_input(self) -> None:
        """Should raise ValueError if seq, post_seq, or pre_seq is None."""
        with pytest.raises(ValueError):
            shift_sequence(None, post_seq="A", pre_seq="A", start=1, end=1)  # type: ignore[arg-type]
        with pytest.raises(ValueError):
            shift_sequence("A", post_seq=None, pre_seq="A", start=1, end=1)  # type: ignore[arg-type]
        with pytest.raises(ValueError):
            shift_sequence("A", post_seq="A", pre_seq=None, start=1, end=1)  # type: ignore[arg-type]


class TestInferVariantSOTerm:
    """Tests for the infer_variant_so_term function."""

    @pytest.mark.parametrize(
        "allele_string, expected",
        [
            ("A/T", ClassSOTerm.SNV),
            ("A/G", ClassSOTerm.SNV),
            ("A/T/G", ClassSOTerm.SNV),
            ("A/T/G/C", ClassSOTerm.SNV),
        ],
    )
    def test_snv(self, allele_string: str, expected: ClassSOTerm) -> None:
        """Test that a SNV variation returns the correct ClassSOTerm."""
        assert infer_variant_so_term(allele_string, ref_correct=True) == expected

    @pytest.mark.parametrize(
        "allele_string, expected",
        [
            ("(CAG)8/(CAG)9", ClassSOTerm.TANDEM_REPEAT),
            ("(AT)5/(AT)7", ClassSOTerm.TANDEM_REPEAT),
        ],
    )
    def test_tandem_repeat(self, allele_string: str, expected: ClassSOTerm) -> None:
        """Test that a tandem repeat variation returns the correct ClassSOTerm."""
        assert infer_variant_so_term(allele_string, ref_correct=True) == expected

    @pytest.mark.parametrize(
        "allele_string, expected",
        [
            ("-/ATT", ClassSOTerm.INSERTION),
            ("-/A", ClassSOTerm.INSERTION),
            ("A/ATT", ClassSOTerm.INDEL),
            ("A/AT", ClassSOTerm.INDEL),
            ("-/LARGEINSERTION", ClassSOTerm.INSERTION),
            ("(LARGEINSERTION)", ClassSOTerm.INSERTION),
        ],
    )
    def test_insertion(self, allele_string: str, expected: ClassSOTerm) -> None:
        """Test that an insertion variation returns the correct ClassSOTerm."""
        assert infer_variant_so_term(allele_string, ref_correct=True) == expected

    @pytest.mark.parametrize(
        "allele_string, expected",
        [
            ("A/-", ClassSOTerm.DELETION),
            ("ATT/-", ClassSOTerm.DELETION),
            ("AT/A", ClassSOTerm.INDEL),
            ("ATT/A", ClassSOTerm.INDEL),
            ("(LARGEDELETION)/-", ClassSOTerm.DELETION),
            ("(LARGEDELETION)", ClassSOTerm.DELETION),
            ("(2345 BP DELETION)/-", ClassSOTerm.DELETION),
        ],
    )
    def test_deletion(self, allele_string: str, expected: ClassSOTerm) -> None:
        """Test that a deletion variation returns the correct ClassSOTerm."""
        assert infer_variant_so_term(allele_string, ref_correct=True) == expected

    @pytest.mark.parametrize(
        "allele_string, expected",
        [
            ("AA/TT", ClassSOTerm.SUBSTITUTION),
            ("ACT/TCG", ClassSOTerm.SUBSTITUTION),
        ],
    )
    def test_substitution(self, allele_string: str, expected: ClassSOTerm) -> None:
        """Test that a substitution variation returns the correct ClassSOTerm."""
        assert infer_variant_so_term(allele_string, ref_correct=True) == expected

    @pytest.mark.parametrize(
        "allele_string, expected",
        [
            ("AA/TTTT", ClassSOTerm.INDEL),
            ("AAAA/TT", ClassSOTerm.INDEL),
        ],
    )
    def test_indel(self, allele_string: str, expected: ClassSOTerm) -> None:
        """Test that an indel variation returns the correct ClassSOTerm."""
        assert infer_variant_so_term(allele_string, ref_correct=True) == expected

    def test_empty_allele_string_defaults_to_sequence_alteration(self) -> None:
        """Test that an empty allele string returns the correct (default) ClassSOTerm."""
        assert infer_variant_so_term("", ref_correct=True) == ClassSOTerm.SEQUENCE_ALTERATION
