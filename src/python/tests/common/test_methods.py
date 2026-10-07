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

from functools import cached_property
from typing import Any

import pytest
from attrs import define, field

from ensembl.vep.common import (
    ClassSOTerm,
    Strand,
    infer_variant_so_term,
    infer_symbolic_variant_so_term,
    format_region_string,
    parse_region_string,
    invalidate,
    invalidates,
    is_resolved,
    resets_attrs,
    reverse_complement,
    shift_sequence,
    trim_sequences,
)


@define(kw_only=True, auto_attribs=True)
class _ResolvedProbe:
    """Small slots-backed fixture for is_resolved tests."""

    start: int
    compute_calls: int = 0

    @cached_property
    def _derived(self) -> int:
        """Cached property used to probe resolved vs unresolved state."""
        object.__setattr__(self, "compute_calls", self.compute_calls + 1)
        return self.start + 1

    @property
    def derived(self) -> int:
        """Public accessor for _derived."""
        return self._derived


@define(kw_only=True, auto_attribs=True, slots=False)
class _CacheInvalidationProbe:
    """Probe cached-property invalidation and attribute reset behavior."""

    start: int = field(on_setattr=invalidates("_length"))
    end: int = field(on_setattr=resets_attrs(("_minimised", False)))
    _minimised: bool = field(default=False)

    @cached_property
    def _length(self) -> int:
        """Calculate the inclusive length."""
        return self.end - self.start + 1

    @cached_property
    def _range(self) -> str:
        """Calculate the range notation."""
        return f"{self.start}..{self.end}"

    @property
    def length(self) -> int:
        """Return the inclusive length."""
        return self._length

    @property
    def minimised(self) -> bool:
        """Return whether the probe has been minimised."""
        return self._minimised

    def minimise(self) -> None:
        """Mark the probe as minimised."""
        self._minimised = True

    @property
    def range(self) -> str:
        """Return the range notation."""
        return self._range


class _InheritedCacheInvalidationProbe(_CacheInvalidationProbe):
    """Probe cached-property invalidation inherited from a parent class."""


class TestIsResolved:
    """Tests for the is_resolved function."""

    def test_unresolved_does_not_trigger_computation(self) -> None:
        """is_resolved should report unresolved state without computing the value."""

        probe = _ResolvedProbe(start=3)
        assert probe.compute_calls == 0
        assert is_resolved(probe, "_derived") is False
        assert probe.compute_calls == 0

    def test_resolved_after_first_access(self) -> None:
        """is_resolved should report resolved state after the cached property is computed."""

        probe = _ResolvedProbe(start=3)
        assert probe.derived == 4
        assert probe.compute_calls == 1
        assert is_resolved(probe, "_derived") is True
        assert probe.compute_calls == 1

    def test_raises_for_non_descriptor_name(self) -> None:
        """is_resolved should raise for names that are not descriptors on the class."""

        probe = _ResolvedProbe(start=3)
        with pytest.raises(AttributeError):
            is_resolved(probe, "_does_not_exist")


class TestCacheInvalidation:
    """Tests for the cached property invalidation functions."""

    def test_on_setattr_invalidation(self) -> None:
        """Tests the invalidates function invalidates the named properties on dependent attributes set."""

        # Evaluate cached_property attributes and values on initialization
        probe = _CacheInvalidationProbe(start=1, end=10)
        assert (
            "_length" not in probe.__dict__
        ), "Protected length attribute was initialized before first attribute call"
        assert (
            "_range" not in probe.__dict__
        ), "Protected range attribute was initialized before first attribute call"
        assert probe.length == 10, "Unexpected initial length after first attribute call"
        assert probe.range == "1..10", "Unexpected initial range after first attribute call"
        assert "_length" in probe.__dict__, "Protected length attribute not cached after first attribute call"
        assert "_range" in probe.__dict__, "Protected range attribute not cached after first attribute call"
        assert (
            probe._length == 10  # pylint: disable=protected-access
        ), "Unexpected cached protected length attribute value"
        assert (
            probe._range == "1..10"  # pylint: disable=protected-access
        ), "Unexpected cached protected range attribute value"

        # Evaluate attributes invalidation through on_setattr
        probe.start = 2
        assert "_length" not in probe.__dict__, "Protected length attribute was not invalidated."
        assert probe.length == 9, "Updated length should be 9"
        assert (
            probe._length == 9  # pylint: disable=protected-access
        ), "Protected length attribute was not reloaded with new value after invalidation"
        assert (
            probe._range == "1..10"  # pylint: disable=protected-access
        ), "Unexpected change to protected range attribute through on_setattr"

    def test_direct_invalidation(self) -> None:
        """Tests the invalidate function invalidates the named properties directly."""

        probe = _CacheInvalidationProbe(start=1, end=10)
        assert probe.length == 10, "Unexpected initial length"
        assert probe.range == "1..10", "Unexpected initial range"

        invalidate(probe, "_range")
        assert (
            "_length" in probe.__dict__ and probe._length == 10  # pylint: disable=protected-access
        ), "Protected length attribute should not be affected by range invalidation."
        assert "_range" not in probe.__dict__, "Protected range attribute did not invalidate"

    def test_on_setattr_inheritance(self) -> None:
        """Tests the invalidates function invalidates the named properties
        on dependent attributes set of inherited attributes."""

        inherited_probe = _InheritedCacheInvalidationProbe(start=1, end=10)
        assert inherited_probe.length == 10, "Unexpected initial length value"
        assert inherited_probe.range == "1..10", "Unexpected initial range value"
        assert (
            "_length" in inherited_probe.__dict__
        ), "Protected length attribute not cached after first attribute call"
        assert (
            inherited_probe._length == 10  # pylint: disable=protected-access
        ), "Unexpected (cached) value in protected length attribute"
        inherited_probe.start = 2
        assert (
            "_length" not in inherited_probe.__dict__
        ), "Protected length attribute did not invalidate on dependent attribute change"
        assert (
            inherited_probe._length == 9  # pylint: disable=protected-access
        ), "Unexpected new value in protected length attribute after invalidation"
        assert (
            inherited_probe._range == "1..10"  # pylint: disable=protected-access
        ), "Protected range attribute should not change on invalidation of length"


class TestAttributeReset:
    """Tests for the attribute reset functions."""

    def test_on_setattr_reset(self) -> None:
        """Tests dependent attribute updates correctly reset the named properties."""

        # Evaluate attributes and values on initialization
        probe = _CacheInvalidationProbe(start=1, end=10)
        assert "_minimised" in probe.__dict__, "Protected minimised attribute was not initialized"
        assert (
            probe._minimised is False  # pylint: disable=protected-access
        ), "Unexpected initial minimised value"

        probe.minimise()
        assert (
            probe._minimised is True  # pylint: disable=protected-access
        ), "Minimised attribute not updated on minimise"

        # trigger reset and evaluate attributes and values
        probe.end = 20
        assert (
            "_minimised" in probe.__dict__
        ), "Protected minimised attribute was invalidated rather than reset"
        assert (
            probe._minimised is False  # pylint: disable=protected-access
        ), "Minimised attribute did not reset on end update"

    def test_on_setattr_inheritance(self) -> None:
        """Tests dependent attribute updates correctly reset the named inherited properties"""

        # Evaluate attributes and values on initialization
        inherited_probe = _InheritedCacheInvalidationProbe(start=1, end=10)
        assert "_minimised" in inherited_probe.__dict__, "Protected minimised attribute was not initialized"
        assert (
            inherited_probe._minimised is False  # pylint: disable=protected-access
        ), "Unexpected initial minimised value"

        inherited_probe.minimise()
        assert (
            inherited_probe._minimised is True  # pylint: disable=protected-access
        ), "Minimised attribute not updated on minimise"

        # trigger reset and evaluate attributes and values
        inherited_probe.end = 20
        assert (
            "_minimised" in inherited_probe.__dict__
        ), "Protected minimised attribute was invalidated rather than reset"
        assert (
            inherited_probe._minimised is False  # pylint: disable=protected-access
        ), "Minimised attribute did not reset on end update"


class TestTrimSequences:
    """Tests for the trim_sequences function."""

    def test_trim_sequences_unknown_strand(self) -> None:
        """Should return untrimmed alleles if strand is unknown."""
        ref, alt, start, end, changed = trim_sequences("ACGT", "ACCT", start=10, strand=Strand.UNKNOWN)
        assert (ref, alt, start, end, changed) == ("ACGT", "ACCT", 10, 13, False)

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


class TestParseRegionString:
    """Tests for the parse_region_string function."""

    @pytest.mark.parametrize(
        ("region", "expected"),
        [
            pytest.param("chr1:10-20", ("chr1", 10, 20), id="name-start-end"),
            pytest.param("chr1:10", ("chr1", 10, None), id="name-start"),
            pytest.param("chr1", ("chr1", None, None), id="name-only"),
            pytest.param("HLA:DRB1:10-20", ("HLA:DRB1", 10, 20), id="name-with-colon"),
        ],
    )
    def test_splits_a_region(self, region: str, expected: tuple[str, int | None, int | None]) -> None:
        """Returns the sequence-region name with its optional bounds."""

        assert parse_region_string(region) == expected

    def test_round_trips_with_format_region_string(self) -> None:
        """Parses back what format_region_string builds."""

        assert parse_region_string(format_region_string(seq_region_name="1", start=5, end=9)) == ("1", 5, 9)

    @pytest.mark.parametrize(
        ("region", "error_match"),
        [
            pytest.param("", None, id="empty"),
            pytest.param("chr1:x-20", None, id="not-a-number"),
            pytest.param("chr1:0-20", "1-based", id="not-positive"),
        ],
    )
    def test_rejects_an_unparseable_region(self, region: str, error_match: str | None) -> None:
        """Raises when the region is empty or its coordinates are not integers."""

        with pytest.raises(ValueError, match=error_match):
            parse_region_string(region)


class TestFormatRegionString:
    """Tests for the format_region_string function."""

    def test_existing_region(self) -> None:
        """Returns an explicit region unchanged."""

        assert format_region_string("chr1:10-20") == "chr1:10-20"

    def test_seq_region_name_only(self) -> None:
        """Formats a whole-sequence region."""

        assert format_region_string(seq_region_name="chr1") == "chr1"

    def test_seq_region_name_start_end(self) -> None:
        """Formats a bounded 1-based region."""

        assert format_region_string(seq_region_name="chr1", start=10, end=20) == "chr1:10-20"

    def test_seq_region_name_start(self) -> None:
        """Formats a region with only start coordinate."""

        assert format_region_string(seq_region_name="chr1", start=10) == "chr1:10"

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"region": "chr1:10-20", "seq_region_name": "chr1", "start": None, "end": None},
            {"region": None, "seq_region_name": None, "start": None, "end": None},
            {"region": None, "seq_region_name": "chr1", "start": None, "end": 10},
            {"region": None, "seq_region_name": "chr1", "start": 0, "end": 20},
            {"region": None, "seq_region_name": "chr1", "start": 10, "end": 0},
        ],
    )
    def test_rejects_invalid_input(self, kwargs: dict[str, Any]) -> None:
        """Rejects ambiguous or invalid region inputs."""

        with pytest.raises(ValueError):
            format_region_string(**kwargs)


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


class TestInferSymbolicVariantSOTerm:
    """Tests for the infer_symbolic_variant_so_term function."""

    @pytest.mark.parametrize(
        "allele_string, expected",
        [
            ("<INS>", ClassSOTerm.INSERTION),
            ("<DEL>", ClassSOTerm.DELETION),
            ("<INV>", ClassSOTerm.INVERSION),
            ("<DUP>", ClassSOTerm.DUPLICATION),
            ("BND", ClassSOTerm.CHROMOSOME_BREAKPOINT),
        ],
    )
    def test_basic_symbolic_types(self, allele_string: str, expected: ClassSOTerm) -> None:
        """Map basic symbolic types to their SO terms."""
        assert infer_symbolic_variant_so_term(allele_string) == expected

    @pytest.mark.parametrize(
        "allele_string, expected",
        [
            ("<INS:ME>", ClassSOTerm.MOBILE_ELEMENT_INSERTION),
            ("<INS:ME:ALU>", ClassSOTerm.ALU_INSERTION),
            ("<INS:ME:HERV>", ClassSOTerm.HERV_INSERTION),
            ("<INS:ME:L1>", ClassSOTerm.LINE1_INSERTION),
            ("<INS:ME:SVA>", ClassSOTerm.SVA_INSERTION),
            ("<DEL:ME>", ClassSOTerm.MOBILE_ELEMENT_DELETION),
            ("<DEL:ME:ALU>", ClassSOTerm.ALU_DELETION),
            ("<DEL:ME:HERV>", ClassSOTerm.HERV_DELETION),
            ("<DEL:ME:L1>", ClassSOTerm.LINE1_DELETION),
            ("<DEL:ME:SVA>", ClassSOTerm.SVA_DELETION),
            ("<INS:ME:OTHER>", ClassSOTerm.MOBILE_ELEMENT_INSERTION),
        ],
    )
    def test_mobile_element_types(self, allele_string: str, expected: ClassSOTerm) -> None:
        """Map generic and recognized mobile element types."""
        assert infer_symbolic_variant_so_term(allele_string) == expected

    @pytest.mark.parametrize(
        "allele_string, expected",
        [
            ("<DEL>/<DUP>", ClassSOTerm.COPY_NUMBER_VARIATION),
            ("<DUP:TANDEM>", ClassSOTerm.TANDEM_DUPLICATION),
            ("<CNV:TR>", ClassSOTerm.TANDEM_REPEAT),
            ("<CN0>", ClassSOTerm.DELETION),
            ("<CN=0>", ClassSOTerm.DELETION),
            ("<CN2>", ClassSOTerm.DUPLICATION),
            ("<CN=2>", ClassSOTerm.DUPLICATION),
            ("<CN1>", ClassSOTerm.COPY_NUMBER_VARIATION),
            ("<CN0>,<CN2>", ClassSOTerm.COPY_NUMBER_VARIATION),
            ("<CNV>", ClassSOTerm.COPY_NUMBER_VARIATION),
        ],
    )
    def test_copy_number_and_tandem_types(self, allele_string: str, expected: ClassSOTerm) -> None:
        """Map copy-number and tandem symbolic types."""
        assert infer_symbolic_variant_so_term(allele_string) == expected

    @pytest.mark.parametrize(
        "allele_string",
        [
            "N]2:321682]",
            "[2:321682[N",
            ".N",
            "N.",
        ],
    )
    def test_breakend_types(self, allele_string: str) -> None:
        """Recognize bracket and single-breakend notation."""
        assert infer_symbolic_variant_so_term(allele_string) == ClassSOTerm.CHROMOSOME_BREAKPOINT

    @pytest.mark.parametrize("allele_string", ["INS/INS", "<INS>,<INS>"])
    def test_repeated_types_are_normalized(self, allele_string: str) -> None:
        """Normalize repeated representations of the same symbolic type."""
        assert infer_symbolic_variant_so_term(allele_string) == ClassSOTerm.INSERTION

    def test_symbolic_allele_takes_precedence_over_svtype(self) -> None:
        """Retain the more specific subtype from a VCF-compliant ALT allele."""
        assert infer_symbolic_variant_so_term("<INS:ME:ALU>", svtype="INS") == ClassSOTerm.ALU_INSERTION

    @pytest.mark.parametrize("allele_string", ["", "<CPX>"])
    def test_svtype_is_used_for_missing_or_non_standard_allele(self, allele_string: str) -> None:
        """Use INFO/SVTYPE when ALT does not contain a supported VCF symbolic prefix."""
        assert infer_symbolic_variant_so_term(allele_string, svtype="DEL") == ClassSOTerm.DELETION

    @pytest.mark.parametrize("allele_string", ["", "<CPX>", "TRA"])
    def test_unsupported_type_returns_none(self, allele_string: str) -> None:
        """Return None when no supported structural variant type can be inferred."""
        assert infer_symbolic_variant_so_term(allele_string) is None


class TestReverseComplement:
    """Tests for reverse_complement."""

    def test_reverse_complement(self) -> None:
        """Reverse-complement a simple nucleotide sequence."""

        assert reverse_complement("ACGTN") == "NACGT"
        assert reverse_complement("acgtn") == "nacgt"
