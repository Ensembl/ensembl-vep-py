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

# pylint: disable=redefined-outer-name
"""Unit tests for the VariationAllele class."""

from __future__ import annotations

import logging
import pytest
from attrs import evolve

from ensembl.vep.common import CoordSystemLevel, Strand
from ensembl.vep.models.variation import Variation, VariationAllele
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import OverlapConsequence, TranscriptVariationAllele
from ensembl.vep.models.sequence import CoordSystem, Slice
from tests.models.conftest import DummySequenceSource

logger = logging.getLogger(__name__)


@pytest.fixture
def genomic_slice(dummy_sequence_source: DummySequenceSource) -> Slice:
    """Create a slice object for testing."""
    dummy_sequence_source.data = {"1": "ATCGATCGATCGG"}
    return Slice(
        coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
        seq_region_name="1",
        start=1,
        end=13,
        strand=Strand.FORWARD,
        sequence_source=dummy_sequence_source,
    )


@pytest.fixture
def variation(genomic_slice: Slice) -> Variation:
    """Returns a Variation object for testing."""
    return Variation(
        stable_id="rs123456",
        start=5,
        end=5,
        strand=Strand.FORWARD,
        slice=genomic_slice,
        allele_string="A/T",
    )


@pytest.fixture
def variation_allele(genomic_slice: Slice) -> VariationAllele:
    """Returns a VariationAllele object for testing."""
    return VariationAllele(
        start=5, end=5, strand=Strand.FORWARD, slice=genomic_slice, reference_allele="A", alternate_allele="T"
    )


@pytest.fixture
def trancript_variation_allele(variation_allele: VariationAllele) -> TranscriptVariationAllele:
    """Returns a TranscriptVariationAllele object for testing."""
    return TranscriptVariationAllele(
        feature=Transcript(start=1, end=100, strand=Strand.FORWARD, stable_id="ENST00000000001"),
        overlapping_feature=variation_allele,
    )


class TestVariationAlleleProperties:
    """Tests property accessors and defaults."""

    def test_initialization(self, variation_allele: VariationAllele) -> None:
        """Test that a VariationAllele can be initialized with the correct attributes."""
        assert variation_allele.start == 5
        assert variation_allele.end == 5
        assert variation_allele.strand == Strand.FORWARD
        assert variation_allele.reference_allele == "A"
        assert variation_allele.alternate_allele == "T"
        assert variation_allele.variation is None

    def test_initialization_without_keywords(self) -> None:
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            VariationAllele(
                5,  # type: ignore[arg-type, call-arg, misc]
                end=5,
                strand=Strand.FORWARD,
                slice=None,
                reference_allele="A",
                alternate_allele="T",
            )

    def test_variation_can_be_assigned(self, variation_allele: VariationAllele, variation: Variation) -> None:
        """Test that an allele can reference its owning variation."""
        variation_allele.variation = variation

        assert variation_allele.variation is variation

    def test_minimised(self, variation_allele: VariationAllele) -> None:
        """Test the minimised property before and after calling minimise."""
        assert variation_allele.minimised is False
        variation_allele.minimise()
        assert variation_allele.minimised is True

        # Check that minimised_locus is set
        assert variation_allele.minimised_locus is not None
        assert variation_allele.original_locus is not None

        # Check that minimised is immutable by external code
        with pytest.raises(
            AttributeError, match="property 'minimised' of 'VariationAllele' object has no setter"
        ):
            variation_allele.minimised = False  # type: ignore

    def test_shifted(self, variation_allele: VariationAllele) -> None:
        """Test the shifted property before and after calling shift."""
        assert variation_allele.shifted is False
        variation_allele.shift()
        assert variation_allele.shifted is True

        # Check that shifted_locus is set
        assert variation_allele.shifted_locus is not None
        assert variation_allele.original_locus is not None

        # Check that shifted is immutable by external code
        with pytest.raises(
            AttributeError, match="property 'shifted' of 'VariationAllele' object has no setter"
        ):
            variation_allele.shifted = False  # type: ignore

    def test_overlap_consequences(
        self,
        variation_allele: VariationAllele,
        trancript_variation_allele: TranscriptVariationAllele,
    ) -> None:
        """Test overlap_consequences property."""
        variation_allele.transcript_variation_alleles = []
        assert variation_allele.overlap_consequences == []

        tva1 = trancript_variation_allele
        tva1.overlap_consequences = [OverlapConsequence(term="missense_variant", accession="SO:0001583")]
        tva2 = evolve(trancript_variation_allele)
        tva2.overlap_consequences = [OverlapConsequence(term="synonymous_variant", accession="SO:0001819")]
        variation_allele.transcript_variation_alleles += [tva1, tva2]
        assert variation_allele.overlap_consequences == [
            OverlapConsequence(term="missense_variant", accession="SO:0001583"),
            OverlapConsequence(term="synonymous_variant", accession="SO:0001819"),
        ]


class TestVariationAlleleMinimise:
    """Tests for the minimise method of VariationAllele."""

    def test_minimise_raises_if_missing_alleles(self, variation_allele: VariationAllele) -> None:
        """Test minimise raises if reference or alternate variation_allele is None."""
        a = VariationAllele(
            start=1,
            end=1,
            strand=Strand.FORWARD,
            slice=variation_allele.slice,
            reference_allele=None,
            alternate_allele="T",
        )
        with pytest.raises(ValueError, match="Missing reference or alternate allele for minimisation"):
            a.minimise()
        b = VariationAllele(
            start=1,
            end=1,
            strand=Strand.FORWARD,
            slice=variation_allele.slice,
            reference_allele="A",
            alternate_allele=None,
        )
        with pytest.raises(ValueError, match="Missing reference or alternate allele for minimisation"):
            b.minimise()

    def test_minimise_sets_locus_correctly(self, variation_allele: VariationAllele) -> None:
        """Test that minimise sets minimised_locus and original_locus correctly."""
        variation_allele.minimise()
        assert variation_allele.minimised_locus is not None
        assert variation_allele.minimised_locus.start == variation_allele.start
        assert variation_allele.minimised_locus.end == variation_allele.end
        assert (
            variation_allele.minimised_locus.allele_string
            == f"{variation_allele.reference_allele}/{variation_allele.alternate_allele}"
        )
        assert variation_allele.original_locus is not None

    def test_minimise_forced(
        self, variation_allele: VariationAllele, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Test that calling minimise with force=True re-minimises even if already minimised."""
        variation_allele.minimise()
        original_locus = variation_allele.original_locus
        minimised_locus = variation_allele.minimised_locus

        with caplog.at_level(logging.DEBUG):
            variation_allele.minimise()
            assert "Variation allele is already minimised. Use force=True to re-minimise." in caplog.text

        variation_allele.minimise(force=True)
        assert variation_allele.minimised is True
        assert variation_allele.original_locus == original_locus
        assert variation_allele.minimised_locus == minimised_locus


class TestVariationAlleleShift:
    """Tests for the shift method of VariationAllele."""

    def test_shift_raises_if_missing_alleles(self, variation_allele: VariationAllele) -> None:
        """Test shift raises if reference or alternate variation_allele is None."""
        a = VariationAllele(
            start=1,
            end=1,
            strand=Strand.FORWARD,
            slice=variation_allele.slice,
            reference_allele=None,
            alternate_allele="T",
        )
        with pytest.raises(ValueError, match="Missing reference or alternate allele for shifting"):
            a.shift()
        b = VariationAllele(
            start=1,
            end=1,
            strand=Strand.FORWARD,
            slice=variation_allele.slice,
            reference_allele="A",
            alternate_allele=None,
        )
        with pytest.raises(ValueError, match="Missing reference or alternate allele for shifting"):
            b.shift()

    def test_shift_sets_locus_correctly(self, variation_allele: VariationAllele) -> None:
        """Test that shift sets shifted_locus and original_locus correctly."""
        variation_allele.shift()
        assert variation_allele.shifted_locus is not None
        assert variation_allele.shifted_locus.start == variation_allele.start
        assert variation_allele.shifted_locus.end == variation_allele.end
        assert (
            variation_allele.shifted_locus.allele_string
            == f"{variation_allele.reference_allele}/{variation_allele.alternate_allele}"
        )
        assert variation_allele.original_locus is not None


class TestVariationAlleleMinimiseAndShift:
    """Tests for interaction between minimise and shift methods."""

    def test_minimise_independent(self, variation_allele: VariationAllele) -> None:
        """Test that minimise can be called independently and respective flags are set."""
        variation_allele.minimise()
        assert variation_allele.minimised is True
        assert variation_allele.shifted is False

    def test_shift_independent(self, variation_allele: VariationAllele) -> None:
        """Test that shift can be called independently and respective flags are set."""
        variation_allele.shift()
        assert variation_allele.shifted is True
        assert variation_allele.minimised is False

    def test_minimise_no_effect_then_shift_no_effect(self, variation_allele: VariationAllele) -> None:
        """Test that calling shift after minimise retains the original locus and sets flags correctly.
        In this testcase:
         * Minimise should not change anything
         * Shift should not change anything
         => Shift should not reset the minimised state."""
        variation_allele.minimise()
        assert variation_allele.minimised is True, "minimised was not set on minimise"
        original_locus = variation_allele.original_locus
        assert original_locus is not None, "original_locus was not set on minimise"
        minimised_start = variation_allele.start
        minimised_end = variation_allele.end
        minimised_alt_seq = variation_allele.alternate_allele
        minimised_ref_seq = variation_allele.reference_allele

        variation_allele.shift()
        assert variation_allele.shifted is True, "shifted was not set on shift"
        assert minimised_start == variation_allele.start, "start was changed by shift"
        assert minimised_end == variation_allele.end, "end was changed by shift"
        assert minimised_alt_seq == variation_allele.alternate_allele, "alt allele was changed by shift"
        assert minimised_ref_seq == variation_allele.reference_allele, "ref allele was changed by shift"

        # minimise should still be True as shift did not change positions nor alt/ref sequences
        assert variation_allele.minimised is True, "minimised was reset on shift"
        assert variation_allele.original_locus == original_locus, "original_locus was changed on shift"

    def test_minimise_no_effect_then_shift_wo_seq_change(self, genomic_slice: Slice) -> None:
        """Test that calling shift after minimise retains the original locus and sets flags correctly.
        In this testcase:
         * Minimise should not change anything
         * Shift should change the position but keep alt/ref sequences
         => Shift should not reset the minimised state (as alt/ref sequences were not changed)."""
        variation_allele2 = VariationAllele(
            start=12,
            end=12,
            strand=Strand.FORWARD,
            slice=genomic_slice,
            reference_allele="G",
            alternate_allele="T",
        )
        variation_allele2.minimise()
        assert variation_allele2.minimised is True, "minimised was not set on minimise"
        original_locus = variation_allele2.original_locus
        assert original_locus is not None, "original_locus was not set on minimise"
        minimised_locus = variation_allele2.minimised_locus
        assert minimised_locus is not None, "minimised_locus returns None after minimise"
        assert minimised_locus.start == original_locus.start, "start was changed by minimise"
        assert minimised_locus.end == original_locus.end, "end was changed by minimise"
        assert (
            minimised_locus.allele_string == original_locus.allele_string
        ), "ref/alt allele string was changed by minimise"
        minimised_alt_seq = variation_allele2.alternate_allele
        minimised_ref_seq = variation_allele2.reference_allele

        variation_allele2.shift()
        assert variation_allele2.shifted is True, "shifted was not set on shift"
        assert minimised_locus.start != variation_allele2.start, "start was not changed by shift"
        assert minimised_locus.end != variation_allele2.end, "end was not changed by shift"
        assert minimised_alt_seq == variation_allele2.alternate_allele, "alt allele was changed by shift"
        assert minimised_ref_seq == variation_allele2.reference_allele, "ref allele was changed by shift"

        # minimise should still be True as shift did not change alt/ref sequences, only positions
        assert variation_allele2.minimised is True, "minimised was reset on shift"
        assert variation_allele2.original_locus == original_locus, "original_locus was changed on shift"

    @pytest.mark.parametrize(
        (
            "slice_start",
            "slice_seq",
            "variation_allele_start",
            "variation_allele_end",
            "ref_allele",
            "alt_allele",
        ),
        [
            (10, "TAGCTAGCTAGCTA", 5, 4, "", "TAGC"),
            (20, "ATATCTCTACT", 4, 3, "", "TC"),
        ],
    )
    def test_minimise_no_effect_then_shift_w_seq_change(
        self,
        slice_start: int,
        slice_seq: str,
        variation_allele_start: int,
        variation_allele_end: int,
        ref_allele: str,
        alt_allele: str,
        dummy_sequence_source: DummySequenceSource,
    ) -> None:
        """Test that calling shift after minimise retains the original locus and sets flags correctly.
        In this testcase:
         * Minimise should not change anything
         * Shift should change the position and ref/alt allele sequence
         => Shift should reset the minimised state (as ref/alt sequences was changed)."""
        dummy_sequence_source.data = {"FAKE_CHROMOSOME": ("N" * (slice_start - 1)) + slice_seq}
        repeat_slice = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="FAKE_CHROMOSOME",
            start=slice_start,
            end=slice_start + len(slice_seq) - 1,
            strand=Strand.FORWARD,
            sequence_source=dummy_sequence_source,
        )
        variation_allele2 = VariationAllele(
            start=variation_allele_start,
            end=variation_allele_end,
            strand=Strand.FORWARD,
            slice=repeat_slice,
            reference_allele=ref_allele,
            alternate_allele=alt_allele,
        )
        variation_allele2.minimise()
        assert variation_allele2.minimised is True, "minimised was not set on minimise"
        original_locus = variation_allele2.original_locus
        assert original_locus is not None, "original_locus was not set on minimise"
        minimised_locus = variation_allele2.minimised_locus
        assert minimised_locus is not None, "minimised_locus returns None after minimise"
        assert minimised_locus.start == original_locus.start, "start was changed by minimise"
        assert minimised_locus.end == original_locus.end, "end was changed by minimise"
        assert (
            minimised_locus.allele_string == original_locus.allele_string
        ), "ref/alt allele string was changed by minimise"

        variation_allele2.shift()
        assert variation_allele2.shifted is True, "shifted was not set on shift"
        assert variation_allele2.original_locus == original_locus, "original_locus was changed on shift"
        shifted_locus = variation_allele2.shifted_locus
        assert shifted_locus is not None, "shifted_locus returns None after shift"
        assert shifted_locus.start != minimised_locus.start, "start was not changed by shift"
        assert shifted_locus.end != minimised_locus.end, "end was not changed by shift"
        assert (
            shifted_locus.allele_string != minimised_locus.allele_string
        ), "ref/alt allele string was not changed by shift"
        assert variation_allele2.minimised is False, "minimised was not reset on shift"

    def test_shift_no_effect_then_minimise_no_effect(self, variation_allele: VariationAllele) -> None:
        """Test that calling minimise after shift retains the original locus and sets flags correctly.
        In this testcase:
         * Shift should not change anything
         * Minimise should not change anything
         => Minimise should not reset the shifted state."""
        assert variation_allele.shifted is False

        variation_allele.shift()
        assert variation_allele.shifted is True
        original_locus = variation_allele.original_locus

        variation_allele.minimise()
        assert variation_allele.minimised is True
        # shift should still be True as minimise does not change anything
        assert variation_allele.shifted is True
        assert variation_allele.original_locus == original_locus

    def test_shift_no_effect_then_minimise_w_seq_change(self, genomic_slice: Slice) -> None:
        """Test that calling minimise after shift retains the original locus and sets flags correctly.
        In this testcase:
         * Shift should not change anything
         * Minimise should change the Allele seqs from GG/G to G/-
         => Minimise should reset the shifted state, as it affected ref/alt seq."""
        variation_allele2 = VariationAllele(
            start=12,
            end=13,
            strand=Strand.FORWARD,
            slice=genomic_slice,
            reference_allele="GG",
            alternate_allele="G",
        )
        assert variation_allele2.shifted is False, "shifted was not initialised as False"
        variation_allele2.shift()
        assert variation_allele2.shifted is True, "shifted was not set on shift"
        original_locus = variation_allele2.original_locus
        assert original_locus is not None, "original_locus was not set on shift"
        shifted_locus = variation_allele2.shifted_locus
        assert shifted_locus is not None, "shifted_locus returns None after shift"
        assert original_locus.start == shifted_locus.start, "start was changed by shift"
        assert original_locus.end == shifted_locus.end, "end was changed by shift"
        assert (
            original_locus.allele_string == shifted_locus.allele_string
        ), "ref/alt allele string was changed by shift"

        variation_allele2.minimise()
        assert variation_allele2.minimised is True, "minimised was not set on minimise"
        minimised_locus = variation_allele2.minimised_locus
        assert minimised_locus is not None, "minimised_locus returns None after minimise"
        assert minimised_locus.start != shifted_locus.start, "start was not changed by minimise"
        assert minimised_locus.end == shifted_locus.end, "end was changed by minimise"
        assert (
            minimised_locus.allele_string != shifted_locus.allele_string
        ), "ref/alt allele string was not changed by minimise"

        # shift should be False as minimise changes coordinate
        # and we do not know if it is still shifted
        assert variation_allele2.shifted is False, "shifted was not reset on minimise"
        assert variation_allele2.original_locus == original_locus, "original_locus was changed by minimise"


class TestVariationAlleleGetAllFeatureOverlaps:
    """Tests for the get_all_feature_overlaps method of VariationAllele."""

    def test_get_all_feature_overlaps_no_overlaps(self, variation_allele: VariationAllele) -> None:
        """Test that get_all_feature_overlaps returns empty list if no overlaps."""
        assert variation_allele.get_all_featureoverlaps() == []

    def test_get_all_feature_overlaps_with_transcript_variation_alleles(
        self,
        variation_allele: VariationAllele,
        trancript_variation_allele: TranscriptVariationAllele,
    ) -> None:
        """Test that get_all_feature_overlaps returns transcript variation alleles."""
        variation_allele.transcript_variation_alleles = [trancript_variation_allele]
        overlaps = variation_allele.get_all_featureoverlaps()

        assert len(overlaps) == 1
        # assert overlaps[0] == trancript_variation_allele

        # Needs to add __eq__ for TranscriptVariationAllele


class TestDisplayName:
    """The name an allele is identified by."""

    def test_an_allele_with_no_variant_takes_its_locus_name(self, variation_allele: VariationAllele) -> None:
        """The locus name is the shape a variant selector is parsed from."""

        assert variation_allele.display_name == "1:5:A>T"

    def test_a_single_allele_record_takes_its_plain_identifier(self, variation: Variation) -> None:
        """One alternate allele needs nothing appended to stay unique."""

        assert variation.variation_alleles[0].display_name == "rs123456"

    def test_a_multi_allelic_record_appends_the_alternate_allele(self, variation: Variation) -> None:
        """One identifier covers every allele, so each name carries its own alternate."""

        multi_allelic = evolve(variation, allele_string="A/T/G")

        names = [allele.display_name for allele in multi_allelic.variation_alleles]

        assert names == ["rs123456:T", "rs123456:G"]

    def test_a_record_without_an_identifier_falls_back_to_the_locus(self, variation: Variation) -> None:
        """A locus is unique per allele already, so no counter is needed."""

        unnamed = evolve(variation, stable_id=None)

        assert unnamed.variation_alleles[0].display_name == "1:5:A>T"

    def test_an_empty_allele_is_spelled_as_a_dash(self, genomic_slice: Slice) -> None:
        """A deletion's empty side is spelled the way a selector spells it."""

        deletion = VariationAllele(
            start=5,
            end=5,
            strand=Strand.FORWARD,
            slice=genomic_slice,
            reference_allele="A",
            alternate_allele="",
        )

        assert deletion.display_name == "1:5:A>-"

    def test_a_missing_sequence_region_degrades(self) -> None:
        """An allele with no slice still names itself rather than raising."""

        unplaced = VariationAllele(
            start=5, end=5, strand=Strand.FORWARD, reference_allele="A", alternate_allele="T"
        )

        assert unplaced.display_name == "?:5:A>T"
