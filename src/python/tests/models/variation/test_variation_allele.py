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
"""Unit tests for the VariationAllele class."""

from __future__ import annotations

import logging
import pytest
from attrs import evolve

from ensembl.vep.common import Strand
from ensembl.vep.models.variation import VariationAllele
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import OverlapConsequence, TranscriptVariationAllele
from ensembl.vep.models.sequence import Slice

logger = logging.getLogger(__name__)


@pytest.fixture
def slice() -> Slice:
    """Create a slice object for testing."""
    return Slice(
        coord_system="chromosome",
        seq_region_name="1",
        seq_region_length=13,
        start=1,
        end=13,
        strand=Strand.FORWARD,
        seq="ATCGATCGATCGG",
    )


@pytest.fixture
def variation_allele(slice: Slice) -> VariationAllele:
    """Returns a VariationAllele object for testing."""
    return VariationAllele(start=5, end=5, strand=Strand.FORWARD, slice=slice, reference_allele="A", alternate_allele="T")


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
        assert variation_allele.strand == 1
        assert variation_allele.reference_allele == "A"
        assert variation_allele.alternate_allele == "T"

    def test_initialization_without_keywords(self) -> None:
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            VariationAllele(start=5, end=5, strand=Strand.FORWARD, slice=None, reference_allele="A", alternate_allele="T")

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
            variation_allele.minimised = False

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
            variation_allele.shifted = False

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
        assert variation_allele.shifted_locus.start == variation_allele.start
        assert variation_allele.shifted_locus.end == variation_allele.end
        assert (
            variation_allele.shifted_locus.allele_string
            == f"{variation_allele.reference_allele}/{variation_allele.alternate_allele}"
        )
        assert variation_allele.original_locus is not None


class TestVariationAlleleMinimiseAndShift:
    """Tests for interaction between minimise and shift methods."""

    def test_minimise_and_shift_independent(self, variation_allele: VariationAllele) -> None:
        """Test that minimise and shift can be called independently and set their respective flags."""
        variation_allele.minimise()
        assert variation_allele.minimised is True
        assert variation_allele.shifted is False
        variation_allele2 = VariationAllele(
            start=5, end=5, strand=Strand.FORWARD, slice=variation_allele.slice, reference_allele="A", alternate_allele="T"
        )
        variation_allele2.shift()
        assert variation_allele2.shifted is True
        assert variation_allele2.minimised is False

    def test_minimise_then_shift(self, variation_allele: VariationAllele) -> None:
        """Test that calling shift after minimise retains the original locus and sets flags correctly."""
        variation_allele.minimise()
        original_locus = variation_allele.original_locus
        variation_allele.shift()
        assert variation_allele.shifted is True
        # minimise should still be True as shift does not change minimised status
        assert variation_allele.minimised is True
        assert variation_allele.original_locus == original_locus

        variation_allele2 = VariationAllele(
            start=12,
            end=12,
            strand=Strand.FORWARD,
            slice=variation_allele.slice,
            reference_allele="G",
            alternate_allele="T",
        )
        variation_allele2.minimise()
        original_locus2 = variation_allele2.original_locus
        variation_allele2.shift()
        assert variation_allele2.shifted is True
        # minimise should should be False as shift changes coordinate and we do not know if it is still minimised
        assert variation_allele2.minimised is False
        assert variation_allele2.original_locus == original_locus2

    def test_shift_then_minimise(self, variation_allele: VariationAllele) -> None:
        """Test that calling minimise after shift retains the original locus and sets flags correctly."""
        variation_allele.shift()
        original_locus = variation_allele.original_locus
        variation_allele.minimise()
        assert variation_allele.minimised is True
        # shift should still be True as minimise does not change anything
        assert variation_allele.shifted is True
        assert variation_allele.original_locus == original_locus

        variation_allele2 = VariationAllele(
            start=12,
            end=12,
            strand=Strand.FORWARD,
            slice=variation_allele.slice,
            reference_allele="GG",
            alternate_allele="G",
        )
        variation_allele2.shift()
        original_locus2 = variation_allele2.original_locus
        variation_allele2.minimise()
        assert variation_allele2.minimised is True
        # shift should should be False as minimise changes coordinate and we do not know if it is still shifted
        assert variation_allele2.shifted is False
        assert variation_allele2.original_locus == original_locus2


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
