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
"""Unit tests for the StructuralVariationAllele class."""

from __future__ import annotations

import logging
import pytest
from attrs import evolve

from ensembl.vep.common import CoordSystemLevel, Strand
from ensembl.vep.models.variation import StructuralVariationAllele
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import OverlapConsequence, TranscriptStructuralVariationAllele
from ensembl.vep.models.sequence import CoordSystem, Slice

logger = logging.getLogger(__name__)


@pytest.fixture
def genomic_slice() -> Slice:
    """Create a slice object for testing."""
    return Slice(
        coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
        seq_region_name="1",
        seq_region_length=1000,
        start=1,
        end=1000,
        strand=Strand.FORWARD,
    )


@pytest.fixture
def structural_variation_allele(genomic_slice: Slice) -> StructuralVariationAllele:
    """Returns a StructuralVariationAllele object for testing."""
    return StructuralVariationAllele(
        start=5,
        end=150,
        strand=Strand.FORWARD,
        slice=genomic_slice,
        reference_allele="-",
        alternate_allele="DEL",
    )


@pytest.fixture
def trancript_structural_variation_allele(
    structural_variation_allele: StructuralVariationAllele,
) -> TranscriptStructuralVariationAllele:
    """Returns a TranscriptStructuralVariationAllele object for testing."""
    return TranscriptStructuralVariationAllele(
        feature=Transcript(start=1, end=200, strand=Strand.FORWARD, stable_id="ENST00000000001"),
        overlapping_feature=structural_variation_allele,
    )


class TestStructuralVariationAlleleProperties:
    """Tests property accessors and defaults."""

    def test_initialization(self, structural_variation_allele: StructuralVariationAllele) -> None:
        """Test that a StructuralVariationAllele can be initialized with the correct attributes."""
        assert structural_variation_allele.start == 5, "start attribute initialized incorrectly"
        assert structural_variation_allele.end == 150, "end attribute initialized incorrectly"
        assert (
            structural_variation_allele.strand == Strand.FORWARD
        ), "strand attribute initialized incorrectly"
        assert (
            structural_variation_allele.reference_allele == "-"
        ), "reference_allele attribute initialized incorrectly"
        assert (
            structural_variation_allele.alternate_allele == "DEL"
        ), "alternate_allele attribute initialized incorrectly"

    def test_initialization_without_keywords(self) -> None:
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            StructuralVariationAllele(
                5,  # type: ignore[arg-type, call-arg, misc]
                end=150,
                strand=Strand.FORWARD,
                slice=None,
                reference_allele="-",
                alternate_allele="DEL",
            )

    def test_minimise_warning(
        self, structural_variation_allele: StructuralVariationAllele, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Test the minimise method emits a warning."""
        # check warning from minimise
        with caplog.at_level(logging.WARNING):
            structural_variation_allele.minimise()
            assert (
                "Structural Variation alleles are currently not minimised." in caplog.text
            ), "Incorrect warning on minimise"

        # check extra warning from minimise for force argument
        with caplog.at_level(logging.WARNING):
            structural_variation_allele.minimise(force=True)
            assert (
                "force=True doesn't change this behaviour." in caplog.text
            ), "Incorrect warning on minimise(force=True)"

    def test_minimised(self, structural_variation_allele: StructuralVariationAllele) -> None:
        """Test the minimised property before and after calling minimise."""
        assert structural_variation_allele.minimised is False, "minimised attribute incorrectly initialized"
        structural_variation_allele.minimise()
        assert (
            structural_variation_allele.minimised is True
        ), "minimised attribute not set on call to minimise"

        # Check that minimised_locus is set
        assert (
            structural_variation_allele.minimised_locus is not None
        ), "minimised_locus not set on call to minimise"
        assert (
            structural_variation_allele.original_locus is not None
        ), "original_locus not set on call to minimise"

        # Check that minimised is immutable by external code
        with pytest.raises(
            AttributeError, match="property 'minimised' of 'StructuralVariationAllele' object has no setter"
        ):
            structural_variation_allele.minimised = False  # type: ignore

    def test_shift_warning(
        self, structural_variation_allele: StructuralVariationAllele, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Test the shift method emits a warning."""
        with caplog.at_level(logging.WARNING):
            structural_variation_allele.shift()
            assert (
                "Structural Variation alleles are currently not shifted." in caplog.text
            ), "Incorrect warning on shift"

        # check extra warning from minimise for force argument
        with caplog.at_level(logging.WARNING):
            structural_variation_allele.shift(force=True)
            assert (
                "force=True doesn't change this behaviour." in caplog.text
            ), "Incorrect warning on shift(force=True)"

    def test_shifted(self, structural_variation_allele: StructuralVariationAllele) -> None:
        """Test the shifted property before and after calling shift."""
        assert structural_variation_allele.shifted is False, "shifted attribute incorrectly initialized"
        structural_variation_allele.shift()
        assert structural_variation_allele.shifted is True, "shifted attribute not set on call to minimise"

        # Check that shifted_locus is set
        assert (
            structural_variation_allele.shifted_locus is not None
        ), "shifted_locus not set on call to minimise"
        assert (
            structural_variation_allele.original_locus is not None
        ), "original_locus not set on call to shift"

        # Check that shifted is immutable by external code
        with pytest.raises(
            AttributeError, match="property 'shifted' of 'StructuralVariationAllele' object has no setter"
        ):
            structural_variation_allele.shifted = False  # type: ignore


class TestStructuralVariationAlleleGetAllFeatureOverlaps:
    """Tests for the get_all_feature_overlaps method of StructuralVariationAllele."""

    def test_get_all_feature_overlaps_no_overlaps(
        self, structural_variation_allele: StructuralVariationAllele
    ) -> None:
        """Test that get_all_feature_overlaps returns empty list if no overlaps."""
        assert (
            structural_variation_allele.get_all_featureoverlaps() == []
        ), "Incorrect return value from get_all_featureoverlaps"

    def test_get_all_feature_overlaps_with_transcript_structural_variation_alleles(
        self,
        structural_variation_allele: StructuralVariationAllele,
        trancript_structural_variation_allele: TranscriptStructuralVariationAllele,
    ) -> None:
        """Test that get_all_feature_overlaps returns transcript variation alleles."""
        structural_variation_allele.transcript_variation_alleles = [trancript_structural_variation_allele]
        overlaps = structural_variation_allele.get_all_featureoverlaps()

        assert len(overlaps) == 1, "List of incorrect length returned from get_all_featureoverlaps"

        # Needs to add __eq__ for TranscriptStructuralVariationAllele


class TestStructuralVariationAlleleOverlapConsequences:
    """Tests for the overlap_consequences property"""

    def test_overlap_consequences(
        self,
        structural_variation_allele: StructuralVariationAllele,
        trancript_structural_variation_allele: TranscriptStructuralVariationAllele,
    ) -> None:
        """Test overlap_consequences property."""
        structural_variation_allele.transcript_variation_alleles = []
        assert (
            structural_variation_allele.overlap_consequences == []
        ), "Value for attribute overlap_consequences is not empty"

        tva1 = trancript_structural_variation_allele
        tva1.overlap_consequences = [OverlapConsequence(term="feature_ablation", accession="SO:0001879")]
        tva2 = evolve(trancript_structural_variation_allele)
        tva2.overlap_consequences = [OverlapConsequence(term="frameshift_truncation", accession="SO:0001910")]
        structural_variation_allele.transcript_variation_alleles += [tva1, tva2]

        assert structural_variation_allele.transcript_variation_alleles == [
            tva1,
            tva2,
        ], "Incorrect value for attribute transcript_variation_alleles"

        assert structural_variation_allele.overlap_consequences == [
            OverlapConsequence(term="feature_ablation", accession="SO:0001879"),
            OverlapConsequence(term="frameshift_truncation", accession="SO:0001910"),
        ], "Incorrect value for attribute overlap_consequences"
