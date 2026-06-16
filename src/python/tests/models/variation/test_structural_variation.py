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
"""Unit tests for the StructuralVariation class."""

from __future__ import annotations

import pytest

from ensembl.vep.common import ClassSOTerm, CoordSystemLevel, Strand
from ensembl.vep.models.variation import StructuralVariation, StructuralVariationAllele, OverlapConsequence
from ensembl.vep.models.sequence import CoordSystem, Slice


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
def structural_variation(genomic_slice: Slice) -> StructuralVariation:
    """Returns a deletion StructuralVariation object for testing."""
    return StructuralVariation(
        start=5,
        end=150,
        strand=Strand.FORWARD,
        slice=genomic_slice,
        allele_string="A/DEL",
    )


@pytest.fixture
def structural_variation_no_alleles(genomic_slice: Slice) -> StructuralVariation:
    """Returns a StructuralVariation object with no allele_string for testing."""
    return StructuralVariation(
        start=5,
        end=150,
        strand=Strand.FORWARD,
        slice=genomic_slice,
    )


@pytest.fixture
def insertion_structural_variation(genomic_slice: Slice) -> StructuralVariation:
    """Returns an insertion StructuralVariation object for testing."""
    return StructuralVariation(
        start=5,
        end=5,
        strand=Strand.FORWARD,
        slice=genomic_slice,
        allele_string="-/INS",
    )


@pytest.fixture
def tandem_repeat_structural_variation(genomic_slice: Slice) -> StructuralVariation:
    """Returns a tandem_repeat StructuralVariation object for testing."""
    return StructuralVariation(
        start=5,
        end=5,
        strand=Strand.FORWARD,
        slice=genomic_slice,
        allele_string="(CAG)9/(CAG)12",
    )


@pytest.fixture
def multiallelic_structural_variation(genomic_slice: Slice) -> StructuralVariation:
    """Returns a tandem_repeat StructuralVariation object for testing."""
    return StructuralVariation(
        start=5,
        end=5,
        strand=Strand.FORWARD,
        slice=genomic_slice,
        allele_string="(CAG)9/(CAG)12/(CAG)15",
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
def overlap_consequence() -> OverlapConsequence:
    """Returns an OverlapConsequence object for testing."""
    return OverlapConsequence(term="transcript_ablation", accession="SO:0001893")


class TestStructuralVariationProperties:
    """Tests property accessors and default values."""

    def test_initialization(self, structural_variation: StructuralVariation) -> None:
        """Test that a StructuralVariation can be initialized with the correct attributes."""
        assert structural_variation.start == 5, "start attribute incorrectly initialized"
        assert structural_variation.end == 150, "end attribute incorrectly initialized"
        assert structural_variation.strand == Strand.FORWARD, "strand attribute incorrectly initialized"
        assert (
            structural_variation.allele_string == "A/DEL"
        ), "allele_string attribute incorrectly initialized"

    def test_initialization_without_keywords(self) -> None:
        """Test Variation initialization without keyword raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = StructuralVariation(
                100,
                end=250,
                strand=Strand.FORWARD,
                slice=None,
                allele_string="A/DEL",  # type: ignore[misc]
            )

    def test_defaults(self, structural_variation: StructuralVariation) -> None:
        """Test default values for optional attributes."""
        assert structural_variation.stable_id is None, "stable_id attribute incorrect default"
        assert structural_variation.is_somatic is False, "is_somatic attribute incorrect default"
        assert structural_variation.copy_number is None, "copy_number attribute incorrect default"
        assert structural_variation.inner_start is None, "inner_start attribute incorrect default"
        assert structural_variation.outer_start is None, "outer_start attribute incorrect default"
        assert structural_variation.inner_end is None, "inner_end attribute incorrect default"
        assert structural_variation.outer_end is None, "outer_end attribute incorrect default"
        assert structural_variation.variation_alleles == [], "variation_alleles attribute incorrect default"

    def test_sva_assigment(
        self,
        structural_variation: StructuralVariation,
        structural_variation_allele: StructuralVariationAllele,
    ) -> None:
        """Test assigning to variation_alleles."""
        structural_variation.variation_alleles = [structural_variation_allele]
        assert structural_variation.variation_alleles == [
            structural_variation_allele
        ], "variation_alleles attribute incorrectly assigned"

    def test_most_severe_consequence_immutable(
        self, structural_variation: StructuralVariation, overlap_consequence: OverlapConsequence
    ) -> None:
        """Test that most_severe_consequence cannot be directly set."""
        with pytest.raises(
            AttributeError,
            match="property 'most_severe_consequence' of 'StructuralVariation' object has no setter",
        ):
            structural_variation.most_severe_consequence = overlap_consequence  # type: ignore[misc]

    def test_class_so_term_immutable(self, structural_variation: StructuralVariation) -> None:
        """Test that class_so_term cannot be directly set."""
        with pytest.raises(
            AttributeError, match="property 'class_so_term' of 'StructuralVariation' object has no setter"
        ):
            structural_variation.class_so_term = ClassSOTerm.SEQUENCE_ALTERATION  # type: ignore[misc]


class TestClassSOTerm:
    """Tests for the class_so_term cached property."""

    def test_correct_term_deletion(self, structural_variation: StructuralVariation) -> None:
        """Test StructuralVariation object for deletion produces correct ClassSOTerm"""
        assert (
            structural_variation.class_so_term == ClassSOTerm.DELETION
        ), "class_so_term for deletion is incorrect"

    def test_correct_term_insertion(self, insertion_structural_variation: StructuralVariation) -> None:
        """Test StructuralVariation object for deletion produces correct ClassSOTerm"""
        assert (
            insertion_structural_variation.class_so_term == ClassSOTerm.INSERTION
        ), "class_so_term for insertion is incorrect"

    def test_correct_term_tandem_repeat(
        self, tandem_repeat_structural_variation: StructuralVariation
    ) -> None:
        """Test StructuralVariation object for deletion produces correct ClassSOTerm"""
        assert (
            tandem_repeat_structural_variation.class_so_term == ClassSOTerm.TANDEM_REPEAT
        ), "class_so_term for tandem repeat is incorrect"


class TestMostSevereConsequence:
    """Tests for the most_severe_consequence cached property."""

    def test_no_csqs(self, structural_variation_no_alleles: StructuralVariation) -> None:
        """Test that most severe consequence is None with no allele string"""
        assert (
            structural_variation_no_alleles.most_severe_consequence is None
        ), "most_severe_consequence for empty allele string is incorrect"
