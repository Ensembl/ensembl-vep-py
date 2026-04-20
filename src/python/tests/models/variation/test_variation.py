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
"""Unit tests for the Variation class."""

from __future__ import annotations

import pytest

from ensembl.vep.common import ClassSOTerm, Locus, Strand
from ensembl.vep.models.variation import Variation, OverlapConsequence
from ensembl.vep.models.sequence import Slice


@pytest.fixture
def genomic_slice() -> Slice:
    """Create a slice object for testing."""
    return Slice(
        coord_system="chromosome",
        seq_region_name="1",
        seq_region_length=12,
        start=1,
        end=12,
        strand=Strand.FORWARD,
        seq="ATCGATCGATCG",
    )


@pytest.fixture
def variation(genomic_slice: Slice) -> Variation:
    """Returns a bi-allelic Variation object for testing."""
    return Variation(
        start=5,
        end=5,
        strand=Strand.FORWARD,
        slice=genomic_slice,
        allele_string="A/T",
    )


@pytest.fixture
def multiallelic_variation(genomic_slice: Slice) -> Variation:
    """Returns a multiallelic Variation object for testing."""
    return Variation(
        start=5,
        end=5,
        strand=Strand.FORWARD,
        slice=genomic_slice,
        allele_string="A/T/G/C",
    )


@pytest.fixture
def insertion_variation(genomic_slice: Slice) -> Variation:
    """Returns an insertion Variation object for testing."""
    return Variation(
        start=5,
        end=5,
        strand=Strand.FORWARD,
        slice=genomic_slice,
        allele_string="A/ATT",
    )


@pytest.fixture
def deletion_variation(genomic_slice: Slice) -> Variation:
    """Returns a deletion Variation object for testing."""
    return Variation(
        start=5,
        end=7,
        strand=Strand.FORWARD,
        slice=genomic_slice,
        allele_string="ATC/A",
    )

@pytest.fixture
def overlap_consequence() -> OverlapConsequence:
    """Returns an OverlapConsequence object for testing."""
    return OverlapConsequence(term="missense_variant", accession="SO:0001583")


class TestVariationProperties:
    """Tests property accessors and default values."""

    def test_initialization(self, variation: Variation) -> None:
        """Test that a Variation can be initialized with the correct attributes."""
        assert variation.start == 5
        assert variation.end == 5
        assert variation.strand == 1
        assert variation.allele_string == "A/T"

    def test_initialization_without_keywords(self) -> None:
        """Test Variation initialization without keyword raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = Variation(start=100, end=100, strand=Strand.FORWARD, slice=None, allele_string="A/T")

    def test_defaults(self, variation: Variation) -> None:
        """Test default values for optional attributes."""
        assert variation.stable_id is None
        assert variation.is_somatic is False

    def test_variation_alleles_immutable(self, variation: Variation) -> None:
        """Test that variation_alleles cannot be directly set."""
        with pytest.raises(
            AttributeError, match="property 'variation_alleles' of 'Variation' object has no setter"
        ):
            variation.variation_alleles = []  # type: ignore[misc]

    def test_most_severe_consequence_immutable(
        self, variation: Variation, overlap_consequence: OverlapConsequence
    ) -> None:
        """Test that most_severe_consequence cannot be directly set."""
        with pytest.raises(
            AttributeError, match="property 'most_severe_consequence' of 'Variation' object has no setter"
        ):
            variation.most_severe_consequence = overlap_consequence  # type: ignore[misc]

    def test_class_so_term_immutable(self, variation: Variation) -> None:
        """Test that class_so_term cannot be directly set."""
        with pytest.raises(
            AttributeError, match="property 'class_so_term' of 'Variation' object has no setter"
        ):
            variation.class_so_term = ClassSOTerm.SEQUENCE_ALTERATION  # type: ignore[misc]


class TestVariationAlleles:
    """Tests for the variation_alleles cached property."""

    def test_single_alt_allele(self, variation: Variation) -> None:
        """Test that a bi-allelic variation returns the correct VariationAllele."""
        alleles = variation.variation_alleles

        assert len(alleles) == 1
        assert alleles[0].reference_allele == "A"
        assert alleles[0].alternate_allele == "T"

    def test_multiple_alt_alleles(self, multiallelic_variation: Variation) -> None:
        """Test that a multi-allelic variation returns the correct VariationAlleles."""
        alleles = multiallelic_variation.variation_alleles

        assert len(alleles) == 3
        assert alleles[0].alternate_allele == "T"
        assert alleles[1].alternate_allele == "G"
        assert alleles[2].alternate_allele == "C"
        # All share the same reference
        assert all(a.reference_allele == "A" for a in alleles)

    def test_insertion_alleles(self, insertion_variation: Variation) -> None:
        """Test that an insertion variation returns the correct VariationAllele."""
        alleles = insertion_variation.variation_alleles

        assert len(alleles) == 1
        assert alleles[0].reference_allele == "A"
        assert alleles[0].alternate_allele == "ATT"

    def test_deletion_alleles(self, deletion_variation: Variation) -> None:
        """Test that a deletion variation returns the correct VariationAllele."""
        alleles = deletion_variation.variation_alleles

        assert len(alleles) == 1
        assert alleles[0].reference_allele == "ATC"
        assert alleles[0].alternate_allele == "A"

    def test_alleles_inherit_coordinates(self, variation: Variation) -> None:
        """Test that variation alleles inherit coordinates from the variation."""
        alleles = variation.variation_alleles

        assert alleles[0].start == 5
        assert alleles[0].end == 5
        assert alleles[0].strand == 1
        assert alleles[0].slice == variation.slice

    def test_none_allele_string_returns_empty(self, genomic_slice: Slice) -> None:
        """Test that an undefined allele string (None) returns an empty list of alleles."""
        var = Variation(
            start=100,
            end=100,
            strand=Strand.FORWARD,
            slice=genomic_slice,
            allele_string=None,
        )
        alleles = var.variation_alleles
        assert alleles == []


class TestVariationMinimisation:
    """Tests for the minimise_alleles of Variation."""

    def test_minimise_called_when_flag_set(self) -> None:
        """Test that the minimise_alleles method is called when the minimise_alleles flag is set."""
        var = Variation(
            start=100,
            end=102,
            strand=Strand.FORWARD,
            slice=None,
            allele_string="AAT/AT",
        )
        var.minimise_alleles()
        for variantion_allele in var.variation_alleles:
            assert variantion_allele.minimised is True
            assert variantion_allele.reference_allele == "A"
            assert variantion_allele.alternate_allele == ""
            assert variantion_allele.start == 101
            assert variantion_allele.end == 101
            assert variantion_allele.original_locus == Locus(start=100, end=102, allele_string="AAT/AT")
            assert variantion_allele.minimised_locus == Locus(start=101, end=101, allele_string="A/")


class TestVariationShifting:
    """Tests for the shift_alleles of Variation."""

    def test_shift_alleles(self, genomic_slice: Slice) -> None:
        """Test that the shift_alleles method shifts the alleles correctly."""
        var = Variation(
            start=2,
            end=5,
            strand=Strand.FORWARD,
            slice=genomic_slice,
            allele_string="TCGA/-",
        )
        var.shift_alleles()
        for variantion_allele in var.variation_alleles:
            assert variantion_allele.shifted is True
            assert variantion_allele.reference_allele == "TCGA"
            assert variantion_allele.alternate_allele == "-"
            assert variantion_allele.start == 6
            assert variantion_allele.end == 9
            assert variantion_allele.original_locus == Locus(start=2, end=5, allele_string="TCGA/-")
            assert variantion_allele.shifted_locus == Locus(start=6, end=9, allele_string="TCGA/-")


class TestMostSevereConsequence:
    """Tests for the most_severe_consequence cached property."""

    pass


class TestClassSOTerm:
    """Tests for the class_so_term cached property."""

    @pytest.mark.parametrize(
        "allele_string, expected",
        [
            ("A/T", ClassSOTerm.SNV),
            ("A/G", ClassSOTerm.SNV),
            ("A/T/G", ClassSOTerm.SNV),
            ("A/T/G/C", ClassSOTerm.SNV),
        ],
    )
    def test_snv(self, genomic_slice: Slice, allele_string: str, expected: ClassSOTerm) -> None:
        """Test that a SNV variation returns the correct ClassSOTerm."""
        var = Variation(start=100, end=100, strand=Strand.FORWARD,
                        slice=genomic_slice, allele_string=allele_string)
        assert var.class_so_term == expected

    @pytest.mark.parametrize(
        "allele_string, expected",
        [
            ("(CAG)8/(CAG)9", ClassSOTerm.TANDEM_REPEAT),
            ("(AT)5/(AT)7", ClassSOTerm.TANDEM_REPEAT),
        ],
    )
    def test_tandem_repeat(self, genomic_slice: Slice, allele_string: str, expected: ClassSOTerm) -> None:
        """Test that a tandem repeat variation returns the correct ClassSOTerm."""
        var = Variation(start=100, end=100, strand=Strand.FORWARD,
                        slice=genomic_slice, allele_string=allele_string)
        assert var.class_so_term == expected

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
    def test_insertion(self, genomic_slice: Slice, allele_string: str, expected: ClassSOTerm) -> None:
        """Test that an insertion variation returns the correct ClassSOTerm."""
        var = Variation(start=100, end=100, strand=Strand.FORWARD,
                        slice=genomic_slice, allele_string=allele_string)
        assert var.class_so_term == expected

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
    def test_deletion(self, genomic_slice: Slice, allele_string: str, expected: ClassSOTerm) -> None:
        """Test that a deletion variation returns the correct ClassSOTerm."""
        var = Variation(start=100, end=102, strand=Strand.FORWARD,
                        slice=genomic_slice, allele_string=allele_string)
        assert var.class_so_term == expected

    @pytest.mark.parametrize(
        "allele_string, expected",
        [
            ("AA/TT", ClassSOTerm.SUBSTITUTION),
            ("ACT/TCG", ClassSOTerm.SUBSTITUTION),
        ],
    )
    def test_substitution(self, genomic_slice: Slice, allele_string: str, expected: ClassSOTerm) -> None:
        """Test that a substitution variation returns the correct ClassSOTerm."""
        var = Variation(start=100, end=101, strand=Strand.FORWARD,
                        slice=genomic_slice, allele_string=allele_string)
        assert var.class_so_term == expected

    @pytest.mark.parametrize(
        "allele_string, expected",
        [
            ("AA/TTTT", ClassSOTerm.INDEL),
            ("AAAA/TT", ClassSOTerm.INDEL),
        ],
    )
    def test_indel(self, genomic_slice: Slice, allele_string: str, expected: ClassSOTerm) -> None:
        """Test that an indel variation returns the correct ClassSOTerm."""
        var = Variation(start=100, end=103, strand=Strand.FORWARD,
                        slice=genomic_slice, allele_string=allele_string)
        assert var.class_so_term == expected

    def test_empty_allele_string_defaults_to_sequence_alteration(self, genomic_slice: Slice) -> None:
        """Test that an empty allele string returns the correct (default) ClassSOTerm."""
        var = Variation(start=100, end=100, strand=Strand.FORWARD, slice=genomic_slice, allele_string="")
        assert var.class_so_term == ClassSOTerm.SEQUENCE_ALTERATION

    def test_none_allele_string_defaults_to_sequence_alteration(self, genomic_slice: Slice) -> None:
        """Test that an undefined (None) allele string returns the correct (default) ClassSOTerm."""
        var = Variation(start=100, end=100, strand=Strand.FORWARD, slice=genomic_slice, allele_string=None)
        assert var.class_so_term == ClassSOTerm.SEQUENCE_ALTERATION
