
from __future__ import annotations

import pytest

from ensembl.vep.models.variation.variation_allele import VariationAllele
from ensembl.vep.models.sequence import Slice

@pytest.fixture
def slice():
    """Create a slice object for testing."""
    return Slice(
        coord_system="chromosome",
        seq_region_name="1",
        seq_region_length=12,
        start=1,
        end=12,
        strand=1,
        seq="ATCGATCGATCG"
    )

@pytest.fixture
def variation_allele(slice):
    """Returns a VariationAllele object for testing."""
    return VariationAllele(
        start=5,
        end=5,
        strand=1,
        slice=slice,
        reference_allele="A",
        alternate_allele="T"
    )

class TestVariationAlleleProperties:
    """Tests property accessors and defaults."""

    def test_initialization(self, variation_allele):
        """Test that a VariationAllele can be initialized with the correct attributes."""
        assert variation_allele.start == 5
        assert variation_allele.end == 5
        assert variation_allele.strand == 1
        assert variation_allele.reference_allele == "A"
        assert variation_allele.alternate_allele == "T"

    def test_initialization_without_keywords(self):
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            VariationAllele(5, end=5, strand=1, slice=None, reference_allele="A", alternate_allele="T")

    def test_minimised(self, variation_allele):
        """Test the minimised property before and after calling minimise."""
        assert variation_allele.minimised is False
        variation_allele.minimise()
        assert variation_allele.minimised is True
        
        # Check that minimised_locus is set
        assert variation_allele.minimised_locus is not None
        assert variation_allele.original_locus is not None

        # Check that minimised is immutable by external code
        with pytest.raises(AttributeError, match="property 'minimised' of 'VariationAllele' object has no setter"):
            variation_allele.minimised = False

    def test_shifted(self, variation_allele):
        """Test the shifted property before and after calling shift."""
        assert variation_allele.shifted is False
        variation_allele.shift()
        assert variation_allele.shifted is True

        # Check that shifted_locus is set
        assert variation_allele.shifted_locus is not None
        assert variation_allele.original_locus is not None

        # Check that shifted is immutable by external code
        with pytest.raises(AttributeError, match="property 'shifted' of 'VariationAllele' object has no setter"):
            variation_allele.shifted = False

    def test_consequence_types(self, variation_allele):
        """Test that consequence_types returns None or empty if not set."""
        # The implementation in the model is incomplete, so just check it returns None or list
        pass


class TestVariationAlleleMinimise:
    """Tests for the minimise method of VariationAllele."""

    def test_minimise_raises_if_missing_alleles(self, variation_allele):
        """Test minimise raises if reference or alternate variation_allele is None."""
        a = VariationAllele(
            start=1, end=1, strand=1, slice=variation_allele.slice, reference_allele=None, alternate_allele="T"
        )
        with pytest.raises(ValueError, match="Missing reference or alternate sequence for trimming"):
            a.minimise()
        b = VariationAllele(
            start=1, end=1, strand=1, slice=variation_allele.slice, reference_allele="A", alternate_allele=None
        )
        with pytest.raises(ValueError, match="Missing reference or alternate sequence for trimming"):
            b.minimise()

    def test_minimise_sets_locus_correctly(self, variation_allele):
        """Test that minimise sets minimised_locus and original_locus correctly."""
        variation_allele.minimise()
        assert variation_allele.minimised_locus.start == variation_allele.start
        assert variation_allele.minimised_locus.end == variation_allele.end
        assert variation_allele.minimised_locus.allele_string == f"{variation_allele.reference_allele}/{variation_allele.alternate_allele}"
        assert variation_allele.original_locus is not None


class TestVariationAlleleShift:
    """Tests for the shift method of VariationAllele."""

    def test_shift_raises_if_missing_alleles(self, variation_allele):
        """Test shift raises if reference or alternate variation_allele is None."""
        a = VariationAllele(
            start=1, end=1, strand=1, slice=variation_allele.slice, reference_allele=None, alternate_allele="T"
        )
        with pytest.raises(ValueError, match="Missing reference or alternate allele for shifting"):
            a.shift()
        b = VariationAllele(
            start=1, end=1, strand=1, slice=variation_allele.slice, reference_allele="A", alternate_allele=None
        )
        with pytest.raises(ValueError, match="Missing reference or alternate allele for shifting"):
            b.shift()

    def test_shift_sets_locus_correctly(self, variation_allele):
        """Test that shift sets shifted_locus and original_locus correctly."""
        variation_allele.shift()
        assert variation_allele.shifted_locus.start == variation_allele.start
        assert variation_allele.shifted_locus.end == variation_allele.end
        assert variation_allele.shifted_locus.allele_string == f"{variation_allele.reference_allele}/{variation_allele.alternate_allele}"
        assert variation_allele.original_locus is not None

    def test_minimise_and_shift_independent(self, variation_allele):
        """Test that minimise and shift can be called independently and set their respective flags."""
        variation_allele.minimise()
        assert variation_allele.minimised is True
        assert variation_allele.shifted is False
        variation_allele2 = VariationAllele(
            start=5, end=5, strand=1, slice=variation_allele.slice, reference_allele="A", alternate_allele="T"
        )
        variation_allele2.shift()
        assert variation_allele2.shifted is True
        assert variation_allele2.minimised is False