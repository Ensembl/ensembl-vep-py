import pytest

from ensembl.vep.common import Locus

class TestLocus:
    """Tests for the Locus frozen attrs class."""

    def test_initialization(self):
        """Test that a Locus can be initialized with the correct attributes."""
        locus = Locus(start=100, end=200, allele_string="A/T")
        assert locus.start == 100
        assert locus.end == 200
        assert locus.allele_string == "A/T"

    def test_immutability(self):
        """Test that the attributes of Locus cannot be modified after initialization."""
        locus = Locus(start=100, end=200, allele_string="A/T")
        with pytest.raises(AttributeError):
            locus.start = 200

    def test_equality(self):
        """Test that two Locus instances with the same attributes are considered equal."""
        locus1 = Locus(start=100, end=200, allele_string="A/T")
        locus2 = Locus(start=100, end=200, allele_string="A/T")
        assert locus1 == locus2

    def test_inequality(self):
        """Test that two Locus instances with different attributes are considered not equal."""
        locus1 = Locus(start=100, end=200, allele_string="A/T")
        locus2 = Locus(start=100, end=200, allele_string="A/G")
        assert locus1 != locus2
