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
