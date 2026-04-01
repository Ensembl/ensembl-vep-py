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
"""Unit tests for OverlapConsequence class."""

import pytest

from ensembl.vep.models.overlap import OverlapConsequence
from ensembl.vep.common import ConsequenceImpact
from ensembl.vep.models import Feature

@pytest.fixture
def overlap_consequence():
    return OverlapConsequence(term="missense_variant", accession="SO:0001583")

class TestOverlapConsequenceInitialization:
    """Test OverlapConsequence initialization."""

    def test_initialization_without_keywords(self, overlap_consequence):
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = OverlapConsequence(overlap_consequence)


class TestOverlapConsequenceProperties:
    """Test OverlapConsequence properties."""

    def test_core_properties(self, overlap_consequence):
        """Test that term and accession properties are set correctly."""
        assert overlap_consequence.term == "missense_variant"
        assert overlap_consequence.accession == "SO:0001583"

    def test_default_properties(self, overlap_consequence):
        """Test that default properties are set correctly."""
        assert not overlap_consequence.is_default
        assert overlap_consequence.rank == 0
        assert overlap_consequence.tier == 0
        assert overlap_consequence.impact == ConsequenceImpact.LOW
        assert overlap_consequence.feature_type is type(Feature)
        assert overlap_consequence.overlapping_feature_type is type(Feature)