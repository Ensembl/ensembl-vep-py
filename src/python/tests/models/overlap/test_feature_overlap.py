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
"""Unit tests for FeatureOverlap class."""

import pytest

from ensembl.vep.common import Strand
from ensembl.vep.models import Feature
from ensembl.vep.models.overlap import FeatureOverlap, OverlapConsequence


@pytest.fixture
def feature() -> Feature:
    """Feature fixture for testing."""
    return Feature(start=1, end=100, strand=Strand.FORWARD)


@pytest.fixture
def overlapping_feature() -> Feature:
    """Overlapping Feature fixture for testing."""
    return Feature(start=50, end=150, strand=Strand.FORWARD)


@pytest.fixture
def feature_overlap(feature: Feature, overlapping_feature: Feature) -> FeatureOverlap:
    """FeatureOverlap fixture for testing."""
    return FeatureOverlap(feature=feature, overlapping_feature=overlapping_feature)


class TestFeatureOverlapInitialization:
    """Test FeatureOverlap initialization."""

    def test_initialization_without_keywords(self, feature: Feature, overlapping_feature: Feature) -> None:
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = FeatureOverlap(feature, overlapping_feature)  # type: ignore[misc]


class TestFeatureOverlapProperties:
    """Test FeatureOverlap properties."""

    def test_core_properties(self, feature_overlap: FeatureOverlap) -> None:
        """Test that feature and overlapping_feature properties are set correctly."""
        assert feature_overlap.feature == feature_overlap.feature
        assert feature_overlap.overlapping_feature == feature_overlap.overlapping_feature

    def test_most_severe_consequence(self, feature_overlap: FeatureOverlap) -> None:
        """Test that most_severe_consequence property can be accessed."""
        assert feature_overlap.most_severe_consequence is None

        feature_overlap.overlap_consequences += [
            OverlapConsequence(term="missense_variant", accession="SO:0001583", rank=2)
        ]
        assert feature_overlap.most_severe_consequence == feature_overlap.overlap_consequences[0]

        feature_overlap.overlap_consequences += [
            OverlapConsequence(term="synonymous_variant", accession="SO:0001819", rank=1)
        ]
        assert feature_overlap.most_severe_consequence == feature_overlap.overlap_consequences[0]

        # We cannot do feature_overlap.overlap_consequences.append()
        # as it would not trigger on_setattr, workaround?

    def test_consequence_terms(self, feature_overlap: FeatureOverlap) -> None:
        """Test that consequence_terms property can be accessed."""
        assert feature_overlap.consequence_terms == set()

        feature_overlap.overlap_consequences += [
            OverlapConsequence(term="missense_variant", accession="SO:0001583", rank=2)
        ]
        assert feature_overlap.consequence_terms == {"missense_variant"}

        feature_overlap.overlap_consequences += [
            OverlapConsequence(term="synonymous_variant", accession="SO:0001819", rank=1)
        ]
        assert feature_overlap.consequence_terms == {"missense_variant", "synonymous_variant"}
