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
from ensembl.vep.models.overlap.interval_overlap import IntervalOverlap
from ensembl.vep.models.variation.variation import Variation

_TEST_INTERVALS = [
    (9, 10,),
    (10, 10,),
    (10, 11,),
    (10, 12,),
    (11, 11,),
    (11, 12,),
    (12, 14,),
    (13, 13,),
    (14, 15,),
    (14, 16,),
    (15, 15,),
    (15, 16,),
    (16, 16,),
    (16, 17,),
]

@pytest.fixture
def variations_for_overlaps() -> list[Variation]:
    """Create a list of Variation objects for testing."""
    var_list = []
    for interval in _TEST_INTERVALS:
        var_list.append(Variation(start=interval[0], end=interval[1], strand=Strand.FORWARD))

    return var_list


@pytest.fixture
def variations_for_nearest_interval() -> list[Variation]:
    """Create a list of Variation objects for testing the ``get_nearest_features`` method.
    
    Three intervals: (9,10), (10,10), (15,15)
    For the interval (13,13) there is a unique nearest interval to the right [(15,15)].
    For the interval (12,12) there are two intervals to the left whose ends are the same distance 
    away [(9,10), (10,10)].
    """
    var_list = []
    for interval in [_TEST_INTERVALS[idx] for idx in (0, 1, 10)]:
        var_list.append(Variation(start=interval[0], end=interval[1], strand=Strand.FORWARD))

    return var_list


class TestIntervalOverlapProperties:
    """Tests property accessors and default values."""

    def test_initialization(self, variations_for_overlaps: list[Variation]) -> None:
        """Test that a IntervalOverlap can be initialized with the correct attributes."""
        interval_overlap = IntervalOverlap()
        assert interval_overlap.features == []

        interval_overlap = IntervalOverlap(features=variations_for_overlaps)
        assert interval_overlap.features == variations_for_overlaps

    def test_initialization_without_keywords(self) -> None:
        """Test IntervalOverlap initialization without keyword raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = IntervalOverlap(
                100  # type: ignore[misc]
            )


class TestIntervalOverlapMethods:
    """Tests overlap methods."""

    def test_add_features(self, variations_for_overlaps: list[Variation]) -> None:
        """Test the ``add`` method."""
        interval_overlap = IntervalOverlap()
        interval_overlap.add_features(variations_for_overlaps)
        # added features are added to the .features list and to the AIList in _intervals
        assert interval_overlap.features == variations_for_overlaps
        assert interval_overlap._intervals.size == len(variations_for_overlaps)  # pylint: disable=protected-access
        
        # added features are added to the end of the .features list so that they can be retrieved by index
        new_var = Variation(start=12, end=12, strand=Strand.FORWARD)
        interval_overlap.add_features([new_var])
        variations_for_overlaps.append(new_var)
        assert interval_overlap.features == variations_for_overlaps
        assert interval_overlap._intervals.size == len(variations_for_overlaps)  # pylint: disable=protected-access

        # add a feature that is already in the list
        new_var = Variation(start=12, end=12, strand=Strand.FORWARD)
        interval_overlap.add_features([new_var])
        variations_for_overlaps.append(new_var)
        assert interval_overlap.features == variations_for_overlaps
        assert interval_overlap._intervals.size == len(variations_for_overlaps)  # pylint: disable=protected-access

    def test_get_overlapping_features(self, variations_for_overlaps: list[Variation]) -> None:
        """Test the ``get_overlapping_features`` method."""
        interval_overlap = IntervalOverlap()
        interval_overlap.add_features(variations_for_overlaps)

        overlapping_vars = interval_overlap.get_overlapping_features(start=10, end=10)
        expected_overlaps = [variations_for_overlaps[idx] for idx in range(4)]
        assert len(overlapping_vars) == len(expected_overlaps)
        assert overlapping_vars == expected_overlaps

        overlapping_vars = interval_overlap.get_overlapping_features(start=11, end=15)
        expected_overlaps = [variations_for_overlaps[idx] for idx in range(2,12)]
        assert len(overlapping_vars) == len(expected_overlaps)
        assert overlapping_vars == expected_overlaps

    def test_get_nearest_features(self, variations_for_nearest_interval: list[Variation]) -> None:
        """Test the ``get_nearest_features`` method."""
        interval_overlap = IntervalOverlap()
        interval_overlap.add_features(variations_for_nearest_interval)

        nearest = interval_overlap.get_nearest_features(start=13, end=13, k=1)
        assert nearest == [variations_for_nearest_interval[2]]

        nearest = interval_overlap.get_nearest_features(start=12, end=12, k=2)
        assert nearest == variations_for_nearest_interval[:2]
