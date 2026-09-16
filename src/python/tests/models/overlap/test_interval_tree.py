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

import logging
import pytest

from ensembl.vep.common import Strand
from ensembl.vep.models import Feature
from ensembl.vep.models.overlap.interval_tree import IntervalTree
from ensembl.vep.models.variation.variation import Variation

# logger = logging.getLogger(__name__)

_TEST_INTERVALS = [
    (
        9,
        10,
    ),
    (
        10,
        10,
    ),
    (
        10,
        11,
    ),
    (
        10,
        12,
    ),
    (
        11,
        11,
    ),
    (
        11,
        12,
    ),
    (
        12,
        14,
    ),
    (
        13,
        13,
    ),
    (
        14,
        15,
    ),
    (
        14,
        16,
    ),
    (
        15,
        15,
    ),
    (
        15,
        16,
    ),
    (
        16,
        16,
    ),
    (
        16,
        17,
    ),
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


class TestIntervalTreeProperties:
    """Tests property accessors and default values."""

    def test_initialization(self, variations_for_overlaps: list[Variation]) -> None:
        """Test that a IntervalTree can be initialized with the correct attributes."""
        interval_tree: IntervalTree = IntervalTree()
        assert interval_tree.features == []

        interval_tree = IntervalTree(features=variations_for_overlaps)
        assert interval_tree.features == variations_for_overlaps
        assert interval_tree._intervals.size == len(  # pylint: disable=protected-access
            variations_for_overlaps
        )

    def test_initialization_without_keywords(self) -> None:
        """Test IntervalTree initialization without keyword raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = IntervalTree([100])  # type: ignore[misc]


class TestIntervalTreeMethods:
    """Tests overlap methods."""

    def test_add_features(self, variations_for_overlaps: list[Variation]) -> None:
        """Test the ``add`` method."""
        interval_tree: IntervalTree = IntervalTree()
        interval_tree.add_features(variations_for_overlaps)
        # added features are added to the .features list and to the AIList in _intervals
        assert interval_tree.features == variations_for_overlaps
        assert interval_tree._intervals.size == len(  # pylint: disable=protected-access
            variations_for_overlaps
        )

        # added features are added to the end of the .features list so that they can be retrieved by index
        new_var = Variation(start=12, end=12, strand=Strand.FORWARD)
        interval_tree.add_features([new_var])
        variations_for_overlaps.append(new_var)
        assert interval_tree.features == variations_for_overlaps
        assert interval_tree._intervals.size == len(  # pylint: disable=protected-access
            variations_for_overlaps
        )

        # add a feature that is already in the list
        # it gets added with a different index
        new_var = Variation(start=12, end=12, strand=Strand.FORWARD)
        interval_tree.add_features([new_var])
        variations_for_overlaps.append(new_var)
        assert interval_tree.features == variations_for_overlaps
        assert interval_tree._intervals.size == len(  # pylint: disable=protected-access
            variations_for_overlaps
        )

    def test_get_overlapping_features(self, variations_for_overlaps: list[Variation]) -> None:
        """Test the ``get_overlapping_features`` method."""
        interval_tree: IntervalTree = IntervalTree()
        interval_tree.add_features(variations_for_overlaps)

        overlapping_vars = interval_tree.get_overlapping_features(start=10, end=10)
        expected_overlaps = [variations_for_overlaps[idx] for idx in range(4)]
        assert len(overlapping_vars) == len(expected_overlaps)
        assert overlapping_vars == expected_overlaps

        overlapping_vars = interval_tree.get_overlapping_features(start=11, end=15)
        expected_overlaps = [variations_for_overlaps[idx] for idx in range(2, 12)]
        assert len(overlapping_vars) == len(expected_overlaps)
        assert overlapping_vars == expected_overlaps

    def test_get_overlapping_features_returns_empty_list_without_matches(
        self, variations_for_overlaps: list[Variation]
    ) -> None:
        """Test the ``get_overlapping_features`` method when no intervals overlap."""
        interval_tree: IntervalTree = IntervalTree()
        interval_tree.add_features(variations_for_overlaps)

        assert interval_tree.get_overlapping_features(start=1, end=8) == []
        assert interval_tree.get_overlapping_features(start=18, end=20) == []

    def test_get_overlapping_features_returns_empty_list_without_features(self) -> None:
        """Test the ``get_overlapping_features`` method before any intervals are added."""
        interval_tree: IntervalTree = IntervalTree()

        assert interval_tree.get_overlapping_features(start=1, end=1) == []

    def test_get_overlapping_features_by_feature(self, variations_for_overlaps: list[Variation]) -> None:
        """Test the ``get_overlapping_features_by_feature`` method."""
        interval_tree: IntervalTree = IntervalTree()
        interval_tree.add_features(variations_for_overlaps)
        query_feature = Feature(start=10, end=10, strand=Strand.FORWARD)

        overlapping_vars = interval_tree.get_overlapping_features_by_feature(query_feature)

        assert overlapping_vars == [variations_for_overlaps[idx] for idx in range(4)]

    def test_get_nearest_features(self, variations_for_nearest_interval: list[Variation]) -> None:
        """Test the ``get_nearest_features`` method."""
        interval_tree: IntervalTree = IntervalTree()
        interval_tree.add_features(variations_for_nearest_interval)

        nearest = interval_tree.get_nearest_features(start=13, end=13, k=1)
        assert nearest == [variations_for_nearest_interval[2]]

        nearest = interval_tree.get_nearest_features(start=12, end=12, k=2)
        assert nearest == variations_for_nearest_interval[:2]

    def test_get_nearest_features_by_feature(self, variations_for_nearest_interval: list[Variation]) -> None:
        """Test the ``get_nearest_features_by_feature`` method."""
        interval_tree: IntervalTree = IntervalTree()
        interval_tree.add_features(variations_for_nearest_interval)
        query_feature = Feature(start=13, end=13, strand=Strand.FORWARD)

        nearest = interval_tree.get_nearest_features_by_feature(query_feature, k=1)

        assert nearest == [variations_for_nearest_interval[2]]

    @pytest.mark.parametrize(
        ("start", "end", "message"),
        [
            (0, 1, "start must be 1-based"),
            (2, 1, "start must be less than or equal to end"),
        ],
    )
    def test_add_features_with_invalid_interval_raises(
        self, start: int, end: int, message: str, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Test the ``add_features`` method raises an error for invalid intervals."""
        interval_tree: IntervalTree = IntervalTree()
        variation = Variation(start=start, end=end, strand=Strand.FORWARD)

        with pytest.raises(ValueError, match=message):
            interval_tree._check_interval(  # pylint: disable=protected-access
                start=variation.start, end=variation.end
            )

        with caplog.at_level(logging.WARNING):
            interval_tree.add_features([variation])

            assert "Feature has invalid coordinates" in caplog.text, "Incorrect warning on invalid interval"

        # the feature is still added to the ``features`` attribute, but not to the ``_intervals``
        # easier than altering the indexes
        assert interval_tree.features == [variation]
        assert interval_tree._intervals.size == 0  # pylint: disable=protected-access

    @pytest.mark.parametrize(
        ("method_name", "start", "end", "message"),
        [
            ("get_overlapping_features", 0, 1, "start must be 1-based"),
            ("get_overlapping_features", 2, 1, "start must be less than or equal to end"),
            ("get_nearest_features", 0, 1, "start must be 1-based"),
            ("get_nearest_features", 2, 1, "start must be less than or equal to end"),
        ],
    )
    def test_query_methods_with_invalid_interval_raise(
        self, method_name: str, start: int, end: int, message: str
    ) -> None:
        """Test query methods raise an error for invalid intervals."""
        interval_tree: IntervalTree = IntervalTree()
        method = getattr(interval_tree, method_name)

        with pytest.raises(ValueError, match=message):
            method(start=start, end=end)
