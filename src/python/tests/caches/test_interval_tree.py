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
"""Tests for interval queries, interval trees and the cache that shares them."""

# pylint: disable=redefined-outer-name,protected-access

from __future__ import annotations

import logging
import pickle
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor

import pytest
from attrs import define, field

from ensembl.vep.caches.interval_tree import (
    FeatureIntervalTreeCache,
    Interval,
    IntervalTree,
    IntervalTreeCache,
)
from ensembl.vep.common import SourceType, Strand
from ensembl.vep.models import Feature
from ensembl.vep.models.variation import Variation
from ensembl.vep.sources.feature_source import FeatureSource

_TEST_INTERVALS = [(9, 10), (10, 10), (14, 13), (15, 15)]


@define(kw_only=True, auto_attribs=True)
class _Span:
    """Minimal interval object, held where a real feature could not be pickled."""

    start: int
    end: int


@define(kw_only=True, auto_attribs=True)
class _DummyFeatureSource(FeatureSource[Feature, Feature]):
    """Minimal feature source for cache identity tests, whose record is its feature."""

    features: list[Feature] = field(factory=list)

    def get_all_records(self) -> list[Feature]:
        """Return all stored features."""

        return self.features

    def iter_all_records(self) -> Iterator[Feature]:
        """Return an iterator over stored features."""

        return iter(self.features)

    def get_records_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
        **kwargs: object,
    ) -> list[Feature]:
        """Return all stored features."""

        return self.features

    def iter_records_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
        **kwargs: object,
    ) -> Iterator[Feature]:
        """Return an iterator over stored features."""

        return iter(self.features)

    def records_to_features(self, records: list[Feature]) -> list[Feature]:
        """Return the records, which are already features."""

        return list(records)

    def close(self) -> None:
        """Clear stored features."""

        self.features = []


@pytest.fixture
def variations_for_overlaps() -> list[Variation]:
    """Create variations covering touching, nested and reversed intervals."""
    return [Variation(start=start, end=end, strand=Strand.FORWARD) for start, end in _TEST_INTERVALS]


@pytest.fixture
def variations_for_nearest_interval() -> list[Variation]:
    """Create variations at (9,10), (10,10) and (15,15).

    (13,13) has one nearest interval to its right; (12,12) has two equally near to its left.
    """
    return [
        Variation(start=_TEST_INTERVALS[index][0], end=_TEST_INTERVALS[index][1], strand=Strand.FORWARD)
        for index in (0, 1, 3)
    ]


@pytest.fixture
def interval_tree(variations_for_overlaps: list[Variation]) -> IntervalTree[Variation]:
    """Create an interval tree holding the overlap fixtures."""
    return IntervalTree(features=variations_for_overlaps)


@pytest.fixture
def nearest_interval_tree(variations_for_nearest_interval: list[Variation]) -> IntervalTree[Variation]:
    """Create an interval tree holding the proximity fixtures."""
    return IntervalTree(features=variations_for_nearest_interval)


class TestInterval:
    """Interval coordinate validation."""

    def test_identifies_insertion_point(self) -> None:
        """Identify an inverted one-base span as an insertion."""
        interval = Interval(start=10, end=9)

        assert interval.is_insertion
        assert interval.lookup_start == 9
        assert interval.lookup_end == 10

    @pytest.mark.parametrize(
        ("start", "end"),
        [
            pytest.param(0, 0, id="zero_start"),
            pytest.param(10, 8, id="more_than_one_inverted_base"),
        ],
    )
    def test_rejects_invalid_coordinates(self, start: int, end: int) -> None:
        """Reject coordinates outside the normal or insertion representation."""
        with pytest.raises(ValueError):
            Interval(start=start, end=end)


class TestIntervalTreeProperties:
    """Property accessors and default values."""

    def test_initialization(
        self, interval_tree: IntervalTree[Variation], variations_for_overlaps: list[Variation]
    ) -> None:
        """Index every feature supplied at construction."""
        empty_tree: IntervalTree[Variation] = IntervalTree()

        assert empty_tree.features == []
        assert interval_tree.features == variations_for_overlaps
        assert interval_tree._intervals.size == len(variations_for_overlaps)

    def test_initialization_without_keywords(self) -> None:
        """Refuse positional construction."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = IntervalTree([100])  # type: ignore[misc]


class TestIntervalTreeAddFeatures:
    """Adding features after construction."""

    def test_add_features_appends_and_indexes(self, variations_for_overlaps: list[Variation]) -> None:
        """Index added features under positions that follow the ones already held."""
        interval_tree: IntervalTree[Variation] = IntervalTree()
        interval_tree.add_features(variations_for_overlaps)

        assert interval_tree.features == variations_for_overlaps
        assert interval_tree._intervals.size == len(variations_for_overlaps)

        added = Variation(start=12, end=12, strand=Strand.FORWARD)
        interval_tree.add_features([added])

        assert interval_tree.features == [*variations_for_overlaps, added]
        assert interval_tree._intervals.size == len(variations_for_overlaps) + 1

    def test_add_features_indexes_an_equal_feature_separately(self) -> None:
        """Give an equal feature its own position rather than merging it."""
        first = Variation(start=12, end=12, strand=Strand.FORWARD)
        second = Variation(start=12, end=12, strand=Strand.FORWARD)
        interval_tree: IntervalTree[Variation] = IntervalTree(features=[first])
        interval_tree.add_features([second])

        assert interval_tree.get_overlapping_features(start=12, end=12) == [first, second]

    def test_add_features_warns_on_invalid_coordinates(self, caplog: pytest.LogCaptureFixture) -> None:
        """Warn about a feature that cannot be indexed, and hold it anyway."""
        interval_tree: IntervalTree[Variation] = IntervalTree()
        variation = Variation(start=0, end=1, strand=Strand.FORWARD)

        with caplog.at_level(logging.WARNING):
            interval_tree.add_features([variation])

        assert "Feature has invalid coordinates" in caplog.text
        assert interval_tree.features == [variation]
        assert interval_tree._intervals.size == 0


class TestIntervalTreeQueries:
    """Querying a tree by interval or by feature."""

    def test_get_overlapping_features(
        self, interval_tree: IntervalTree[Variation], variations_for_overlaps: list[Variation]
    ) -> None:
        """Return the features meeting an interval, in the order they were added."""
        assert interval_tree.get_overlapping_features(start=10, end=10) == variations_for_overlaps[:2]
        assert interval_tree.get_overlapping_features(start=11, end=15) == variations_for_overlaps[2:4]

    def test_get_overlapping_positions(self, interval_tree: IntervalTree[Variation]) -> None:
        """Return the positions of the matching features, ascending."""
        assert interval_tree.get_overlapping_positions(start=10, end=10) == [0, 1]

    def test_get_overlapping_features_returns_empty_list_without_matches(
        self, interval_tree: IntervalTree[Variation]
    ) -> None:
        """Return nothing when no feature meets the interval."""
        assert interval_tree.get_overlapping_features(start=1, end=8) == []
        assert interval_tree.get_overlapping_features(start=18, end=20) == []

    def test_get_overlapping_features_returns_empty_list_without_features(self) -> None:
        """Return nothing when nothing has been indexed."""
        interval_tree: IntervalTree[Variation] = IntervalTree()

        assert interval_tree.get_overlapping_features(start=1, end=1) == []

    def test_get_overlapping_features_by_feature(
        self, interval_tree: IntervalTree[Variation], variations_for_overlaps: list[Variation]
    ) -> None:
        """Query with another feature's own span."""
        query = Feature(start=10, end=10, strand=Strand.FORWARD)

        assert interval_tree.get_overlapping_features_by_feature(query) == variations_for_overlaps[:2]

    def test_get_nearest_features(
        self,
        nearest_interval_tree: IntervalTree[Variation],
        variations_for_nearest_interval: list[Variation],
    ) -> None:
        """Return the nearest features to an interval."""
        assert nearest_interval_tree.get_nearest_features(start=13, end=13, k=1) == [
            variations_for_nearest_interval[2]
        ]
        assert (
            nearest_interval_tree.get_nearest_features(start=12, end=12, k=2)
            == variations_for_nearest_interval[:2]
        )

    def test_get_nearest_features_by_feature(
        self,
        nearest_interval_tree: IntervalTree[Variation],
        variations_for_nearest_interval: list[Variation],
    ) -> None:
        """Query proximity with another feature's own span."""
        query = Feature(start=13, end=13, strand=Strand.FORWARD)

        assert nearest_interval_tree.get_nearest_features_by_feature(query, k=1) == [
            variations_for_nearest_interval[2]
        ]

    @pytest.mark.parametrize(
        "method_name",
        ["get_overlapping_features", "get_overlapping_positions", "get_nearest_features"],
    )
    def test_query_methods_reject_an_interval_that_is_not_1_based(self, method_name: str) -> None:
        """Refuse a query that starts before the first base."""
        interval_tree: IntervalTree[Variation] = IntervalTree()
        method = getattr(interval_tree, method_name)

        with pytest.raises(ValueError, match="start must be 1-based"):
            method(start=0, end=1)

    def test_reversed_interval_is_read_as_an_insertion_point(
        self, interval_tree: IntervalTree[Variation], variations_for_overlaps: list[Variation]
    ) -> None:
        """Swap a reversed query interval rather than refusing it."""
        assert interval_tree.get_overlapping_features(start=10, end=9) == variations_for_overlaps[:2]

    def test_round_trips_through_pickle(self) -> None:
        """Survive a process boundary, leaving the lock behind.

        The features are plain intervals: a ``Feature`` holds a lock of its own and never pickles.
        """
        held = _Span(start=10, end=20)
        interval_tree: IntervalTree[_Span] = IntervalTree(features=[held])

        restored = pickle.loads(pickle.dumps(interval_tree))

        assert restored.get_overlapping_features(start=15, end=15) == [held]


class TestFeatureIntervalTreeCache:
    """Sharing one tree per region and source."""

    def test_cache_implements_the_cache_protocol(self) -> None:
        """Satisfy the cache contract."""
        assert isinstance(FeatureIntervalTreeCache(), IntervalTreeCache)

    def test_reuses_tree_for_same_identity(self) -> None:
        """Serve identical requests from one tree."""
        source = _DummyFeatureSource(source_type=SourceType.FILE)
        cache = FeatureIntervalTreeCache()

        assert cache.get_interval_tree("chr1", 1, 100, source) is cache.get_interval_tree(
            "chr1", 1, 100, source
        )

    def test_key_distinguishes_region_coordinates(self) -> None:
        """Make region name and coordinates part of tree identity."""
        source = _DummyFeatureSource(source_type=SourceType.FILE)
        cache = FeatureIntervalTreeCache()

        trees = [
            cache.get_interval_tree("chr1", 1, 100, source),
            cache.get_interval_tree("chr2", 1, 100, source),
            cache.get_interval_tree("chr1", 2, 100, source),
            cache.get_interval_tree("chr1", 1, 101, source),
        ]

        assert len({id(tree) for tree in trees}) == 4

    def test_key_distinguishes_sources(self) -> None:
        """Keep equal regions from different sources apart."""
        cache = FeatureIntervalTreeCache()
        first_source = _DummyFeatureSource(source_type=SourceType.FILE)
        second_source = _DummyFeatureSource(source_type=SourceType.FILE)

        first = cache.get_interval_tree("chr1", 1, 100, first_source)
        second = cache.get_interval_tree("chr1", 1, 100, second_source)

        assert first is not second

    def test_clear_stops_reusing_a_cached_tree(self) -> None:
        """Drop entries even while a consumer still holds a tree."""
        source = _DummyFeatureSource(source_type=SourceType.FILE)
        cache = FeatureIntervalTreeCache()
        first = cache.get_interval_tree("chr1", 1, 100, source)

        cache.clear()

        assert cache.get_interval_tree("chr1", 1, 100, source) is not first

    @pytest.mark.parametrize(("start", "end"), [(0, 1), (2, 1)])
    def test_rejects_invalid_coordinates(self, start: int, end: int) -> None:
        """Refuse a span that is not a 1-based ordered interval."""
        source = _DummyFeatureSource(source_type=SourceType.FILE)

        with pytest.raises(ValueError, match="slice start"):
            FeatureIntervalTreeCache().get_interval_tree("chr1", start, end, source)

    def test_concurrent_requests_share_one_tree(self) -> None:
        """Converge concurrent misses on one canonical tree."""
        source = _DummyFeatureSource(source_type=SourceType.FILE)
        cache = FeatureIntervalTreeCache()

        with ThreadPoolExecutor(max_workers=4) as executor:
            trees = list(
                executor.map(lambda _index: cache.get_interval_tree("chr1", 1, 100, source), range(20))
            )

        assert all(tree is trees[0] for tree in trees)


class TestIntervalTreeAddMethods:
    """Tests methods for adding features."""

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

    def test_add_features_warnings(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test that ``add_features`` warns on an invalid interval"""
        interval_tree: IntervalTree = IntervalTree()
        variation = Variation(start=0, end=1, strand=Strand.FORWARD)
        with caplog.at_level(logging.WARNING):
            interval_tree.add_features([variation])
            assert "Feature has invalid coordinates" in caplog.text, "Incorrect warning on invalid interval"

        # the feature is still added to the ``features`` attribute, but not to the ``_intervals``
        # easier than altering the indexes
        assert interval_tree.features == [variation]
        assert interval_tree._intervals.size == 0  # pylint: disable=protected-access


class TestIntervalTreeHelperMethods:
    """Tests for internal helper methods."""

    def test_check_interval(self) -> None:
        """Test the ``_check_interval`` method."""
        interval_tree: IntervalTree = IntervalTree()
        variation = Variation(start=0, end=1, strand=Strand.FORWARD)
        # check raises error for start < 0
        with pytest.raises(ValueError, match="start must be 1-based"):
            interval_tree._check_interval(  # pylint: disable=protected-access
                start=variation.start, end=variation.end
            )

        # check start and end are reversed if start > end
        insertion = Variation(start=2, end=1, strand=Strand.FORWARD)
        start, end = interval_tree._check_interval(  # pylint: disable=protected-access
            start=insertion.start, end=insertion.end
        )
        assert start == insertion.end
        assert end == insertion.start

    def test_add_features_to_intervals(self, variations_for_overlaps: list[Variation]) -> None:
        """Test the ``_add_features_to_intervals`` method"""
        interval_tree: IntervalTree = IntervalTree()
        interval_tree._add_features_to_intervals(  # pylint: disable=protected-access
            features=variations_for_overlaps, start_index=0
        )
        assert interval_tree._intervals.size == len(  # pylint: disable=protected-access
            variations_for_overlaps
        )


class TestIntervalTreeQueryMethods:
    """Tests methods for querying the tree."""

    def test_get_overlapping_features(
        self, interval_tree: IntervalTree, variations_for_overlaps: list[Variation]
    ) -> None:
        """Test the ``get_overlapping_features`` method."""
        overlapping_vars = interval_tree.get_overlapping_features(start=10, end=10)
        expected_overlaps = variations_for_overlaps[:2]
        assert len(overlapping_vars) == len(expected_overlaps)
        assert overlapping_vars == expected_overlaps

        overlapping_vars = interval_tree.get_overlapping_features(start=11, end=15)
        expected_overlaps = variations_for_overlaps[2:4]
        assert len(overlapping_vars) == len(expected_overlaps)
        assert overlapping_vars == expected_overlaps

    def test_get_overlapping_features_returns_empty_list_without_matches(
        self, interval_tree: IntervalTree
    ) -> None:
        """Test the ``get_overlapping_features`` method when no intervals overlap."""
        assert interval_tree.get_overlapping_features(start=1, end=8) == []
        assert interval_tree.get_overlapping_features(start=18, end=20) == []

    def test_get_overlapping_features_returns_empty_list_without_features(self) -> None:
        """Test the ``get_overlapping_features`` method before any intervals are added."""
        interval_tree: IntervalTree = IntervalTree()
        assert interval_tree.get_overlapping_features(start=1, end=1) == []

    def test_get_overlapping_features_by_feature(
        self, interval_tree: IntervalTree, variations_for_overlaps: list[Variation]
    ) -> None:
        """Test the ``get_overlapping_features_by_feature`` method."""
        query_feature = Feature(start=10, end=10, strand=Strand.FORWARD)
        overlapping_vars = interval_tree.get_overlapping_features_by_feature(query_feature)
        assert overlapping_vars == variations_for_overlaps[:2]

    def test_get_nearest_features(
        self, interval_tree_nearest_intervals: IntervalTree, variations_for_nearest_interval: list[Variation]
    ) -> None:
        """Test the ``get_nearest_features`` method."""
        nearest = interval_tree_nearest_intervals.get_nearest_features(start=13, end=13, k=1)
        assert nearest == [variations_for_nearest_interval[2]]

        nearest = interval_tree_nearest_intervals.get_nearest_features(start=12, end=12, k=2)
        assert nearest == variations_for_nearest_interval[:2]

    def test_get_nearest_features_by_feature(
        self, interval_tree_nearest_intervals: IntervalTree, variations_for_nearest_interval: list[Variation]
    ) -> None:
        """Test the ``get_nearest_features_by_feature`` method."""
        query_feature = Feature(start=13, end=13, strand=Strand.FORWARD)
        nearest = interval_tree_nearest_intervals.get_nearest_features_by_feature(query_feature, k=1)
        assert nearest == [variations_for_nearest_interval[2]]

    @pytest.mark.parametrize(
        ("method_name", "start", "end", "message"),
        [
            ("get_overlapping_features", 0, 1, "start must be 1-based"),
            ("get_nearest_features", 0, 1, "start must be 1-based"),
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


@pytest.fixture
def interval_tree_nearest_intervals(variations_for_nearest_interval: list[Variation]) -> IntervalTree:
    """Create an IntervalTree with the ``variations_for_nearest_interval`` fixture already loaded"""
    return IntervalTree(features=variations_for_nearest_interval)
