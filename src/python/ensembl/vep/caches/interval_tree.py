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
"""Augmented interval-list feature indexing.

``IntervalTree`` wraps ``ailist``, so interval intersection runs in C over the whole batch rather
than in Python per feature. ``FeatureIntervalTreeCache`` shares one tree per region and source,
the way ``zone_map.py`` pairs a structure with the cache that reuses it. ``Interval`` is the query
value the rest of the package builds, and keeps the insertion-point representation the matchers
depend on.
"""

from __future__ import annotations

import logging
from collections import OrderedDict
from threading import RLock
from typing import TYPE_CHECKING, Any, Generic, Protocol, TypeVar, runtime_checkable

from ailist import AIList
from attrs import define, field, frozen

if TYPE_CHECKING:
    # Both of these packages import ``caches``; keeping their imports out of runtime keeps that
    # one-way, so nothing that holds an interval has to be imported to build one.
    from ensembl.vep.matcher.base import IntervalLike
    from ensembl.vep.sources.feature_source import FeatureSource

__all__ = [
    "FeatureIntervalTreeCache",
    "Interval",
    "IntervalTree",
    "IntervalTreeCache",
    "IntervalTreeKey",
]

logger = logging.getLogger(__name__)

FeatureT = TypeVar("FeatureT", bound="IntervalLike")
"""The feature type an interval tree holds; anything with inclusive coordinates."""


@frozen(kw_only=True, auto_attribs=True)
class Interval:
    """Inclusive genomic interval or interbase insertion point.

    ``start=end+1`` denotes an insertion before ``start``. The interval index
    expands that point to its two flanking bases only for candidate retrieval.
    Matchers retain the zero-width insertion semantics.

    Attributes:
        start (int): Inclusive start coordinate.
        end (int): Inclusive end coordinate, or ``start - 1`` for an insertion.
    """

    start: int
    end: int

    def __attrs_post_init__(self) -> None:
        """Validate the interval representation."""

        if self.start < 1:
            raise ValueError("interval start must be at least 1")
        if self.end < self.start - 1:
            raise ValueError("interval end must not be less than start - 1")

    @property
    def is_insertion(self) -> bool:
        """Return whether the interval is a zero-width insertion point."""

        return self.start == self.end + 1

    @property
    def lookup_start(self) -> int:
        """Return the inclusive start used for candidate retrieval."""

        return max(1, self.end) if self.is_insertion else self.start

    @property
    def lookup_end(self) -> int:
        """Return the inclusive end used for candidate retrieval."""

        return self.start if self.is_insertion else self.end


@define(kw_only=True, auto_attribs=True)
class IntervalTree(Generic[FeatureT]):
    """Feature index for overlap and proximity queries.

    Attributes:
        features (list[FeatureT]): The indexed features, in the order they were added.
    """

    features: list[FeatureT] = field(factory=list)
    _intervals: AIList = field(factory=AIList, init=False, repr=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    def __attrs_post_init__(self) -> None:
        """Index the features supplied at construction."""

        self._add_features_to_intervals(self.features, start_index=0)

    def __getstate__(self) -> dict[str, Any]:
        """Return the tree's state without the lock, which no lock can be pickled through."""

        return {"features": self.features, "_intervals": self._intervals}

    def __setstate__(self, state: dict[str, Any]) -> None:
        """Restore the tree with a lock of its own."""

        object.__setattr__(self, "features", state["features"])
        object.__setattr__(self, "_intervals", state["_intervals"])
        object.__setattr__(self, "_lock", RLock())

    def add_features(self, features: list[FeatureT]) -> None:
        """Add features to the index.

        Args:
            features (list[FeatureT]): The features to add.
        """

        self._add_features_to_intervals(features, start_index=len(self.features))
        self.features.extend(features)

    def get_overlapping_features(self, start: int, end: int) -> list[FeatureT]:
        """Return every indexed feature overlapping an inclusive interval.

        Args:
            start (int): The 1-based inclusive start coordinate.
            end (int): The inclusive end coordinate.

        Returns:
            list[FeatureT]: The overlapping features, in the order they were added.

        Raises:
            ValueError: If the start coordinate is not 1-based.
        """

        return [self.features[index] for index in self.get_overlapping_positions(start, end)]

    def get_overlapping_features_by_feature(self, feature: IntervalLike) -> list[FeatureT]:
        """Return every indexed feature overlapping another feature's span.

        Args:
            feature (IntervalLike): The feature whose span is queried.

        Returns:
            list[FeatureT]: The overlapping features, in the order they were added.

        Raises:
            ValueError: If the feature's coordinates are not 1-based.
        """

        return self.get_overlapping_features(start=feature.start, end=feature.end)

    def get_overlapping_positions(self, start: int, end: int) -> list[int]:
        """Return the positions in ``features`` of every feature overlapping an interval.

        A caller querying several intervals deduplicates and orders positions far more cheaply
        than the features themselves, which is why the positions are public.

        Args:
            start (int): The 1-based inclusive start coordinate.
            end (int): The inclusive end coordinate.

        Returns:
            list[int]: The matching positions, ascending.

        Raises:
            ValueError: If the start coordinate is not 1-based.
        """

        with self._lock:
            start, end = self._check_interval(start, end)
            return sorted(int(index) for index in self._intervals.intersect_ids(start, end + 1))

    def get_nearest_features(self, start: int, end: int, k: int = 10) -> list[FeatureT]:
        """Return the ``k`` indexed features nearest an inclusive interval.

        Args:
            start (int): The 1-based inclusive start coordinate.
            end (int): The inclusive end coordinate.
            k (int): How many features to return.

        Returns:
            list[FeatureT]: The nearest features, nearest first.

        Raises:
            ValueError: If the start coordinate is not 1-based.
        """

        with self._lock:
            start, end = self._check_interval(start, end)
            closest = self._intervals.closest(start, end + 1, k=k)
            return [self.features[index] for index in closest.ids]

    def get_nearest_features_by_feature(self, feature: IntervalLike, k: int = 10) -> list[FeatureT]:
        """Return the ``k`` indexed features nearest another feature's span.

        Args:
            feature (IntervalLike): The feature whose span is queried.
            k (int): How many features to return.

        Returns:
            list[FeatureT]: The nearest features, nearest first.

        Raises:
            ValueError: If the feature's coordinates are not 1-based.
        """

        return self.get_nearest_features(start=feature.start, end=feature.end, k=k)

    def _add_features_to_intervals(self, features: list[FeatureT], start_index: int) -> None:
        """Index each feature under its own coordinates, skipping any that are invalid.

        Args:
            features (list[FeatureT]): The features to index.
            start_index (int): The position the first feature is indexed under.
        """

        for offset, feature in enumerate(features):
            try:
                start, end = self._check_interval(feature.start, feature.end)
            except ValueError:
                logger.warning("Feature has invalid coordinates")
                continue
            self._intervals.add(start, end + 1, start_index + offset)

    @staticmethod
    def _check_interval(start: int, end: int) -> tuple[int, int]:
        """Return the interval as an ordered 1-based pair.

        Args:
            start (int): The 1-based start coordinate.
            end (int): The end coordinate.

        Returns:
            tuple[int, int]: The start and end, swapped when given reversed.

        Raises:
            ValueError: If the start coordinate is not 1-based.
        """

        if start < 1:
            raise ValueError("start must be 1-based: cannot be less than 1")
        if start > end:
            return end, start
        return start, end


@frozen(kw_only=True, auto_attribs=True)
class IntervalTreeKey:
    """Identity of one cached interval tree.

    Attributes:
        seq_region_name (str): The sequence region the tree covers.
        start (int): The 1-based inclusive start of the covered span.
        end (int): The inclusive end of the covered span.
        source_id (int | None): Identity of the feature source the tree was filled from.
    """

    seq_region_name: str
    start: int
    end: int
    source_id: int | None


@runtime_checkable
class IntervalTreeCache(Protocol):
    """Region-based interval-tree cache."""

    def get_interval_tree(
        self,
        seq_region_name: str,
        start: int,
        end: int,
        source: FeatureSource,
    ) -> IntervalTree:
        """Return the interval tree for one region and source."""

        raise NotImplementedError

    def clear(self) -> None:
        """Clear cached state."""

        raise NotImplementedError


@define(kw_only=True, auto_attribs=True)
class FeatureIntervalTreeCache:
    """Shared interval trees keyed by region and source.

    Nothing is evicted yet; the cache is bounded by how many regions a run touches.
    """

    _cache: OrderedDict[IntervalTreeKey, IntervalTree] = field(factory=OrderedDict, init=False, repr=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    def get_interval_tree(
        self,
        seq_region_name: str,
        start: int,
        end: int,
        source: FeatureSource,
    ) -> IntervalTree:
        """Return the interval tree for one region and source, building it on a miss.

        Args:
            seq_region_name (str): The sequence region name.
            start (int): The 1-based inclusive start coordinate.
            end (int): The inclusive end coordinate.
            source (FeatureSource): The source whose features the tree holds.

        Returns:
            IntervalTree: The cached or newly built tree.

        Raises:
            ValueError: If the span is not a 1-based ordered interval.
        """

        if start < 1:
            raise ValueError("slice start must be 1-based")
        if end < start:
            raise ValueError("slice start must be less than or equal to slice end")

        key = IntervalTreeKey(
            seq_region_name=seq_region_name,
            start=start,
            end=end,
            source_id=id(source),
        )
        with self._lock:
            cached = self._cache.get(key)
            if cached is not None:
                return cached
            tree: IntervalTree = IntervalTree()
            self._cache[key] = tree

        return tree

    def clear(self) -> None:
        """Clear every cached tree."""

        with self._lock:
            self._cache.clear()
