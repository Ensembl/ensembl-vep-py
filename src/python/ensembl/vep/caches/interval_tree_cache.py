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
"""IntervalTree cache."""

from __future__ import annotations

from collections import OrderedDict
from threading import RLock
from typing import Protocol, runtime_checkable

from attrs import define, field, frozen

from ensembl.vep.caches.interval_tree import IntervalTree
from ensembl.vep.sources.feature_source import FeatureSource


@frozen(kw_only=True, auto_attribs=True)
class IntervalTreeKey:
    """Identity of a cached IntervalTree object."""

    seq_region_name: str
    start: int
    end: int
    source_id: int | None


@runtime_checkable
class IntervalTreeCache(Protocol):
    """Interface for objects used to hold intervals in a data structure used to do overlaps."""

    def get_interval_tree(
        self,
        seq_region_name: str,
        start: int,
        end: int,
        source: FeatureSource,
    ) -> IntervalTree:
        """Return an object for looking up overlaps of intervals, constructing it on a cache miss."""

        raise NotImplementedError

    def clear(self) -> None:
        """Clear cached state."""

        raise NotImplementedError


@define(kw_only=True, auto_attribs=True)
class FeatureIntervalTreeCache:
    """Share interval tree objects.

    Currently we never remove a tree from the cache, but could update this if it is a memory hotspot.
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
        """Return an object for looking up overlaps of intervals.

        Args:
            seq_region_name: The sequence region name.
            start: 1-based inclusive start coordinate.
            end: 1-based inclusive end coordinate.
            source: The FeatureSource for the objects.
        Returns:
            The requested interval tree, or None if unavailable.
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
        """Clear all live cache entries."""

        with self._lock:
            self._cache.clear()
