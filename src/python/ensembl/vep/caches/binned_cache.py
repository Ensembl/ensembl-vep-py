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
"""UCSC-style binned feature cache."""

from __future__ import annotations

from collections.abc import Hashable, Iterable
from typing import Generic, Protocol, TypeVar, runtime_checkable

from attrs import define, field, frozen

from ensembl.vep.caches.interval_tree import Interval

__all__ = ["BinnedCache", "BinnedFeatureCache"]

ObjectT = TypeVar("ObjectT")
_BIN_LEVELS = ((1, 26), (9, 23), (73, 20), (585, 17))


@runtime_checkable
class BinnedCache(Protocol[ObjectT]):
    """Cache of objects keyed by sequence region and interval."""

    def add(
        self,
        seq_region_name: str,
        value: ObjectT,
        interval: Interval,
        *,
        key: Hashable | None = None,
    ) -> None:
        """Cache one object under its interval."""

        raise NotImplementedError

    def get_objects(self, seq_region_name: str, intervals: Iterable[Interval]) -> list[ObjectT]:
        """Return cached objects that meet any query interval."""

        raise NotImplementedError

    def clear(self) -> None:
        """Clear cached state."""

        raise NotImplementedError

    def __len__(self) -> int:
        """Return how many objects the cache retains."""

        raise NotImplementedError


@frozen(kw_only=True, auto_attribs=True)
class _CachedFeature(Generic[ObjectT]):
    """Store one feature and its precise interval."""

    value: ObjectT
    interval: Interval


@define(kw_only=True, auto_attribs=True)
class BinnedFeatureCache(Generic[ObjectT]):
    """Source-feature cache using UCSC hierarchical bins.

    A feature is stored in its smallest enclosing bin. Cache membership is
    partitioned by sequence-region name. An optional stable key prevents a
    source from adding the same logical feature after overlapping fetches.
    """

    _bins: dict[str, dict[int, list[_CachedFeature[ObjectT]]]] = field(
        factory=dict,
        init=False,
        repr=False,
    )
    _keys: dict[str, set[Hashable]] = field(factory=dict, init=False, repr=False)

    def add(
        self,
        seq_region_name: str,
        value: ObjectT,
        interval: Interval,
        *,
        key: Hashable | None = None,
    ) -> None:
        """Cache a feature under its smallest enclosing bin.

        Args:
            seq_region_name: Sequence region containing the feature.
            value: Cached source feature.
            interval: Feature interval in sequence-region coordinates.
            key: Stable source identity used to ignore a repeated logical feature.
        """

        feature_key = id(value) if key is None else key
        region_keys = self._keys.setdefault(seq_region_name, set())
        if feature_key in region_keys:
            return
        region_keys.add(feature_key)
        region_bins = self._bins.setdefault(seq_region_name, {})
        cached_feature = _CachedFeature(value=value, interval=interval)
        bin_number = self._bin_for_interval(interval)
        region_bins.setdefault(bin_number, []).append(cached_feature)

    def get_objects(self, seq_region_name: str, intervals: Iterable[Interval]) -> list[ObjectT]:
        """Return cached candidates that meet any query interval."""

        region_bins = self._bins.get(seq_region_name)
        if region_bins is None:
            return []

        found: list[ObjectT] = []
        seen: set[int] = set()
        for interval in intervals:
            for bin_number in self._query_bins(interval):
                for cached_feature in region_bins.get(bin_number, []):
                    if not self._intersects(cached_feature.interval, interval):
                        continue
                    feature_id = id(cached_feature.value)
                    if feature_id in seen:
                        continue
                    seen.add(feature_id)
                    found.append(cached_feature.value)
        return found

    def clear(self) -> None:
        """Clear cached features."""

        self._bins.clear()
        self._keys.clear()

    def __len__(self) -> int:
        """Return the number of cached logical features."""

        return sum(len(keys) for keys in self._keys.values())

    @classmethod
    def _bin_for_interval(cls, interval: Interval) -> int:
        """Return the smallest UCSC bin enclosing an inclusive interval."""

        start = interval.lookup_start - 1
        end = interval.lookup_end - 1
        for offset, shift in reversed(_BIN_LEVELS):
            if start >> shift == end >> shift:
                return offset + (start >> shift)
        return 0

    @classmethod
    def _query_bins(cls, interval: Interval) -> tuple[int, ...]:
        """Return bins that may contain a query interval."""

        start = interval.lookup_start - 1
        end = interval.lookup_end - 1
        bins = [0]
        for offset, shift in _BIN_LEVELS:
            bins.extend(range(offset + (start >> shift), offset + (end >> shift) + 1))
        return tuple(bins)

    @staticmethod
    def _intersects(left: Interval, right: Interval) -> bool:
        """Return whether candidate-retrieval bounds intersect."""

        return left.lookup_start <= right.lookup_end and left.lookup_end >= right.lookup_start
