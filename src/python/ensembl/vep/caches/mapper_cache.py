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
"""Cache expanded mapper indexes used by consequence annotation.

``MapperIndex`` accelerates repeated coordinate conversion but varies substantially in size, so
the LRU cache uses a byte budget rather than an entry count.
"""

from __future__ import annotations

import sys
from collections import OrderedDict
from collections.abc import Callable
from threading import RLock
from typing import TYPE_CHECKING, Protocol, TypeAlias, runtime_checkable

from attrs import define, field

from ensembl.vep.models.sequence import MapperIndex

if TYPE_CHECKING:
    from ensembl.vep.models.genic import Transcript

__all__ = ["MapperIndexCache", "MapperIndexLRUCache", "MapperIndexLoader"]

MapperIndexLoader: TypeAlias = Callable[[], MapperIndex]


@runtime_checkable
class MapperIndexCache(Protocol):
    """Define storage used to reuse mapper indexes."""

    def get_index(self, transcript: Transcript, loader: MapperIndexLoader) -> MapperIndex:
        """Return the cached mapper index, loading it on a miss.

        Args:
            transcript (Transcript): Transcript used to key the cache by stable ID.
            loader (MapperIndexLoader): Loader used when the transcript is not cached.

        Returns:
            MapperIndex: Cached or newly loaded mapper index.
        """
        raise NotImplementedError

    def clear(self) -> None:
        """Clear cached state."""
        raise NotImplementedError


def _estimate_mapper_index_bytes(mapper_index: MapperIndex) -> int:
    """Estimate the heap space used by a mapper index.

    Args:
        mapper_index (MapperIndex): Mapper index to size.

    Returns:
        int: Estimated size in bytes.
    """
    offset_to_position = mapper_index.offset_to_position
    position_to_offset = mapper_index.position_to_offset
    list_bytes = sys.getsizeof(offset_to_position) + sum(
        sys.getsizeof(position) for position in offset_to_position
    )
    dict_bytes = sys.getsizeof(position_to_offset) + sum(
        sys.getsizeof(position) + sys.getsizeof(offset) for position, offset in position_to_offset.items()
    )
    return list_bytes + dict_bytes


@define(kw_only=True, auto_attribs=True, eq=False)
class MapperIndexLRUCache:
    """Store mapper indexes in an instance-scoped, byte-budget LRU cache.

    Entries use transcript stable IDs because ``Transcript`` is not hashable. ``eq=False``
    preserves identity semantics for mutable shared state.

    Attributes:
        enabled (bool): Whether caching is enabled.
        max_bytes (int): Maximum estimated size of cached entries in bytes.
        hit_count (int): Number of lookups served from the cache.
        miss_count (int): Number of lookups that invoked ``loader``.
    """

    enabled: bool = True
    max_bytes: int = 64 * 1024 * 1024
    hit_count: int = field(default=0, init=False)
    miss_count: int = field(default=0, init=False)
    _cache: OrderedDict[str, MapperIndex] = field(factory=OrderedDict, init=False, repr=False)
    _sizes: dict[str, int] = field(factory=dict, init=False, repr=False)
    _total_bytes: int = field(default=0, init=False, repr=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    def get_index(self, transcript: Transcript, loader: MapperIndexLoader) -> MapperIndex:
        """Return the cached mapper index, loading and storing it on a miss.

        Args:
            transcript (Transcript): Transcript used to key the cache.
            loader (MapperIndexLoader): Loader used when the transcript is not cached.

        Returns:
            MapperIndex: Cached or newly loaded mapper index.
        """
        if not self.enabled:
            self.miss_count += 1
            return loader()

        key = transcript.stable_id
        with self._lock:
            cached = self._cache.get(key)
            if cached is not None:
                self._cache.move_to_end(key)
                self.hit_count += 1
                return cached

        mapper_index = loader()
        size = _estimate_mapper_index_bytes(mapper_index)

        with self._lock:
            if key in self._sizes:
                # Another thread built and inserted the same key first; keep this build
                # (deterministic given the same transcript) but do not double-count its bytes.
                self._total_bytes -= self._sizes[key]
            self._cache[key] = mapper_index
            self._cache.move_to_end(key)
            self._sizes[key] = size
            self._total_bytes += size
            self.miss_count += 1
            while self._total_bytes > self.max_bytes and self._cache:
                evicted_key, _ = self._cache.popitem(last=False)
                self._total_bytes -= self._sizes.pop(evicted_key, 0)

        return mapper_index

    def clear(self) -> None:
        """Clear the cache contents."""
        with self._lock:
            self._cache.clear()
            self._sizes.clear()
            self._total_bytes = 0
