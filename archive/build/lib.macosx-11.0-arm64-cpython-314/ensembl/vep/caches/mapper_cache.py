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
"""Byte-budget LRU cache for ``Mapper``'s O(cds_length) accelerator (``MapperIndex``).

``Mapper`` never builds or owns a ``MapperIndex`` itself, its own bisect index is already small
at every scale (~52 KiB even for a titin-scale CDS), so nothing requires a cache for correctness or
bounded memory. A ``MapperIndex`` only closes the remaining lookup gap between bisect (O(log exon count))
and a fully-expanded dict (O(1)), which matters on the two per-variant hot paths that call
``genomic_to_offset``/ ``offset_to_genomic`` in a loop (e.g., ``facts/position.py``'s ``_cds_footprint``,
``caches/codon_cache.py``'s codon-window builder).

Mirrors ``CodonEditLRUCache``'s shape with one deliberate difference:
``CodonEditLRUCache.max_entries`` is count-based, which works because every codon edit is roughly
the same size. A ``MapperIndex``'s size varies by orders of magnitude across transcripts (a few KB
for an ordinary gene, ~15 MB for a titin-scale one). A count cap does not bound memory here, so
eviction is driven by a byte budget instead, tracked via a running total.
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
    """Interface for ``MapperIndex`` caches."""

    def get_index(self, transcript: Transcript, loader: MapperIndexLoader) -> MapperIndex:
        """Returns the cached accelerator for ``transcript``, building it via ``loader`` on a miss.

        Args:
            transcript (Transcript): The transcript to key the cache on (by ``stable_id``).
            loader (MapperIndexLoader): Builds the index on a miss - typically
                ``transcript.mapper.build_index``. Not called on a hit.

        Returns:
            MapperIndex: The cached or newly built accelerator.
        """
        raise NotImplementedError

    def clear(self) -> None:
        """Clear cached state."""
        raise NotImplementedError


def _estimate_mapper_index_bytes(mapper_index: MapperIndex) -> int:
    """Estimate a ``MapperIndex``'s heap footprint: container overhead plus every boxed int.

    Args:
        mapper_index (MapperIndex): The accelerator to size.

    Returns:
        int: Estimated bytes, list and dict combined.
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
    """A byte-budget LRU cache for ``MapperIndex``, keyed on ``Transcript.stable_id``.

    Instance-scoped: construct one and pass it to every annotator that should share it.

    Keyed by ``stable_id`` rather than the ``Transcript`` object itself: ``Transcript`` defines a
    custom ``__eq__`` (on ``stable_id``/``version``) without ``unsafe_hash=True``, so attrs drops
    its auto-generated ``__hash__`` - it is not safely usable as a dict key.

    ``eq=False`` keeps identity semantics: this is a shared, mutable instance, not a value object,
    so two caches with equal counters must not compare equal.

    Attributes:
        enabled (bool): Whether the cache is enabled. Default is True.
        max_bytes (int): The maximum total estimated size of cached entries, in bytes. Default is
            64 MiB.
        hit_count (int): Number of get_index calls served from the cache so far.
        miss_count (int): Number of get_index calls that invoked ``loader`` so far.
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
        """Returns the cached accelerator for ``transcript``, building and storing it via
        ``loader`` on a miss.

        Args:
            transcript (Transcript): The transcript to key the cache on.
            loader (MapperIndexLoader): Builds the index on a miss. Not called on a hit.

        Returns:
            MapperIndex: The cached or newly built accelerator.
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
