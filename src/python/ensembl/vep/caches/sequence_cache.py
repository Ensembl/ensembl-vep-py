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
"""Sequence cache."""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Callable
from threading import RLock
from typing import Protocol, TypeAlias, runtime_checkable

from attrs import define, field

SequenceSpanLoader: TypeAlias = Callable[[str, int, int], str | None]


# TBD: should this be an ABC instead of a Protocol? It would be more explicit, but less flexible.
@runtime_checkable
class SequenceCache(Protocol):
    """Interface for sequence caches."""

    def get_sequence(
        self,
        seq_region_name: str,
        start: int,
        end: int,
        loader: SequenceSpanLoader,
    ) -> str | None:
        """Return sequence for a span.

        Args:
            seq_region_name (str): The sequence region name.
            start (int): 1-based inclusive start coordinate.
            end (int): 1-based inclusive end coordinate.
            loader (SequenceSpanLoader): Loader used when sequence is not cached.

        Returns:
            str | None: The requested sequence, or None if unavailable.
        """

        raise NotImplementedError

    def clear(self) -> None:
        """Clear cached state."""

        raise NotImplementedError


@define(kw_only=True, auto_attribs=True)
class SequenceBinnedLRUCache:
    """An in-memory LRU cache for fixed-size sequence chunks (or bins).

    Attributes:
        enabled (bool): Whether the cache is enabled. Default is True.
        chunk_power (int): The power of two for the chunk size. Default is 18 (chunk size of 2^18 = 262,144).
        max_chunks (int): The maximum number of chunks to keep in the cache. Default is 1024.
    """

    enabled: bool = True
    chunk_power: int = 18
    max_chunks: int = 1024
    _cache: OrderedDict[tuple[str, int], str] = field(factory=OrderedDict, init=False, repr=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    @property
    def chunk_size(self) -> int:
        """Return the cache chunk size."""

        return 1 << self.chunk_power

    @property
    def max_span(self) -> int:
        """Return the largest sequence that can be stored in the cache."""

        return self.chunk_size * self.max_chunks

    def chunk_start(self, position: int) -> int:
        """Return the start coordinate of the chunk containing the given position.

        Position is 1-based inclusive. The returned chunk start is also 1-based inclusive.

        Args:
            position (int): position of sequence in coordinate space.

        Returns:
            int: chunk start coordinate.
        """

        return ((position - 1) >> self.chunk_power << self.chunk_power) + 1

    def get_sequence(
        self,
        seq_region_name: str,
        start: int,
        end: int,
        loader: SequenceSpanLoader,
    ) -> str | None:
        """Return sequence for a span using binned chunk caching when useful.

        Args:
            seq_region_name (str): The sequence region name.
            start (int): Start position of the sequence span.
            end (int): End position of the sequence span.
            loader (SequenceSpanLoader): Loader used when sequence is not cached.

        Returns:
            str | None: The requested sequence, or None if unavailable.
        """

        requested_length = end - start + 1
        if not self.enabled or requested_length > self.max_span:
            return loader(seq_region_name, start, end)

        chunk_start = self.chunk_start(start)
        sequence = ""

        for current_start in range(chunk_start, end + 1, self.chunk_size):
            current_end = min(end, current_start + self.chunk_size - 1)
            chunk = self._get_or_load(
                (seq_region_name, current_start),
                lambda: self._load_chunk(
                    seq_region_name,
                    current_start,  # pylint: disable=cell-var-from-loop
                    loader,
                ),
            )
            local_start = max(start, current_start) - current_start
            local_end = current_end - current_start + 1
            sequence += chunk[local_start:local_end]

        return sequence

    def _load_chunk(
        self,
        seq_region_name: str,
        chunk_start: int,
        loader: SequenceSpanLoader,
    ) -> str:
        """Load a full cached sequence chunk through the supplied span loader.

        Args:
            seq_region_name (str): The sequence region name.
            chunk_start (int): Start coordinate of the chunk.
            loader (SequenceSpanLoader): Loader used to retrieve uncached sequence.

        Returns:
            str: The loaded chunk sequence.
        """

        chunk_end = chunk_start + self.chunk_size - 1
        chunk = loader(seq_region_name, chunk_start, chunk_end)
        if chunk is None:
            return "N" * self.chunk_size

        return chunk

    def _get_or_load(self, key: tuple[str, int], loader: Callable[[], str]) -> str:
        """Return a cached chunk or populate it by calling loader.

        Args:
            key (tuple[str, int]): Cache key made of sequence region name and chunk start.
            loader (Callable[[], str]): Loader used on cache miss.

        Returns:
            str: The cached or newly loaded chunk.
        """

        if not self.enabled:
            return loader()

        with self._lock:
            cached = self._cache.get(key)
            if cached is not None:
                self._cache.move_to_end(key)
                return cached

        chunk = loader()

        with self._lock:
            self._cache[key] = chunk
            self._cache.move_to_end(key)
            while len(self._cache) > self.max_chunks:
                self._cache.popitem(last=False)

        return chunk

    def clear(self) -> None:
        """Clear the cache contents."""

        with self._lock:
            self._cache.clear()
