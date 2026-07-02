"""Sequence cache helpers."""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Callable
from threading import RLock

from attrs import define, field


@define(kw_only=True, auto_attribs=True)
class SequenceChunkCache:
    """An in-memory LRU cache for fixed-size sequence chunks."""

    enabled: bool = True
    chunk_power: int = 12
    max_chunks: int = 1024
    _cache: OrderedDict[tuple[str, int], str] = field(factory=OrderedDict, init=False, repr=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    @property
    def chunk_size(self) -> int:
        """Return the cache chunk size."""

        return 1 << self.chunk_power

    def chunk_start(self, position: int) -> int:
        """Return the 1-based start coordinate of the chunk containing position.

        Args:
            position (int): 1-based sequence position.

        Returns:
            int: 1-based chunk start coordinate.
        """

        return ((position - 1) >> self.chunk_power << self.chunk_power) + 1

    def get_or_load(self, key: tuple[str, int], loader: Callable[[], str]) -> str:
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
