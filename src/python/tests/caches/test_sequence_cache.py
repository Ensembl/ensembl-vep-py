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
"""Unit tests for sequence caches."""

from ensembl.vep.caches import SequenceBinnedLRUCache, SequenceCache


class TestSequenceBinnedLRUCache:
    """Tests for the binned LRU sequence cache implementation."""

    def test_cache_implements_the_sequence_cache_interface(self) -> None:
        """Test the concrete cache satisfies the protocol at runtime."""

        cache = SequenceBinnedLRUCache()

        assert isinstance(cache, SequenceCache)

    def test_chunk_size_is_power_of_two(self) -> None:
        """Test chunk_size is derived from chunk_power as a power of two."""

        cache = SequenceBinnedLRUCache(chunk_power=3)

        assert cache.chunk_size == 8

    def test_max_span_uses_chunk_size_and_max_chunks(self) -> None:
        """Test max_span is derived from chunk_size times max_chunks."""

        cache = SequenceBinnedLRUCache(chunk_power=2, max_chunks=3)

        assert cache.chunk_size == 4
        assert cache.max_span == 12

    def test_chunk_start_uses_bitwise_alignment(self) -> None:
        """Test chunk_start aligns positions to the configured chunk boundary."""

        cache = SequenceBinnedLRUCache(chunk_power=2)

        assert cache.chunk_start(1) == 1
        assert cache.chunk_start(4) == 1
        assert cache.chunk_start(5) == 5
        assert cache.chunk_start(7) == 5
        assert cache.chunk_start(9) == 9
        assert cache.chunk_start(10) == 9

    def test_lru_eviction_discards_oldest_chunk(self) -> None:
        """Test the cache evicts the least-recently-used chunk."""

        cache = SequenceBinnedLRUCache(chunk_power=2, max_chunks=2)
        load_calls: list[tuple[str, int, int]] = []
        data = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

        def loader(seq_region_name: str, start: int, end: int) -> str:
            load_calls.append((seq_region_name, start, end))
            return data[start - 1 : end]

        assert cache.get_sequence("chr1", 1, 4, loader) == "ABCD"
        assert cache.get_sequence("chr1", 5, 8, loader) == "EFGH"
        assert cache.get_sequence("chr1", 1, 4, loader) == "ABCD"
        assert cache.get_sequence("chr1", 9, 12, loader) == "IJKL"
        assert cache.get_sequence("chr1", 1, 4, loader) == "ABCD"
        assert cache.get_sequence("chr1", 5, 8, loader) == "EFGH"

        # Evicts chunk 5-8 when requesting 9-12, so last call to 5-8 requires a reload.
        assert load_calls == [
            ("chr1", 1, 4),
            ("chr1", 5, 8),
            ("chr1", 9, 12),
            ("chr1", 5, 8),
        ]

    def test_disabled_cache_always_loads(self) -> None:
        """Test disabled caches bypass storage entirely."""

        cache = SequenceBinnedLRUCache(enabled=False)
        load_calls = 0

        def loader(_seq_region_name: str, start: int, end: int) -> str:
            nonlocal load_calls
            load_calls += 1
            return "ABCDEFGHIJKLMNOP"[start - 1 : end]

        assert cache.get_sequence("chr1", 2, 6, loader) == "BCDEF"
        assert cache.get_sequence("chr1", 2, 6, loader) == "BCDEF"
        assert load_calls == 2

    def test_large_sequence_requests_bypass_cache(self) -> None:
        """Test spans larger than the configured cache capacity load directly."""

        cache = SequenceBinnedLRUCache(chunk_power=2, max_chunks=2)
        load_calls: list[tuple[str, int, int]] = []

        def loader(seq_region_name: str, start: int, end: int) -> str:
            load_calls.append((seq_region_name, start, end))
            return "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[start - 1 : end]

        assert cache.get_sequence("chr1", 1, 10, loader) == "ABCDEFGHIJ"
        assert load_calls == [("chr1", 1, 10)]
        assert len(cache._cache) == 0  # pylint: disable=protected-access

    def test_clear_removes_all_cached_chunks(self) -> None:
        """Test clear() removes all cached chunks."""

        cache = SequenceBinnedLRUCache(chunk_power=2, max_chunks=2)
        load_calls: list[tuple[str, int, int]] = []

        def loader(seq_region_name: str, start: int, end: int) -> str:
            load_calls.append((seq_region_name, start, end))
            return "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[start - 1 : end]

        assert cache.get_sequence("chr1", 1, 4, loader) == "ABCD"
        assert cache.get_sequence("chr1", 5, 8, loader) == "EFGH"
        assert len(cache._cache) == 2  # pylint: disable=protected-access

        cache.clear()
        assert len(cache._cache) == 0  # pylint: disable=protected-access
