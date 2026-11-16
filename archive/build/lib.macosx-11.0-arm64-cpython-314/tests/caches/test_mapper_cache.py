# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Unit tests for mapper index caches."""

from __future__ import annotations

from threading import Barrier, Thread

from ensembl.vep.caches.mapper_cache import (
    MapperIndexCache,
    MapperIndexLRUCache,
    _estimate_mapper_index_bytes,
)
from ensembl.vep.common import Strand
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.sequence import Mapper, MapperIndex


def _build_transcript(stable_id: str) -> Transcript:
    """Create a transcript for mapper cache testing."""
    return Transcript(start=1, end=1, strand=Strand.FORWARD, stable_id=stable_id)


def _build_mapper(exon_count: int, block_length: int) -> Mapper:
    """Create a mapper with a predictable number of equally sized blocks."""
    genomic_spans = [
        (1 + index * (block_length + 100), block_length + index * (block_length + 100))
        for index in range(exon_count)
    ]
    coding_end_offset = sum(end - start + 1 for start, end in genomic_spans) - 1
    return Mapper(
        genomic_spans=genomic_spans,
        strand=Strand.FORWARD,
        coding_start_offset=0,
        coding_end_offset=coding_end_offset,
    )


class TestMapperIndexLRUCache:
    """Tests for the mapper index LRU cache implementation."""

    def test_cache_implements_the_mapper_index_cache_interface(self) -> None:
        """Test the concrete cache satisfies the protocol at runtime."""
        cache = MapperIndexLRUCache()

        assert isinstance(cache, MapperIndexCache)

    def test_cached_index_avoids_reload(self) -> None:
        """Test a cached transcript does not invoke the loader again."""
        cache = MapperIndexLRUCache()
        transcript = _build_transcript("ENST_A")
        mapper = _build_mapper(exon_count=2, block_length=10)
        load_calls = 0

        def loader() -> MapperIndex:
            nonlocal load_calls
            load_calls += 1
            return mapper.build_index()

        first_index = cache.get_index(transcript, loader)
        second_index = cache.get_index(transcript, loader)

        assert first_index is second_index
        assert load_calls == 1
        assert cache.hit_count == 1
        assert cache.miss_count == 1

    def test_disabled_cache_always_loads(self) -> None:
        """Test a disabled cache bypasses storage entirely."""
        cache = MapperIndexLRUCache(enabled=False)
        transcript = _build_transcript("ENST_A")
        mapper = _build_mapper(exon_count=2, block_length=10)
        load_calls = 0

        def loader() -> MapperIndex:
            nonlocal load_calls
            load_calls += 1
            return mapper.build_index()

        cache.get_index(transcript, loader)
        cache.get_index(transcript, loader)

        assert load_calls == 2
        assert cache.hit_count == 0
        assert cache.miss_count == 2

    def test_clear_removes_cached_indexes_and_sizes(self) -> None:
        """Test clear removes cached indexes and their recorded sizes."""
        cache = MapperIndexLRUCache()
        mapper = _build_mapper(exon_count=4, block_length=20)
        transcript = _build_transcript("ENST_A")

        cache.get_index(transcript, mapper.build_index)
        assert cache._total_bytes > 0  # pylint: disable=protected-access

        cache.clear()

        assert cache._cache == {}  # pylint: disable=protected-access
        assert cache._sizes == {}  # pylint: disable=protected-access
        assert cache._total_bytes == 0  # pylint: disable=protected-access


class TestMapperIndexLRUCacheEviction:
    """Tests for byte-budget eviction from the mapper index cache."""

    def test_large_index_evicts_smaller_indexes(self) -> None:
        """Test one large index can evict multiple smaller indexes."""
        small_mapper = _build_mapper(exon_count=2, block_length=5)
        large_mapper = _build_mapper(exon_count=50, block_length=200)
        small_bytes = _estimate_mapper_index_bytes(small_mapper.build_index())
        large_bytes = _estimate_mapper_index_bytes(large_mapper.build_index())
        assert large_bytes > small_bytes * 10

        cache = MapperIndexLRUCache(max_bytes=int(small_bytes * 1.5))
        cache.get_index(_build_transcript("SMALL_1"), small_mapper.build_index)
        cache.get_index(_build_transcript("SMALL_2"), small_mapper.build_index)

        assert "SMALL_1" not in cache._cache  # pylint: disable=protected-access
        assert "SMALL_2" in cache._cache  # pylint: disable=protected-access

        cache.max_bytes = large_bytes
        cache.get_index(_build_transcript("LARGE"), large_mapper.build_index)

        assert "SMALL_2" not in cache._cache  # pylint: disable=protected-access
        assert "LARGE" in cache._cache  # pylint: disable=protected-access

    def test_least_recently_used_index_is_evicted_first(self) -> None:
        """Test a cache hit refreshes an index's eviction order."""
        mapper = _build_mapper(exon_count=3, block_length=10)
        one_index_bytes = _estimate_mapper_index_bytes(mapper.build_index())
        cache = MapperIndexLRUCache(max_bytes=int(one_index_bytes * 2.5))

        cache.get_index(_build_transcript("A"), mapper.build_index)
        cache.get_index(_build_transcript("B"), mapper.build_index)
        cache.get_index(_build_transcript("A"), mapper.build_index)
        cache.get_index(_build_transcript("C"), mapper.build_index)

        assert "A" in cache._cache  # pylint: disable=protected-access
        assert "B" not in cache._cache  # pylint: disable=protected-access
        assert "C" in cache._cache  # pylint: disable=protected-access

    def test_index_larger_than_budget_is_not_retained(self) -> None:
        """Test an index larger than the byte budget is immediately evicted."""
        cache = MapperIndexLRUCache(max_bytes=0)
        transcript = _build_transcript("ENST_A")
        mapper = _build_mapper(exon_count=1, block_length=3)

        first_index = cache.get_index(transcript, mapper.build_index)
        second_index = cache.get_index(transcript, mapper.build_index)

        assert first_index is not second_index
        assert cache.hit_count == 0
        assert cache.miss_count == 2
        assert cache._total_bytes == 0  # pylint: disable=protected-access


class TestMapperIndexLRUCacheRace:
    """Regression test for the double-count race in ``get_index``'s insert path."""

    def test_concurrent_miss_on_one_key_keeps_total_bytes_exact(self) -> None:
        """Two threads racing to build the same key never double-count its bytes."""
        cache = MapperIndexLRUCache()
        mapper = _build_mapper(exon_count=4, block_length=20)
        transcript = _build_transcript("RACE")
        barrier = Barrier(2)

        def loader() -> MapperIndex:
            barrier.wait(timeout=1)
            return mapper.build_index()

        errors: list[BaseException] = []

        def race() -> None:
            try:
                cache.get_index(transcript, loader)
            except BaseException as error:  # pylint: disable=broad-except
                errors.append(error)

        threads = [Thread(target=race) for _ in range(2)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=2)

        assert not errors, errors
        assert all(not thread.is_alive() for thread in threads)
        assert cache.miss_count == 2
        assert list(cache._cache) == ["RACE"]  # pylint: disable=protected-access
        assert cache._total_bytes == cache._sizes["RACE"]  # pylint: disable=protected-access
