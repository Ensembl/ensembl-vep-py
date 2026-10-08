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
"""Unit tests for slice caches."""

import gc
from concurrent.futures import ThreadPoolExecutor
from weakref import ref

import pytest

from ensembl.vep.caches import SliceCache, WeakSliceCache
from ensembl.vep.common import Strand
from tests.models.conftest import DummySequenceSource


class TestWeakSliceCache:
    """Tests for the weak-value slice cache."""

    def test_cache_implements_slice_cache_interface(self) -> None:
        """Test the concrete cache satisfies the cache protocol."""

        assert isinstance(WeakSliceCache(), SliceCache)

    def test_uses_sequence_source_length_for_default_end(self) -> None:
        """Test a default slice spans the complete sequence region."""

        source = DummySequenceSource(data={"chr1": "ACGT"})
        seq_slice = WeakSliceCache().get_slice("chr1", source)

        assert seq_slice.start == 1
        assert seq_slice.end == 4
        assert seq_slice.strand == Strand.FORWARD
        assert seq_slice.sequence_source is source

    def test_requires_resolvable_end(self) -> None:
        """Test slices cannot be built without an explicit or resolvable end."""

        cache = WeakSliceCache()

        with pytest.raises(ValueError, match="no end or sequence source"):
            cache.get_slice("chr1")
        with pytest.raises(ValueError, match="sequence-region length is unknown"):
            cache.get_slice("chr1", DummySequenceSource())

    def test_explicit_end_does_not_require_sequence_source(self) -> None:
        """Test explicit coordinates are sufficient to construct a slice."""

        seq_slice = WeakSliceCache().get_slice("chr1", start=10, end=20)

        assert seq_slice.start == 10
        assert seq_slice.end == 20
        assert seq_slice.sequence_source is None

    def test_reuses_slice_for_same_identity(self) -> None:
        """Test identical requests share a canonical slice."""

        source = DummySequenceSource(data={"chr1": "ACGT"})
        cache = WeakSliceCache()

        first = cache.get_slice("chr1", source)
        second = cache.get_slice("chr1", source)

        assert first is second

    def test_key_distinguishes_coordinates_and_strand(self) -> None:
        """Test range and strand are part of slice identity."""

        source = DummySequenceSource(data={"chr1": "A" * 100})
        cache = WeakSliceCache()

        forward = cache.get_slice("chr1", source, start=1, end=50)
        shifted = cache.get_slice("chr1", source, start=2, end=50)
        shorter = cache.get_slice("chr1", source, start=1, end=49)
        reverse = cache.get_slice("chr1", source, start=1, end=50, strand=Strand.REVERSE)

        assert len({id(forward), id(shifted), id(shorter), id(reverse)}) == 4

    def test_key_distinguishes_sequence_sources(self) -> None:
        """Test equal region names from different sequence sources do not collide."""

        first_source = DummySequenceSource(data={"chr1": "ACGT"})
        second_source = DummySequenceSource(data={"chr1": "ACGT"})
        cache = WeakSliceCache()

        first = cache.get_slice("chr1", first_source)
        second = cache.get_slice("chr1", second_source)

        assert first is not second
        assert first.sequence_source is first_source
        assert second.sequence_source is second_source

    def test_releases_slice_without_external_references(self) -> None:
        """Test the cache does not extend a slice's lifetime."""

        source = DummySequenceSource(data={"chr1": "ACGT"})
        cache = WeakSliceCache()
        seq_slice = cache.get_slice("chr1", source)
        slice_reference = ref(seq_slice)

        del seq_slice
        gc.collect()

        assert slice_reference() is None

    def test_clear_stops_reusing_live_slice(self) -> None:
        """Test clear removes entries even when consumers retain the slices."""

        source = DummySequenceSource(data={"chr1": "ACGT"})
        cache = WeakSliceCache()
        first = cache.get_slice("chr1", source)

        cache.clear()
        second = cache.get_slice("chr1", source)

        assert first is not second

    @pytest.mark.parametrize(("start", "end"), [(0, 1), (2, 1)])
    def test_rejects_invalid_coordinates(self, start: int, end: int) -> None:
        """Test invalid slice coordinates are rejected."""

        with pytest.raises(ValueError, match="slice start"):
            WeakSliceCache().get_slice("chr1", start=start, end=end)

    def test_concurrent_requests_share_slice(self) -> None:
        """Test concurrent misses converge on one canonical slice."""

        source = DummySequenceSource(data={"chr1": "ACGT"})
        cache = WeakSliceCache()

        with ThreadPoolExecutor(max_workers=4) as executor:
            slices = list(executor.map(lambda _index: cache.get_slice("chr1", source), range(20)))

        assert all(seq_slice is slices[0] for seq_slice in slices)
