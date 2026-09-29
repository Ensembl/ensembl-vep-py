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
"""Unit tests for the cross-payload transcript pool."""

from __future__ import annotations

from ensembl.vep.caches.feature_pool import (
    FeaturePool,
    TranscriptKey,
    TranscriptPool,
    _estimate_transcript_bytes,
)
from ensembl.vep.common import Strand
from ensembl.vep.models.genic import Exon, Transcript


def _build_transcript(
    *, stable_id: str = "ENST0001", start: int = 100, end: int = 200, exon_count: int = 1
) -> Transcript:
    """Create a transcript with ``exon_count`` equally sized exons, for pool size testing."""
    transcript = Transcript(start=start, end=end, strand=Strand.FORWARD, stable_id=stable_id)
    span = end - start + 1
    exon_length = max(1, span // exon_count)
    exons = []
    for index in range(exon_count):
        exon_start = start + index * exon_length
        exon_end = min(end, exon_start + exon_length - 1)
        exons.append(
            Exon(start=exon_start, end=exon_end, strand=Strand.FORWARD, stable_id=f"{stable_id}-{index}")
        )
    transcript.assign_without_hooks(spliced_exons=exons)
    return transcript


def _key(transcript: Transcript, *, region: str = "chr1") -> TranscriptKey:
    """Return the pool key for a transcript built by ``_build_transcript``."""
    return (region, transcript.stable_id, transcript.start, transcript.end)


class TestFeaturePoolConformance:
    """TranscriptPool is usable wherever a FeaturePool is expected."""

    def test_pool_satisfies_the_feature_pool_protocol(self) -> None:
        """Test the concrete pool conforms to the generic contract."""

        assert isinstance(TranscriptPool(), FeaturePool)

    def test_a_pool_typed_as_the_protocol_round_trips_a_transcript(self) -> None:
        """Test get/put/clear work through the protocol's own surface."""

        pool: FeaturePool[Transcript] = TranscriptPool()
        transcript = _build_transcript()
        key: TranscriptKey = ("1", transcript.stable_id, transcript.start, transcript.end)

        pool.put(key, transcript)

        assert pool.get(key) is transcript
        pool.clear()
        assert pool.get(key) is None


class TestGetPut:
    """Basic get/put behaviour."""

    def test_get_misses_on_empty_pool(self) -> None:
        """Test an empty pool reports a miss and counts it."""
        pool = TranscriptPool()

        assert pool.get(("chr1", "ENST0001", 100, 200)) is None
        assert pool.miss_count == 1
        assert pool.hit_count == 0

    def test_put_then_get_returns_the_same_object(self) -> None:
        """Test a stored transcript is returned by identity on the next get."""
        pool = TranscriptPool()
        transcript = _build_transcript()

        pool.put(_key(transcript), transcript)

        assert pool.get(_key(transcript)) is transcript
        assert pool.hit_count == 1
        assert len(pool) == 1

    def test_put_same_key_twice_does_not_double_count_bytes(self) -> None:
        """Test re-inserting the same key does not inflate the tracked byte total."""
        pool = TranscriptPool()
        transcript = _build_transcript()

        pool.put(_key(transcript), transcript)
        first_total = pool._total_bytes  # pylint: disable=protected-access
        pool.put(_key(transcript), transcript)
        second_total = pool._total_bytes  # pylint: disable=protected-access

        assert first_total == second_total
        assert len(pool) == 1

    def test_hit_moves_entry_to_most_recently_used(self) -> None:
        """Test a hit protects an entry from being the next eviction victim."""
        first = _build_transcript(stable_id="ENST0001", start=100, end=200)
        second = _build_transcript(stable_id="ENST0002", start=300, end=400)
        pool = TranscriptPool(max_bytes=_estimate(first) + _estimate(second))

        pool.put(_key(first), first)
        pool.put(_key(second), second)
        pool.get(_key(first))  # refresh first's recency
        third = _build_transcript(stable_id="ENST0003", start=500, end=600)
        pool.put(_key(third), third)

        assert pool.get(_key(first)) is first
        assert pool.get(_key(second)) is None
        assert pool.get(_key(third)) is third


class TestByteBudgetEviction:
    """Eviction driven by the byte budget."""

    def test_eviction_drops_oldest_entry_once_over_budget(self) -> None:
        """Test the least-recently-used entry is evicted first."""
        first = _build_transcript(stable_id="ENST0001", start=100, end=200)
        second = _build_transcript(stable_id="ENST0002", start=300, end=400)
        budget = _estimate(first) + _estimate(second) - 1
        pool = TranscriptPool(max_bytes=budget)

        pool.put(_key(first), first)
        pool.put(_key(second), second)

        assert pool.get(_key(first)) is None
        assert pool.get(_key(second)) is second
        assert pool.eviction_count == 1

    def test_single_entry_bigger_than_budget_is_kept_alone(self) -> None:
        """Test one oversized entry is never evicted to leave the pool empty."""
        huge = _build_transcript(stable_id="ENST0001", start=1, end=100_000, exon_count=50)
        pool = TranscriptPool(max_bytes=1)

        pool.put(_key(huge), huge)

        assert pool.get(_key(huge)) is huge
        assert len(pool) == 1

    def test_max_bytes_zero_disables_the_pool(self) -> None:
        """Test ``max_bytes <= 0`` makes every get a miss and every put a no-op."""
        pool = TranscriptPool(max_bytes=0)
        transcript = _build_transcript()

        pool.put(_key(transcript), transcript)

        assert pool.get(_key(transcript)) is None
        assert len(pool) == 0
        # Disabled short-circuits before counting -- a miss here would be misleading, since the
        # pool never attempted a lookup.
        assert pool.miss_count == 0


class TestRegionChange:
    """Dropping entries when the pooled region changes."""

    def test_inserting_a_different_region_drops_the_old_one(self) -> None:
        """Test only one region's entries are kept at a time."""
        pool = TranscriptPool()
        chr1_transcript = _build_transcript(stable_id="ENST0001", start=100, end=200)
        chr2_transcript = _build_transcript(stable_id="ENST0002", start=100, end=200)

        pool.put(
            ("chr1", chr1_transcript.stable_id, chr1_transcript.start, chr1_transcript.end), chr1_transcript
        )
        pool.put(
            ("chr2", chr2_transcript.stable_id, chr2_transcript.start, chr2_transcript.end), chr2_transcript
        )

        assert (
            pool.get(("chr1", chr1_transcript.stable_id, chr1_transcript.start, chr1_transcript.end)) is None
        )
        assert (
            pool.get(("chr2", chr2_transcript.stable_id, chr2_transcript.start, chr2_transcript.end))
            is chr2_transcript
        )
        assert pool.seq_region_eviction_count == 1
        assert len(pool) == 1

    def test_revisiting_a_dropped_region_is_just_a_miss(self) -> None:
        """Test a region can be pooled again later; the earlier drop is not an error state."""
        pool = TranscriptPool()
        chr1_transcript = _build_transcript(stable_id="ENST0001", start=100, end=200)
        chr2_transcript = _build_transcript(stable_id="ENST0002", start=100, end=200)

        pool.put(("chr1", "ENST0001", 100, 200), chr1_transcript)
        pool.put(("chr2", "ENST0002", 100, 200), chr2_transcript)
        pool.put(("chr1", "ENST0001", 100, 200), chr1_transcript)

        assert pool.get(("chr1", "ENST0001", 100, 200)) is chr1_transcript
        assert pool.get(("chr2", "ENST0002", 100, 200)) is None


class TestClear:
    """Resetting pool contents."""

    def test_clear_resets_entries_sizes_and_region(self) -> None:
        """Test clear empties the pool and forgets the pooled region."""
        pool = TranscriptPool()
        transcript = _build_transcript()
        pool.put(_key(transcript), transcript)

        pool.clear()

        assert len(pool) == 0
        assert pool._total_bytes == 0  # pylint: disable=protected-access
        assert pool._seq_region_name is None  # pylint: disable=protected-access
        # A transcript from a different region can now be inserted without triggering a region
        # eviction, proving the pooled region really was forgotten.
        other = _build_transcript(stable_id="ENST0002", start=1, end=50)
        pool.put(("chr9", "ENST0002", 1, 50), other)
        assert pool.seq_region_eviction_count == 0


def _estimate(transcript: Transcript) -> int:
    """Budget arithmetic in tests uses the pool's own size estimate, not a re-derived copy."""
    return _estimate_transcript_bytes(transcript)
