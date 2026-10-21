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
"""Unit tests for transcript zone maps."""

from __future__ import annotations

import tracemalloc
from threading import Barrier, Thread

import pytest

from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.caches import zone_map as zone_map_module
from ensembl.vep.caches.zone_map import MEMO_FACT_NAMES, SERVED_FACT_NAMES, ZoneMap, ZoneMapCache
from ensembl.vep.common import Strand
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import VariationAllele


def _build_transcript(*, stable_id: str = "ENST0001", version: int | None = None) -> Transcript:
    """Create a transcript for zone map cache testing."""
    transcript = Transcript(start=100, end=200, strand=Strand.FORWARD, stable_id=stable_id)
    transcript.version = version
    return transcript


class TestZoneMap:
    """Tests for looking up spans in a transcript zone map."""

    def test_lookup_returns_zero_outside_mapped_positions(self) -> None:
        """Test positions outside the zone map have no fact mask."""
        zone_map = ZoneMap(
            lowest_position=10,
            highest_position=29,
            starts=(10, 20),
            masks=(1, 2),
        )

        assert zone_map.lookup(1, 9) == 0
        assert zone_map.lookup(30, 40) == 0

    def test_lookup_combines_masks_across_zones(self) -> None:
        """Test a span crossing zone boundaries combines every overlapping mask."""
        zone_map = ZoneMap(
            lowest_position=10,
            highest_position=39,
            starts=(10, 20, 30),
            masks=(1, 2, 4),
        )

        assert zone_map.lookup(15, 25) == 3
        assert zone_map.lookup(5, 35) == 7

    def test_lookup_includes_fully_covered_inverted_interval(self) -> None:
        """Test a span must cover both bounds of an inverted interval to include its mask."""
        zone_map = ZoneMap(
            lowest_position=10,
            highest_position=19,
            starts=(10,),
            masks=(1,),
            inverted=((8, 30, 20),),
        )

        assert zone_map.lookup(20, 30) == 8
        assert zone_map.lookup(21, 29) == 0


class TestZoneMapCache:
    """Tests for building and reusing transcript zone maps."""

    def test_same_transcript_identifier_and_version_reuses_zone_map(self) -> None:
        """Test equivalent transcript models share one cached zone map."""
        cache = ZoneMapCache()
        first_transcript = _build_transcript()
        second_transcript = _build_transcript()

        first_zone_map = cache.zone_map(first_transcript)
        second_zone_map = cache.zone_map(second_transcript)

        assert first_zone_map is second_zone_map
        assert cache.maps_built == 1

    def test_different_transcript_version_builds_another_zone_map(self) -> None:
        """Test transcript versions use independent zone map cache entries."""
        cache = ZoneMapCache()
        first_transcript = _build_transcript(version=1)
        second_transcript = _build_transcript(version=2)

        first_zone_map = cache.zone_map(first_transcript)
        second_zone_map = cache.zone_map(second_transcript)

        assert first_zone_map is not second_zone_map
        assert cache.maps_built == 2

    def test_columns_returns_empty_bit_columns_for_no_contexts(self) -> None:
        """Test an empty batch returns a zero-valued column for every served fact."""
        cache = ZoneMapCache()

        columns = cache.columns([])

        assert tuple(columns) == SERVED_FACT_NAMES
        assert all(column == 0 for column in columns.values())
        assert cache.facts_served == 0

    def test_columns_prefills_required_facts_and_updates_counter(self) -> None:
        """Test column lookup also prefills facts required by other facts."""
        transcript = _build_transcript()
        variation_allele = VariationAllele(
            start=150,
            end=150,
            strand=Strand.FORWARD,
            reference_allele="A",
            alternate_allele="T",
        )
        transcript_variation_allele = TranscriptVariationAllele(
            feature=transcript,
            overlapping_feature=variation_allele,
        )
        fact_context = FactContext(overlap=transcript_variation_allele)
        cache = ZoneMapCache()

        columns = cache.columns([fact_context])

        assert columns["overlap_feature"] == 1
        assert set(fact_context.memo) == set(MEMO_FACT_NAMES)
        assert fact_context.memo["overlap_feature"] is True
        assert cache.facts_served == len(SERVED_FACT_NAMES)


class TestZoneMapCacheFlankOnlyOverlaps:
    """Test overlaps that reach a transcript only through its flank."""

    @staticmethod
    def _columns_for(transcript: Transcript, position: int) -> tuple[dict[str, int], ZoneMapCache]:
        """Return the columns of a single-base variant at a position, and the cache used."""
        variation_allele = VariationAllele(
            start=position,
            end=position,
            strand=Strand.FORWARD,
            reference_allele="A",
            alternate_allele="T",
        )
        overlap = TranscriptVariationAllele(feature=transcript, overlapping_feature=variation_allele)
        cache = ZoneMapCache()
        columns = cache.columns([FactContext(overlap=overlap)])
        return columns, cache

    @pytest.mark.parametrize(
        ("strand", "position", "expected_fact_name"),
        [
            pytest.param(Strand.FORWARD, 50, "is_upstream", id="forward-before"),
            pytest.param(Strand.FORWARD, 260, "is_downstream", id="forward-after"),
            pytest.param(Strand.REVERSE, 50, "is_downstream", id="reverse-before"),
            pytest.param(Strand.REVERSE, 260, "is_upstream", id="reverse-after"),
        ],
    )
    def test_flank_fact_is_served_without_building_a_map(
        self, strand: Strand, position: int, expected_fact_name: str
    ) -> None:
        """Test a flank-only variant gets its flank fact and no zone map is built."""
        transcript = Transcript(start=100, end=200, strand=strand, stable_id="ENST0001")

        columns, cache = self._columns_for(transcript, position)

        set_fact_names = {fact_name for fact_name, column in columns.items() if column}
        assert set_fact_names == {expected_fact_name}
        assert cache.maps_built == 0

    def test_variant_beyond_both_flanks_sets_no_fact(self) -> None:
        """Test a variant past the flank distance sets no positional fact."""
        transcript = Transcript(start=100_000, end=100_100, strand=Strand.FORWARD, stable_id="ENST0001")

        columns, cache = self._columns_for(transcript, 10)

        assert not any(columns.values())
        assert cache.maps_built == 0

    def test_variant_next_to_transcript_still_uses_the_map(self) -> None:
        """Test a variant within the margin of the transcript is answered from its zone map."""
        transcript = Transcript(start=100, end=200, strand=Strand.FORWARD, stable_id="ENST0001")

        columns, cache = self._columns_for(transcript, 95)

        assert columns["is_upstream"] == 1
        assert cache.maps_built == 1


class TestZoneMapEstimatedBytes:
    """Tests for the LRU byte-accounting estimate."""

    def test_estimated_bytes_grows_with_zone_count(self) -> None:
        """More zones (more starts) yield a larger estimate."""
        few_zones = ZoneMap(lowest_position=1, highest_position=10, starts=(1,), masks=(1,))
        many_zones = ZoneMap(
            lowest_position=1, highest_position=40, starts=(1, 11, 21, 31), masks=(1, 2, 4, 8)
        )

        assert many_zones.estimated_bytes() > few_zones.estimated_bytes()

    def test_estimated_bytes_grows_with_inverted_interval_count(self) -> None:
        """More inverted intervals also yield a larger estimate."""
        no_inverted = ZoneMap(lowest_position=1, highest_position=10, starts=(1,), masks=(1,))
        one_inverted = ZoneMap(
            lowest_position=1, highest_position=10, starts=(1,), masks=(1,), inverted=((2, 9, 5),)
        )

        assert one_inverted.estimated_bytes() > no_inverted.estimated_bytes()

    def test_estimate_is_within_25_percent_of_tracemalloc_for_a_mid_size_map(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A ~500-zone map's estimate stays within +/-25% of a tracemalloc-measured first build.

        Uses a fresh mask-intern cache so the measurement reflects a cold first build (the
        conservative, worst case), not one that benefits from another test's already-boxed ints.
        """
        monkeypatch.setattr(zone_map_module, "_MASK_INTERN_CACHE", {})
        zone_count = 500
        starts = list(range(1, zone_count * 3, 3))
        distinct_masks = [1 << bit for bit in range(9, 17)]
        masks = [distinct_masks[index % len(distinct_masks)] for index in range(zone_count)]

        tracemalloc.start()
        baseline = tracemalloc.take_snapshot()
        zone_map = ZoneMap(
            lowest_position=starts[0], highest_position=starts[-1] + 100, starts=starts, masks=masks
        )
        after = tracemalloc.take_snapshot()
        tracemalloc.stop()

        actual_bytes = sum(
            stat.size_diff for stat in after.compare_to(baseline, "lineno") if stat.size_diff > 0
        )
        estimated_bytes = zone_map.estimated_bytes()

        assert 0.75 * actual_bytes <= estimated_bytes <= 1.25 * actual_bytes


class TestZoneMapCacheByteEviction:
    """Tests for the byte-budgeted LRU behaviour of ZoneMapCache."""

    @staticmethod
    def _one_map_bytes() -> int:
        """Return the estimated size of the fixture transcript's zone map."""
        return ZoneMap.build(_build_transcript(stable_id="SIZE_PROBE")).estimated_bytes()

    def test_hit_refreshes_lru_order(self) -> None:
        """A cache hit moves its entry to the most-recently-used end, so it survives eviction."""
        one_map_bytes = self._one_map_bytes()
        cache = ZoneMapCache(max_bytes=int(one_map_bytes * 2.5))

        cache.zone_map(_build_transcript(stable_id="A"))
        cache.zone_map(_build_transcript(stable_id="B"))
        cache.zone_map(_build_transcript(stable_id="A"))  # refresh A's LRU position
        cache.zone_map(_build_transcript(stable_id="C"))  # over budget: evicts B, not A

        cached_keys = set(cache._zone_maps)  # pylint: disable=protected-access
        assert ("A", None) in cached_keys
        assert ("B", None) not in cached_keys
        assert ("C", None) in cached_keys
        assert cache.eviction_count == 1

    def test_eviction_drops_the_oldest_entry_and_keeps_the_newest(self) -> None:
        """Exceeding the budget evicts from the front; the map just built survives."""
        one_map_bytes = self._one_map_bytes()
        cache = ZoneMapCache(max_bytes=int(one_map_bytes * 1.5))

        cache.zone_map(_build_transcript(stable_id="OLD"))
        cache.zone_map(_build_transcript(stable_id="NEW"))

        cached_keys = set(cache._zone_maps)  # pylint: disable=protected-access
        assert ("OLD", None) not in cached_keys
        assert ("NEW", None) in cached_keys
        assert cache.eviction_count == 1
        assert cache._total_bytes == one_map_bytes  # pylint: disable=protected-access

    def test_single_map_larger_than_budget_is_kept_alone(self) -> None:
        """A map that alone exceeds the budget stays cached rather than being evicted to empty."""
        cache = ZoneMapCache(max_bytes=1)

        built = cache.zone_map(_build_transcript(stable_id="ONLY"))
        reused = cache.zone_map(_build_transcript(stable_id="ONLY"))

        assert reused is built
        assert list(cache._zone_maps) == [("ONLY", None)]  # pylint: disable=protected-access
        assert cache.eviction_count == 0
        assert cache.hit_count == 1

    def test_clear_resets_entries_sizes_and_total(self) -> None:
        """clear() empties every tracked entry, size and byte total."""
        cache = ZoneMapCache()
        cache.zone_map(_build_transcript(stable_id="A"))
        assert cache._total_bytes > 0  # pylint: disable=protected-access

        cache.clear()

        assert cache._zone_maps == {}  # pylint: disable=protected-access
        assert cache._sizes == {}  # pylint: disable=protected-access
        assert cache._total_bytes == 0  # pylint: disable=protected-access

    def test_concurrent_miss_on_one_key_keeps_total_bytes_exact(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Two threads racing to build the same key never double-count its bytes."""
        cache = ZoneMapCache()
        barrier = Barrier(2)
        original_build = ZoneMap.build

        def synchronized_build(_cls: type[ZoneMap], transcript: Transcript, *distances: int) -> ZoneMap:
            barrier.wait(timeout=1)
            return original_build(transcript, *distances)

        monkeypatch.setattr(ZoneMap, "build", classmethod(synchronized_build))

        errors: list[BaseException] = []

        def race() -> None:
            try:
                cache.zone_map(_build_transcript(stable_id="RACE"))
            except BaseException as error:  # pylint: disable=broad-except
                errors.append(error)

        threads = [Thread(target=race), Thread(target=race)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=2)

        assert not errors, errors
        assert all(not thread.is_alive() for thread in threads)
        key = ("RACE", None)
        assert cache.miss_count == 2
        assert list(cache._zone_maps) == [key]  # pylint: disable=protected-access
        assert cache._total_bytes == cache._sizes[key]  # pylint: disable=protected-access
