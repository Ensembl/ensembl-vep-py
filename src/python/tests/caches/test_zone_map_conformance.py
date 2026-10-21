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

# pylint: disable=redefined-outer-name
"""Conformance tests for zone-map fact preloading."""

from __future__ import annotations

import pytest

from ensembl.vep.annotators.consequence.facts.base import get_fact
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.caches.zone_map import SERVED_FACT_NAMES, ZoneMap, ZoneMapCache
from ensembl.vep.common import Strand
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import VariationAllele

_FLANK_BASES = 5000
_FIXTURE_TRANSCRIPT_STABLE_IDS = ("ENST0001", "ENST0002", "ENST0003", "ENST0004")
_SAMPLE_SPAN_LENGTHS = (2, 3, 5, 8, 17, 50)


def _create_single_base_transcript_variation_allele(
    transcript: Transcript, position: int
) -> TranscriptVariationAllele:
    """Create a single-base transcript variation allele at a genomic position."""
    variation_allele = VariationAllele(
        start=position,
        end=position,
        strand=Strand.FORWARD,
        slice=transcript.slice,
        reference_allele="A",
        alternate_allele="T",
    )
    return TranscriptVariationAllele(feature=transcript, overlapping_feature=variation_allele)


def _create_spanning_transcript_variation_allele(
    transcript: Transcript, span_start: int, span_end: int
) -> TranscriptVariationAllele:
    """Create a deletion spanning the supplied genomic coordinates."""
    variation_allele = VariationAllele(
        start=span_start,
        end=span_end,
        strand=Strand.FORWARD,
        slice=transcript.slice,
        reference_allele="A" * (span_end - span_start + 1),
        alternate_allele="",
    )
    return TranscriptVariationAllele(feature=transcript, overlapping_feature=variation_allele)


def _brute_force_served_facts(overlap: TranscriptVariationAllele) -> dict[str, bool]:
    """Return each served fact computed without a zone map."""
    fact_context = FactContext(overlap=overlap)
    return {fact_name: bool(get_fact(fact_context, fact_name)) for fact_name in SERVED_FACT_NAMES}


def _zone_map_served_facts(
    zone_map_annotator: ZoneMapCache, overlap: TranscriptVariationAllele
) -> dict[str, bool]:
    """Return each served fact prefilled by a zone map."""
    fact_context = FactContext(overlap=overlap)
    zone_map_annotator.prefill(fact_context)
    return {fact_name: bool(fact_context.memo[fact_name]) for fact_name in SERVED_FACT_NAMES}


def _flank_window(transcript: Transcript) -> tuple[int, int]:
    """Return the ``(window_start, window_end)`` genomic span to exhaustively check.

    Args:
        transcript (Transcript): The transcript to compute a window around.

    Returns:
        tuple[int, int]: The transcript's span extended by _FLANK_BASES on each side,
            clipped to the attached slice's own bounds.
    """
    assert transcript.slice is not None
    window_start = max(transcript.slice.start, transcript.start - _FLANK_BASES)
    window_end = min(transcript.slice.end, transcript.end + _FLANK_BASES)
    return window_start, window_end


class TestZoneMapExhaustiveSingleBasePositions:
    """Every base of each fixture transcript's span plus its 5kb flanks, one at a time."""

    @pytest.mark.parametrize("transcript_stable_id", _FIXTURE_TRANSCRIPT_STABLE_IDS)
    def test_every_position_matches_brute_force(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript], transcript_stable_id: str
    ) -> None:
        """Every position matches brute force."""
        transcript = synthetic_transcripts_by_stable_id[transcript_stable_id]
        zone_map_annotator = ZoneMapCache()
        window_start, window_end = _flank_window(transcript)

        mismatches: list[str] = []
        for position in range(window_start, window_end + 1):
            overlap = _create_single_base_transcript_variation_allele(transcript, position)
            expected_facts = _brute_force_served_facts(overlap)
            actual_facts = _zone_map_served_facts(zone_map_annotator, overlap)
            if actual_facts != expected_facts:
                mismatches.append(f"{transcript_stable_id}@{position}: {actual_facts} != {expected_facts}")
        assert not mismatches, "\n".join(mismatches[:20])

    def test_every_served_fact_name_is_exercised_at_least_once(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Sanity check on the test itself: every served fact goes True somewhere in the sweep.

        Guards against a vacuously-passing exhaustive test (e.g. a fact that is never True across
        every fixture transcript's window would still "match" trivially).
        """
        fact_names_seen_true: set[str] = set()
        for transcript_stable_id in _FIXTURE_TRANSCRIPT_STABLE_IDS:
            transcript = synthetic_transcripts_by_stable_id[transcript_stable_id]
            zone_map_annotator = ZoneMapCache()
            window_start, window_end = _flank_window(transcript)
            for position in range(window_start, window_end + 1):
                overlap = _create_single_base_transcript_variation_allele(transcript, position)
                actual_facts = _zone_map_served_facts(zone_map_annotator, overlap)
                fact_names_seen_true.update(name for name, value in actual_facts.items() if value)
        missing = set(SERVED_FACT_NAMES) - fact_names_seen_true
        assert not missing, f"served facts never observed True across any fixture transcript: {missing}"


class TestZoneMapCacheTinyBudgetConformance:
    """The array-backed cache must answer correctly even when a tiny byte budget forces eviction
    and rebuild on nearly every lookup, not just when everything fits."""

    _SAMPLE_COUNT = 150

    def test_alternating_transcripts_with_a_tiny_budget_matches_brute_force(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Alternating between two transcripts evicts and rebuilds the single cached entry each call."""
        first_transcript = synthetic_transcripts_by_stable_id[_FIXTURE_TRANSCRIPT_STABLE_IDS[0]]
        second_transcript = synthetic_transcripts_by_stable_id[_FIXTURE_TRANSCRIPT_STABLE_IDS[1]]
        # Small enough that the two transcripts' maps cannot both fit, so every switch between
        # them evicts and rebuilds -- but big enough that either one alone is kept, not evicted
        # to empty.
        one_map_bytes = max(
            ZoneMap.build(first_transcript).estimated_bytes(),
            ZoneMap.build(second_transcript).estimated_bytes(),
        )
        zone_map_annotator = ZoneMapCache(max_bytes=one_map_bytes + 1)

        first_window_start, first_window_end = _flank_window(first_transcript)
        second_window_start, second_window_end = _flank_window(second_transcript)
        first_stride = max(1, (first_window_end - first_window_start) // self._SAMPLE_COUNT)
        second_stride = max(1, (second_window_end - second_window_start) // self._SAMPLE_COUNT)
        first_positions = range(first_window_start, first_window_end + 1, first_stride)
        second_positions = range(second_window_start, second_window_end + 1, second_stride)

        mismatches: list[str] = []
        call_count = 0
        for first_position, second_position in zip(first_positions, second_positions):
            for transcript, position in (
                (first_transcript, first_position),
                (second_transcript, second_position),
            ):
                overlap = _create_single_base_transcript_variation_allele(transcript, position)
                expected_facts = _brute_force_served_facts(overlap)
                actual_facts = _zone_map_served_facts(zone_map_annotator, overlap)
                call_count += 1
                if actual_facts != expected_facts:
                    mismatches.append(
                        f"{transcript.stable_id}@{position}: {actual_facts} != {expected_facts}"
                    )

        assert not mismatches, "\n".join(mismatches[:20])
        assert zone_map_annotator.eviction_count > 0
        # Strict two-key alternation over a budget that fits only one map: every call rebuilds,
        # not just most of them.
        assert zone_map_annotator.miss_count == call_count


class TestZoneMapMultiBaseSpans:
    """A sample of multi-base spans (deletions), since a variant need not be a single base."""

    @pytest.mark.parametrize("transcript_stable_id", _FIXTURE_TRANSCRIPT_STABLE_IDS)
    def test_sample_of_spans_matches_brute_force(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript], transcript_stable_id: str
    ) -> None:
        """Sample of spans matches brute force."""
        transcript = synthetic_transcripts_by_stable_id[transcript_stable_id]
        zone_map_annotator = ZoneMapCache()
        window_start, window_end = _flank_window(transcript)

        mismatches: list[str] = []
        for span_length in _SAMPLE_SPAN_LENGTHS:
            stride = max(1, span_length * 37)
            for span_start in range(window_start, window_end - span_length + 1, stride):
                span_end = span_start + span_length - 1
                overlap = _create_spanning_transcript_variation_allele(transcript, span_start, span_end)
                expected_facts = _brute_force_served_facts(overlap)
                actual_facts = _zone_map_served_facts(zone_map_annotator, overlap)
                if actual_facts != expected_facts:
                    mismatches.append(
                        f"{transcript_stable_id}@{span_start}-{span_end}: {actual_facts} != {expected_facts}"
                    )
        assert not mismatches, "\n".join(mismatches[:20])
