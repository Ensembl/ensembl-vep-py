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
"""The transcript consequence track and its candidate lookup.

The track holds the whole run path for transcript consequences: the transcript source, the caches
that make a payload's windows cheap, the pairing rule, the annotator, and the alleles the results
are attached to. Pairing and the fetch window are built by their own functions as well, because
explaining one derivation must use the same rule the run uses without wanting any of the caching
around it.
"""

from __future__ import annotations

import logging
from bisect import bisect_left, bisect_right
from collections.abc import Sequence
from typing import Any

from attrs import asdict, define, field

from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.caches import BinnedCache, BinnedFeatureCache, Interval, TranscriptPool
from ensembl.vep.common import Readthrough
from ensembl.vep.common.constants import UNKNOWN_VALUE
from ensembl.vep.config import Config
from ensembl.vep.matcher import OverlapMatcher
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import OverlapConsequence, TranscriptVariationAllele
from ensembl.vep.models.variation import Variation
from ensembl.vep.overlaps import CandidateLookup, TranscriptOverlapFinder
from ensembl.vep.overlaps.base import VariationAlleleBuffer
from ensembl.vep.sources.build import TranscriptSource, build_transcript_source
from ensembl.vep.workflow.tracks.base import (
    AnnotationTrack,
    OverlapAnnotationTrack,
    TrackResources,
    TrackResult,
)
from ensembl.vep.workflow.tracks.registry import register_track

__all__ = [
    "TranscriptCandidateLookup",
    "TranscriptConsequenceTrack",
    "build_transcript_overlap_finder",
    "build_transcript_tracks",
    "transcript_candidate_window",
]

logger = logging.getLogger(__name__)

TRACK_IDENTIFIER = "transcript"
"""The identifier the transcript consequence track is configured under."""


def transcript_candidate_window(config: Config) -> int:
    """Return the bases added each side of an allele before transcripts are fetched.

    The wider of the two distances: which one applies depends on each transcript's strand, and the
    overlap finder is what decides that once they are fetched.
    """
    return max(config.source.upstream_distance, config.source.downstream_distance)


def build_transcript_overlap_finder(
    config: Config,
    candidates: CandidateLookup[Transcript],
) -> TranscriptOverlapFinder:
    """Return the rule pairing transcripts with alleles, over a caller's own candidates.

    The one definition of which ``(allele, transcript)`` pairs exist. A run reaches it through its
    track, which fetches candidates from a cached source; ``vep explain`` reaches it with the
    transcripts it already found, so an explanation is never a second opinion about the pairing.

    Args:
        config (Config): The run's configuration, for the distances and readthrough policy.
        candidates (CandidateLookup[Transcript]): Return the transcripts an input buffer may pair with.

    Returns:
        TranscriptOverlapFinder: The configured finder.
    """
    readthrough = config.annotation.readthrough
    return TranscriptOverlapFinder(
        candidates=candidates,
        matcher=OverlapMatcher(),
        readthrough=Readthrough[readthrough.upper()] if isinstance(readthrough, str) else readthrough,
        upstream_distance=config.source.upstream_distance,
        downstream_distance=config.source.downstream_distance,
    )


@register_track("consequence.transcript")
def build_transcript_tracks(
    config: Config,
    resources: TrackResources,
) -> tuple[AnnotationTrack[Any], ...]:
    """Build the transcript consequence track a run always annotates with.

    Raises:
        ValueError: If the transcript files are ambiguous or absent.
        NotImplementedError: If the configured transcript provider has no reader yet.
    """
    source = build_transcript_source(
        resources.source_config,
        sequence_source=resources.sequence_source,
        synonyms=resources.synonyms,
    )
    lookup = TranscriptCandidateLookup(
        source=source,
        window=transcript_candidate_window(config),
        pool=TranscriptPool(max_bytes=config.execution.transcript_pool_bytes),
    )
    return (
        TranscriptConsequenceTrack(
            identifier=TRACK_IDENTIFIER,
            lookup=lookup,
            overlap_finder=build_transcript_overlap_finder(config, lookup.candidates),
            annotator=build_annotator(
                **asdict(config.annotation),
                upstream_distance=config.source.upstream_distance,
                downstream_distance=config.source.downstream_distance,
            ),
            batch_size=config.annotation.batch_size,
        ),
    )


class _CoverageIndex:
    """Merged genomic query coverage for one sequence region.

    Replaces an ``any(...)`` scan over every previously queried interval with a ``bisect`` over a
    merged, much shorter list. ``contains`` can return True for a query no single old interval
    covered on its own -- when the query lies inside a run of *strictly overlapping* merged
    intervals, it is still a hit. That is safe here: a transcript overlapping the query then
    overlaps one of the intervals in that run, and was already fetched and cached when that
    interval was added.

    Touching intervals (``previous.end + 1 == next.start``) are kept as separate entries rather
    than merged, since only strictly-overlapping merges (a shared coordinate) have been checked
    against this invariant.
    """

    __slots__ = ("_starts", "_ends")

    def __init__(self) -> None:
        self._starts: list[int] = []
        self._ends: list[int] = []

    def contains(self, interval: Interval) -> bool:
        """Return whether some already-added interval fully contains ``interval``."""
        query_start = interval.lookup_start
        query_end = interval.lookup_end
        index = bisect_right(self._starts, query_start) - 1
        return index >= 0 and self._ends[index] >= query_end

    def add(self, interval: Interval) -> None:
        """Merge a newly queried interval into the tracked set."""
        start = interval.lookup_start
        end = interval.lookup_end

        insert_at = bisect_left(self._starts, start)
        if insert_at > 0 and self._ends[insert_at - 1] >= start:
            insert_at -= 1
            start = self._starts[insert_at]
            end = max(end, self._ends[insert_at])

        # insert_at is either the existing entry just merged into start/end above, or a pure
        # insertion point with nothing there yet -- either way, scan forward from it (not past
        # it) for further entries the (possibly widened) interval now reaches.
        merge_end = insert_at
        while merge_end < len(self._starts) and self._starts[merge_end] <= end:
            end = max(end, self._ends[merge_end])
            merge_end += 1

        self._starts[insert_at:merge_end] = [start]
        self._ends[insert_at:merge_end] = [end]


def _coalesced_intervals(intervals: list[Interval]) -> list[Interval]:
    """Merge overlapping or directly adjacent source-query intervals."""
    ordered = sorted(intervals, key=lambda interval: (interval.lookup_start, interval.lookup_end))
    coalesced: list[Interval] = []
    for interval in ordered:
        if not coalesced or interval.lookup_start > coalesced[-1].lookup_end + 1:
            coalesced.append(interval)
            continue
        previous = coalesced[-1]
        coalesced[-1] = Interval(
            start=previous.lookup_start,
            end=max(previous.lookup_end, interval.lookup_end),
        )
    return coalesced


@define(kw_only=True, auto_attribs=True)
class TranscriptCandidateLookup:
    """Cached transcript lookup for an input buffer.

    Attributes:
        source (TranscriptSource): The transcript source every lookup reads.
        window (int): Bases added each side of an allele before querying.
        pool (TranscriptPool): Bounded cache of built transcripts, kept across payloads.
    """

    source: TranscriptSource
    window: int
    pool: TranscriptPool
    _cache: BinnedCache[Transcript] = field(factory=BinnedFeatureCache, init=False)
    _queried_regions: dict[str, _CoverageIndex] = field(factory=dict, init=False)
    _uses_complete_source_fallback: bool = field(default=False, init=False)

    def __attrs_post_init__(self) -> None:
        """Attach the pool to the source that builds transcripts through it.

        ``transcript_pool`` is declared on ``GXFTranscriptSource``, which both members of the
        ``TranscriptSource`` union subclass; naming the assignment keeps that requirement explicit
        rather than accidental via the union.
        """
        self.source.transcript_pool = self.pool

    def begin_payload(self) -> None:
        """Release transcript state retained by the previous payload.

        ``pool`` is deliberately not cleared: it is a bounded, per-worker cache of built
        transcripts meant to persist across payloads -- see its own docstring and the module
        docstring on ``caches/feature_pool.py``.
        """
        self._cache.clear()
        self._queried_regions.clear()

    def log_statistics(self) -> None:
        """Log the transcript pool's running counters, cheap enough to do every payload."""
        logger.debug(
            "transcript pool: hits=%d misses=%d evictions=%d region_evictions=%d poisoned=%s",
            self.pool.hit_count,
            self.pool.miss_count,
            self.pool.eviction_count,
            self.pool.seq_region_eviction_count,
            self.source.transcript_pool_poisoned,
        )

    def close(self) -> None:
        """Close the transcript source."""
        self.source.close()

    def candidates(self, input_buffer: VariationAlleleBuffer) -> list[Transcript]:
        """Return each transcript that can pair with this one-region buffer."""
        candidates: dict[int, Transcript] = {}
        intervals_by_region: dict[str, list[Interval]] = {}
        for allele in input_buffer.objects:
            region = allele.slice.seq_region_name if allele.slice is not None else None
            if region is None:
                continue
            intervals_by_region.setdefault(region, []).append(
                Interval(
                    start=max(1, min(allele.start, allele.end) - self.window),
                    end=max(allele.start, allele.end) + self.window,
                )
            )

        source_intervals = {
            region: _coalesced_intervals(intervals) for region, intervals in intervals_by_region.items()
        }
        if self._uses_complete_source_fallback:
            self._load_complete_source(source_intervals)
        else:
            for region, intervals in source_intervals.items():
                if self._load_region(region, intervals):
                    continue
                self._cache.clear()
                self._queried_regions.clear()
                self._load_complete_source(source_intervals)
                break

        for region, intervals in intervals_by_region.items():
            for transcript in self._cache.get_objects(region, intervals):
                candidates[id(transcript)] = transcript
        return list(candidates.values())

    def _load_region(self, sequence_region: str, intervals: list[Interval]) -> bool:
        """Load indexed transcript graphs and report whether regional access succeeded."""
        if self._uses_complete_source_fallback:
            return False
        coverage = self._queried_regions.setdefault(sequence_region, _CoverageIndex())
        for interval in intervals:
            if coverage.contains(interval):
                continue
            try:
                transcripts = self.source.get_features_by_region(
                    seq_region_name=sequence_region,
                    start=interval.lookup_start,
                    end=interval.lookup_end,
                )
            except OSError:
                logger.warning(
                    "%s has no usable region index; scanning each payload by gene",
                    self.source.file_path,
                )
                self._uses_complete_source_fallback = True
                return False
            for transcript in transcripts:
                self._cache.add(
                    sequence_region,
                    transcript,
                    Interval(start=transcript.start, end=transcript.end),
                    key=(sequence_region, transcript.stable_id, transcript.start, transcript.end),
                )
            coverage.add(interval)
        return True

    def _load_complete_source(self, intervals_by_region: dict[str, list[Interval]]) -> None:
        """Stream an unindexed GXF and retain only features needed by this payload."""
        uncovered: dict[str, list[Interval]] = {}
        for sequence_region, intervals in intervals_by_region.items():
            coverage = self._queried_regions.setdefault(sequence_region, _CoverageIndex())
            for interval in intervals:
                if coverage.contains(interval):
                    continue
                uncovered.setdefault(sequence_region, []).append(interval)
        if not uncovered:
            return

        source_intervals = {
            sequence_region: [(interval.lookup_start, interval.lookup_end) for interval in intervals]
            for sequence_region, intervals in uncovered.items()
        }
        for transcript in self.source.iter_features_by_regions(source_intervals):
            sequence_region = (
                transcript.slice.seq_region_name if transcript.slice is not None else UNKNOWN_VALUE
            )
            intervals = uncovered.get(sequence_region, [])
            if not any(
                transcript.start <= interval.lookup_end and transcript.end >= interval.lookup_start
                for interval in intervals
            ):
                continue
            self._cache.add(
                sequence_region,
                transcript,
                Interval(start=transcript.start, end=transcript.end),
                key=(sequence_region, transcript.stable_id, transcript.start, transcript.end),
            )
        for sequence_region, intervals in uncovered.items():
            coverage = self._queried_regions[sequence_region]
            for interval in intervals:
                coverage.add(interval)
        self._uses_complete_source_fallback = True


@define(kw_only=True, auto_attribs=True)
class TranscriptConsequenceTrack(OverlapAnnotationTrack[OverlapConsequence]):
    """Transcript consequence annotation track.

    Attributes:
        lookup (TranscriptCandidateLookup): The source and caches this track's finder reads through.
    """

    lookup: TranscriptCandidateLookup

    def begin_payload(self, input_buffer: VariationAlleleBuffer) -> None:
        """Release the previous payload's transcripts and fetch this payload's in one pass.

        Priming with the whole payload is what makes each annotation window a cache hit rather
        than another source query.
        """
        self.lookup.begin_payload()
        self.lookup.log_statistics()
        self.lookup.candidates(input_buffer)

    def attach(
        self,
        variations: Sequence[Variation],
        results: Sequence[TrackResult[OverlapConsequence]],
    ) -> None:
        """Attach consequences to the transcript overlaps of each allele.

        Raises:
            TypeError: If a result holds an overlap this track cannot have produced.
        """
        by_allele: dict[int, list[TranscriptVariationAllele]] = {}
        for result in results:
            overlap = result.overlap
            if not isinstance(overlap, TranscriptVariationAllele):
                raise TypeError(
                    f"Track '{self.identifier}' cannot attach a {type(overlap).__name__} overlap."
                )
            overlap.overlap_consequences = list(result.annotations)
            by_allele.setdefault(id(overlap.overlapping_feature), []).append(overlap)
        for variation in variations:
            for allele in variation.variation_alleles:
                allele.transcript_variation_alleles = by_allele.get(id(allele), [])

    def close(self) -> None:
        """Close the transcript source this track reads."""
        self.lookup.close()
