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
"""Worker-local runtime for the default annotation runner."""

from __future__ import annotations

import logging
from bisect import bisect_left, bisect_right
from itertools import batched
from pathlib import Path
from typing import Any

from attrs import frozen

from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.buffer import InputBuffer
from ensembl.vep.buffer.reader import input_header_lines
from ensembl.vep.buffer.records import VCFBufferPayload
from ensembl.vep.caches import BinnedFeatureCache, Interval, TranscriptPool
from ensembl.vep.config import Config
from ensembl.vep.matcher import OverlapMatcher
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import Variation, VariationAllele
from ensembl.vep.output.factory.base import AnnotatedVariant
from ensembl.vep.output.factory.write.base import OutputFactory
from ensembl.vep.parsers.vcf import VCFParser
from ensembl.vep.sources.build import (
    TranscriptSource,
    build_annotation_sources,
    build_chromosome_synonyms,
    build_input_source,
    build_sequence_source,
    transcript_source_of,
)
from ensembl.vep.stats import RunSummaryCollector, RunSummaryDelta
from ensembl.vep.workflow import compile_workflow

__all__ = ["AnnotatedBufferResult", "BufferedWorker", "annotate_buffer", "build_buffered_worker"]

logger = logging.getLogger(__name__)

_VCF_TRANSPORT_HEADER = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


@frozen(kw_only=True, auto_attribs=True)
class AnnotatedBufferResult:
    """Carry one ordered worker result without retaining model objects.

    Attributes:
        index (int): Input payload index.
        variants (tuple[AnnotatedVariant, ...]): Detached annotated variants.
        summary_delta (RunSummaryDelta): Mergeable run-summary counts.
    """

    index: int
    variants: tuple[AnnotatedVariant, ...]
    summary_delta: RunSummaryDelta


def _read_input_headers(path: Path) -> list[str]:
    """Return input headers, or an empty list when the source cannot be read yet."""
    try:
        return input_header_lines(path)
    except (OSError, UnicodeDecodeError):
        return []


def _attach_transcript_pool(source: TranscriptSource, pool: TranscriptPool) -> None:
    """Attach a per-worker transcript pool to a transcript source.

    Both members of the ``TranscriptSource`` union are ``GXFTranscriptSource`` subclasses, which
    is where ``transcript_pool`` is declared; naming this assignment keeps that requirement
    explicit at the one call site that relies on it, rather than accidental via the union.
    """
    source.transcript_pool = pool


class _CoverageIndex:
    """Track already-queried genomic ranges for one sequence region, merged and sorted.

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


class BufferedWorker:
    """Annotate primitive VCF buffers with worker-local sources and workflow state."""

    def __init__(self, config: Config) -> None:
        """Build sources, caches, and workflow after the worker starts."""
        source_config = config.source
        assert source_config.input_file is not None
        synonyms = build_chromosome_synonyms(source_config)
        self.sequence_source = build_sequence_source(source_config)
        self.annotation_sources = build_annotation_sources(
            source_config, sequence_source=self.sequence_source, synonyms=synonyms
        )
        self.variation_source = build_input_source(
            source_config, sequence_source=self.sequence_source, synonyms=synonyms
        )
        headers = _read_input_headers(Path(source_config.input_file))
        self._vcf_transport_header = "\n".join(headers) + "\n" if headers else _VCF_TRANSPORT_HEADER
        self.config = config
        self.transcript_cache: BinnedFeatureCache[Any] = BinnedFeatureCache()
        self.transcript_source = transcript_source_of(self.annotation_sources)
        # Persists across payloads, unlike transcript_cache above: a bounded, per-worker cache of
        # built Transcript models so a gene revisited by a later payload is not rebuilt from GXF
        # records. Never cleared by _begin_payload.
        self.transcript_pool = TranscriptPool(max_bytes=config.execution.transcript_pool_bytes)
        _attach_transcript_pool(self.transcript_source, self.transcript_pool)
        self._queried_regions: dict[str, _CoverageIndex] = {}
        self._uses_complete_source_fallback = False

        annotation = config.annotation
        consequence_annotator = build_annotator(
            annotator=annotation.annotator,
            batch_size=annotation.batch_size,
            early_subsumption_skip=annotation.early_subsumption_skip,
            zone_map=annotation.zone_map,
            zone_map_bytes=annotation.zone_map_bytes,
            codon_cache=annotation.codon_cache,
            codon_cache_entries=annotation.codon_cache_entries,
            mapper_cache=annotation.mapper_cache,
            mapper_cache_bytes=annotation.mapper_cache_bytes,
            readthrough=annotation.readthrough,
            fact_cost_override=annotation.fact_cost_override,
            upstream_distance=self.config.source.upstream_distance,
            downstream_distance=self.config.source.downstream_distance,
        )
        worker_context = {
            "transcript_candidates": self._transcript_candidates,
            "transcript_matcher": OverlapMatcher(),
            "consequence_annotator": consequence_annotator,
        }
        self.workflow = compile_workflow(config, worker_context)
        self.snapshot_factory = OutputFactory()

    def _transcript_candidates(self, input_buffer: InputBuffer) -> list[Any]:
        """Return each transcript that can pair with this one-region buffer."""
        candidates: dict[int, Any] = {}
        intervals_by_region: dict[str, list[Interval]] = {}
        for allele in input_buffer.objects:
            region = allele.slice.seq_region_name if allele.slice is not None else None
            if region is None:
                continue
            # The wider of the two distances: which one applies depends on each transcript's
            # strand, and the overlap finder is what decides that once they are fetched.
            window = max(self.config.source.upstream_distance, self.config.source.downstream_distance)
            intervals_by_region.setdefault(region, []).append(
                Interval(
                    start=max(1, min(allele.start, allele.end) - window),
                    end=max(allele.start, allele.end) + window,
                )
            )

        source_intervals = {
            region: self._coalesced_intervals(intervals) for region, intervals in intervals_by_region.items()
        }
        if self._uses_complete_source_fallback:
            self._load_complete_transcript_source(source_intervals)
        else:
            for region, intervals in source_intervals.items():
                if self._load_transcripts(region, intervals):
                    continue
                self.transcript_cache.clear()
                self._queried_regions.clear()
                self._load_complete_transcript_source(source_intervals)
                break

        for region, intervals in intervals_by_region.items():
            for transcript in self.transcript_cache.get_objects(region, intervals):
                candidates[id(transcript)] = transcript
        return list(candidates.values())

    @staticmethod
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

    def _load_transcripts(self, sequence_region: str, intervals: list[Interval]) -> bool:
        """Load indexed transcript graphs and report whether regional access succeeded."""
        if self._uses_complete_source_fallback:
            return False
        coverage = self._queried_regions.setdefault(sequence_region, _CoverageIndex())
        for interval in intervals:
            if coverage.contains(interval):
                continue
            try:
                transcripts = self.transcript_source.get_features_by_region(
                    seq_region_name=sequence_region,
                    start=interval.lookup_start,
                    end=interval.lookup_end,
                )
            except OSError:
                logger.warning(
                    "%s has no usable region index; scanning each payload by gene",
                    self.transcript_source.file_path,
                )
                self._uses_complete_source_fallback = True
                return False
            for transcript in transcripts:
                self.transcript_cache.add(
                    sequence_region,
                    transcript,
                    Interval(start=transcript.start, end=transcript.end),
                    key=(sequence_region, transcript.stable_id, transcript.start, transcript.end),
                )
            coverage.add(interval)
        return True

    def _load_complete_transcript_source(
        self,
        intervals_by_region: dict[str, list[Interval]],
    ) -> None:
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
        for transcript in self.transcript_source.iter_features_by_regions(source_intervals):
            sequence_region = transcript.slice.seq_region_name if transcript.slice is not None else "?"
            intervals = uncovered.get(sequence_region, [])
            if not any(
                transcript.start <= interval.lookup_end and transcript.end >= interval.lookup_start
                for interval in intervals
            ):
                continue
            self.transcript_cache.add(
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

    def _begin_payload(self) -> None:
        """Release transcript state retained by the previous process payload.

        ``transcript_source.transcript_pool`` is deliberately not cleared here: it is a bounded,
        per-worker cache of built transcripts meant to persist across payloads -- see its own
        docstring and the module docstring on ``caches/feature_pool.py``.
        """
        self.transcript_cache.clear()
        self._queried_regions.clear()

    def annotate(self, payload: VCFBufferPayload) -> AnnotatedBufferResult:
        """Return detached annotations for one primitive record buffer."""
        self._begin_payload()
        records = VCFParser.parse_text(
            self._vcf_transport_header + "\n".join("\t".join(row) for row in payload.columns)
        )
        variations = self.variation_source.records_to_features(records)
        unknown_regions = frozenset(self.sequence_source.unknown_regions)
        alleles = self._annotatable_alleles(variations, unknown_regions)
        annotated_variants: list[AnnotatedVariant] = []
        collector = RunSummaryCollector()

        # Prime source/cache lookup once. Windows retain only their temporary overlap graph.
        self._transcript_candidates(InputBuffer(objects=tuple(alleles)))
        variation_columns = zip(variations, payload.columns, strict=True)
        for group in batched(variation_columns, self.config.execution.annotation_window):
            group_variations = [variation for variation, _columns in group]
            group_alleles = self._annotatable_alleles(group_variations, unknown_regions)
            results = []
            if group_alleles:
                results = self.workflow.annotate(
                    InputBuffer(objects=tuple(group_alleles)),
                    threads=self.config.execution.threads,
                )
            self._attach_results(group_variations, results)
            for variation, columns in group:
                region = variation.slice.seq_region_name if variation.slice is not None else None
                annotated_variants.append(
                    self.snapshot_factory.build_annotated_variant(
                        variation,
                        columns,
                        skipped=region in unknown_regions,
                    )
                )
                collector.observe(variation, sequence_region_name=columns[0])
            self._clear_transient_annotations(group_variations)
        self._log_transcript_pool_stats()
        return AnnotatedBufferResult(
            index=payload.index,
            variants=tuple(annotated_variants),
            summary_delta=collector.delta(),
        )

    def _log_transcript_pool_stats(self) -> None:
        """Log the transcript pool's running counters, cheap enough to do every payload."""
        pool = self.transcript_pool
        logger.debug(
            "transcript pool: hits=%d misses=%d evictions=%d region_evictions=%d poisoned=%s",
            pool.hit_count,
            pool.miss_count,
            pool.eviction_count,
            pool.region_eviction_count,
            self.transcript_source.transcript_pool_poisoned,
        )

    @staticmethod
    def _annotatable_alleles(
        variations: list[Variation],
        unknown_regions: frozenset[str],
    ) -> list[VariationAllele]:
        """Return alleles whose sequence region is available to annotation sources."""
        alleles: list[VariationAllele] = []
        for variation in variations:
            region = variation.slice.seq_region_name if variation.slice is not None else None
            if region not in unknown_regions:
                alleles.extend(variation.variation_alleles)
        return alleles

    @staticmethod
    def _attach_results(variations: list[Variation], results: list[Any]) -> None:
        """Attach one annotation window's track results to its variation alleles."""
        by_allele: dict[int, list[TranscriptVariationAllele]] = {}
        for result in results:
            overlap = result.overlap
            if not isinstance(overlap, TranscriptVariationAllele):
                continue
            overlap.overlap_consequences = list(result.annotations)
            by_allele.setdefault(id(overlap.overlapping_feature), []).append(overlap)
        for variation in variations:
            for allele in variation.variation_alleles:
                allele.transcript_variation_alleles = by_allele.get(id(allele), [])

    @staticmethod
    def _clear_transient_annotations(variations: list[Variation]) -> None:
        """Release a rendered window's overlap graph before annotating the next window."""
        for variation in variations:
            for allele in variation.variation_alleles:
                allele.transcript_variation_alleles = []


def build_buffered_worker(config: Config) -> BufferedWorker:
    """Build one worker-local buffered annotation runtime."""
    return BufferedWorker(config)


def annotate_buffer(
    worker: BufferedWorker,
    payload: VCFBufferPayload,
) -> AnnotatedBufferResult:
    """Annotate one primitive buffer through a worker runtime."""
    return worker.annotate(payload)
