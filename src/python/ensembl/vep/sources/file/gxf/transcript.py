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
"""Shared GXF-based transcript source.

This module contains the common transcript-building workflow used by both file
formats:

* buffer records from the parser in streaming order,
* group the buffered records into transcript-sized blocks,
* build Transcript, Gene, Exon, and Translation objects (introns are derived lazily by
  ``Transcript`` itself), and
* reuse Slice and Gene instances across records that share the same region.

Only the parser-specific record identity and parent resolution rules are left
to the concrete GFF3 and GTF adapters.
"""

from __future__ import annotations

import logging
import os
from abc import ABC, abstractmethod
from collections import defaultdict
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from typing import ClassVar, Generic, Protocol, TypeVar, cast

from attrs import define, evolve, field

from ensembl.vep.caches.feature_pool import FeaturePool
from ensembl.vep.common import (
    CoordSystemLevel,
    ExonPhase,
    STRAND_SYMBOL_TO_ENUM,
    GFF_TO_ENSEMBL_PHASE,
    Strand,
    format_region_string,
    parse_region_string,
)
from ensembl.vep.models.genic import Exon, Gene, Transcript, Translation
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.parsers import Parser
from ensembl.vep.sources.file.file_feature_source import FileFeatureSource
from ensembl.vep.sources.synonym import SequenceRegionResolver

logger = logging.getLogger(__name__)

_CHROMOSOME_COORD_SYSTEM = CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME, default=True)
_PLACEHOLDER_SEQ_REGION_END = 2**31 - 1
_ORDERED_FILE_ASSUMPTION_MESSAGE = (
    "GFF3/GTF streaming iteration assumes records are grouped by gene. A gene record may precede "
    "several transcript records and their exon/CDS children. Files that revisit a gene after "
    "starting another gene require an index for complete assembly."
)

P = TypeVar("P", bound=Parser)


class _GXFRecord(Protocol):
    """Minimal parser record interface shared by GFF3 and GTF records.

    Attributes:
        feature_type (str): Raw feature type reported by the parser.
        seqid (str): Sequence-region name for the record.
        start (int): Genomic start position, 1-based inclusive.
        end (int): Genomic end position, 1-based inclusive.
        strand_symbol (str): Strand symbol from the parser record.
        phase_value (str | None): Phase string for coding features, if present.
        attributes (dict[str, tuple[str, ...]]): Parsed key/value attributes from the record.
    """

    feature_type: str
    seqid: str
    start: int
    end: int
    strand_symbol: str
    phase_value: int | None
    attributes: dict[str, tuple[str, ...]]

    def first_attribute(self, key: str) -> str | None:
        """Return the first attribute value for a key, or None when absent."""


@define(kw_only=True, auto_attribs=True)
class _GXFTranscriptGroup:
    """Grouped parser records describing one transcript.

    ``transcript_id`` and ``gene_id`` are resolved once, while grouping, from the same
    ``get_record_id``/``get_parent_id`` calls that route records into this group in the first
    place. Carrying them here is what lets the builder stop recomputing them per transcript.

    Attributes:
        transcript_id (str): Stable identifier of the transcript record.
        gene_id (str | None): Parent gene identifier declared by the transcript record, if any --
            not necessarily the id of ``gene_record``, which may be absent even when this is set.
        transcript_record (_GXFRecord): Primary transcript record for the group.
        gene_record (_GXFRecord | None): Associated gene record, if one was present in the input.
        exon_records (list[_GXFRecord]): Exon records attached to the transcript.
        cds_records (list[_GXFRecord]): Coding records attached to the transcript.
    """

    transcript_id: str
    gene_id: str | None
    transcript_record: _GXFRecord
    gene_record: _GXFRecord | None
    exon_records: list[_GXFRecord]
    cds_records: list[_GXFRecord]


@define(kw_only=True, auto_attribs=True)
class GXFTranscriptSource(FileFeatureSource[P, Transcript], ABC, Generic[P]):
    """Transcript source backed by GFF3 or GTF file.

    Attributes:
        gene_feature_types (frozenset[str]): Lower-cased feature types that identify a gene record.
        transcript_feature_types (frozenset[str]): Lower-cased feature types that identify a
            transcript record.
        exon_feature_types (frozenset[str]): Lower-cased feature types that identify an exon record.
        coding_feature_types (frozenset[str]): Lower-cased feature types that contribute to the
            transcript translation.
        reference_region_resolver (SequenceRegionResolver | None): Resolve raw GXF names to the
            shared reference names used by model slices.
        source_region_resolver (SequenceRegionResolver | None): Resolve model slice names to names
            accepted by indexed GXF queries.
        transcript_pool (FeaturePool[Transcript] | None): Cross-payload cache of built transcripts,
            assigned by the owning worker; ``None`` disables reuse.
    """

    gene_feature_types: ClassVar[frozenset[str]] = frozenset({"gene"})
    transcript_feature_types: ClassVar[frozenset[str]] = frozenset({"transcript"})
    exon_feature_types: ClassVar[frozenset[str]] = frozenset({"exon"})
    coding_feature_types: ClassVar[frozenset[str]] = frozenset({"cds"})
    adaptor: P
    reference_region_resolver: SequenceRegionResolver | None = field(default=None)
    source_region_resolver: SequenceRegionResolver | None = field(default=None)
    transcript_pool: FeaturePool[Transcript] | None = field(default=None)
    _query_parser: P | None = field(default=None, init=False, repr=False)
    _query_parser_pid: int | None = field(default=None, init=False, repr=False)
    _pool_poisoned: bool = field(default=False, init=False, repr=False)

    @property
    def transcript_pool_poisoned(self) -> bool:
        """Return whether a repeated transcript id has disabled pooling for this source."""
        return self._pool_poisoned

    @property
    def sequence_region_names(self) -> frozenset[str]:
        """Return the tabix index's reference names together with ``##sequence-region`` names.

        A header can declare a contig that has no rows, and so no index entry.
        """
        try:
            with self.adaptor.open(self.file_path, mode="auto") as adaptor:
                index_names = frozenset(adaptor.index_sequence_names)
                headers = adaptor.metadata.get("header", ())
        except OSError:
            return frozenset()

        header_names: set[str] = set()
        if isinstance(headers, (list, tuple)):
            for header in headers:
                fields = str(header).split()
                if len(fields) >= 4 and fields[0] == "##sequence-region":
                    header_names.add(fields[1])
        return index_names | header_names

    def iter_all_features(self) -> Iterator[Transcript]:
        """Yield transcripts by streaming the file sequentially.

        Returns:
            Iterator[Transcript]: Transcripts built from the entire file.
        """

        logger.info(_ORDERED_FILE_ASSUMPTION_MESSAGE)

        with self.adaptor.open(self.file_path, mode="sequential") as adaptor:
            yield from self._iter_features_from_adaptor(cast(P, adaptor))

    def iter_features_by_regions(
        self,
        intervals_by_region: Mapping[str, Sequence[tuple[int, int]]],
    ) -> Iterator[Transcript]:
        """Yield only transcripts intersecting reference-space query intervals.

        This source's own method, not part of the ``FeatureSource`` contract: a format whose record
        already carries a whole feature answers each region on its own, and gains nothing from the
        batch. Here one pass serves every requested interval, so a whole input buffer costs one
        read of the file rather than one per region.

        The record predicate runs before transcript, exon, CDS, and gene models are built, which
        keeps the unindexed fallback bounded by the current payload rather than the source size.

        Args:
            intervals_by_region (Mapping[str, Sequence[tuple[int, int]]]): 1-based inclusive
                ``(start, end)`` intervals, keyed by sequence-region name.

        Returns:
            Iterator[Transcript]: Transcripts intersecting any requested interval.
        """

        logger.info(_ORDERED_FILE_ASSUMPTION_MESSAGE)

        def overlaps_requested_interval(record: _GXFRecord) -> bool:
            sequence_region = record.seqid
            if self.reference_region_resolver is not None:
                sequence_region = self.reference_region_resolver.resolve(sequence_region)
            intervals = intervals_by_region.get(sequence_region, ())
            return any(record.start <= end and record.end >= start for start, end in intervals)

        with self.adaptor.open(self.file_path, mode="sequential") as adaptor:
            yield from self._iter_features_from_adaptor(
                cast(P, adaptor),
                transcript_predicate=overlaps_requested_interval,
            )

    def iter_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> Iterator[Transcript]:
        """Yield complete transcripts intersecting one query region.

        One region is the batch of one: an indexed query cannot serve this, because a transcript's
        exons reach outside the queried region, so the answer is the same sequential pass
        ``iter_features_by_regions`` makes. Callers with an index want ``get_features_by_region``,
        which queries each candidate's full span instead.

        Args:
            region (str | None): The region (e.g., "chr1:1000-2000") to get transcripts from.
            seq_region_name (str | None): The name of the sequence region (e.g., "chr1").
            start (int | None): The start position of the region.
            end (int | None): The end position of the region.

        Returns:
            Iterator[Transcript]: Transcripts intersecting the region.
        """

        try:
            query_region = format_region_string(region, seq_region_name=seq_region_name, start=start, end=end)
            region_name, region_start, region_end = parse_region_string(query_region)
        except ValueError:
            logger.debug("Cannot parse region: %s %s:%s-%s", region, seq_region_name, start, end)
            return

        interval = (region_start or 1, region_end or _PLACEHOLDER_SEQ_REGION_END)
        yield from self.iter_features_by_regions({region_name: [interval]})

    def get_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> list[Transcript]:
        """Return complete transcripts whose spans meet one query region.

        Two queries, because one is not enough: stage 1 finds the transcript rows overlapping the
        region, stage 2 queries the full span of each candidate not already pooled, so exons and
        coding records outside the region are kept. Results are in stage-1 (file) order whether a
        transcript came from the pool or was just built, which is what fixes CSQ order.

        A region without coordinates queries the whole sequence region, and so assembles every
        transcript on it.
        """
        try:
            query_region = format_region_string(region, seq_region_name=seq_region_name, start=start, end=end)
            region_name, region_start, region_end = parse_region_string(query_region)
        except ValueError:
            logger.debug("Cannot parse region: %s %s:%s-%s", region, seq_region_name, start, end)
            return []

        records = self._records_for_region(
            region_name, region_start, region_end, feature_types=self.transcript_feature_types
        )
        candidates: dict[str, tuple[int, int]] = {}
        has_duplicate_candidate = False
        for record in records:
            transcript_id = self.get_record_id(record)
            if transcript_id is None:
                continue
            previous_span = candidates.get(transcript_id)
            span = (record.start, record.end)
            if previous_span is not None and previous_span != span:
                has_duplicate_candidate = True
                self._disable_pool_for_duplicate_id(transcript_id)
            candidates[transcript_id] = span
        if not candidates:
            return []

        pool = None if self._pool_poisoned else self.transcript_pool
        if pool is None or pool.max_bytes <= 0 or has_duplicate_candidate:
            return self._build_candidates(region_name, candidates)
        return self._build_candidates_with_pool(region_name, candidates, pool)

    def _build_candidates(
        self,
        seq_region_name: str,
        candidates: Mapping[str, tuple[int, int]],
    ) -> list[Transcript]:
        """Query every candidate's span and build all of them, without the pool."""

        complete_records: list[_GXFRecord] = []
        for span_start, span_end in self._merged_spans(candidates.values()):
            complete_records.extend(self._records_for_region(seq_region_name, span_start, span_end))
        return self._build_transcripts(
            complete_records, gene_records={}, gene_cache={}, candidate_ids=candidates
        )

    def _build_candidates_with_pool(
        self,
        seq_region_name: str,
        candidates: Mapping[str, tuple[int, int]],
        pool: FeaturePool[Transcript],
    ) -> list[Transcript]:
        """Reuse pooled transcripts, query and build only the missing ones, in candidate order.

        A pool hit skips its own span, so stage-2 rows are also checked for a pooled id at a
        different span.
        """

        pooled: dict[str, Transcript] = {}
        missing: dict[str, tuple[int, int]] = {}
        for transcript_id, span in candidates.items():
            start, end = span
            cached = pool.get((seq_region_name, transcript_id, start, end))
            if cached is None:
                missing[transcript_id] = span
            else:
                pooled[transcript_id] = cached

        if not missing:
            return [pooled[transcript_id] for transcript_id in candidates if transcript_id in pooled]

        complete_records: list[_GXFRecord] = []
        for span_start, span_end in self._merged_spans(missing.values()):
            complete_records.extend(self._records_for_region(seq_region_name, span_start, span_end))

        resurfaced_id = self._first_pooled_id_in_records(complete_records, pooled)
        if resurfaced_id is not None:
            self._disable_pool_for_duplicate_id(resurfaced_id)
            return self._build_candidates(seq_region_name, candidates)

        built = self._build_transcripts(
            complete_records, gene_records={}, gene_cache={}, candidate_ids=missing
        )
        built_by_id: dict[str, Transcript] = {}
        for transcript in built:
            built_by_id[transcript.stable_id] = transcript
            pool.put(
                (seq_region_name, transcript.stable_id, transcript.start, transcript.end),
                transcript,
            )

        combined = {**pooled, **built_by_id}
        return [combined[transcript_id] for transcript_id in candidates if transcript_id in combined]

    def _first_pooled_id_in_records(
        self, records: Iterable[_GXFRecord], pooled: Mapping[str, Transcript]
    ) -> str | None:
        """Return the first transcript id in ``records`` whose span differs from its pooled model."""
        for record in records:
            if record.feature_type.lower() not in self.transcript_feature_types:
                continue
            record_id = self.get_record_id(record)
            if record_id is None:
                continue
            cached = pooled.get(record_id)
            if cached is not None and (record.start != cached.start or record.end != cached.end):
                return record_id
        return None

    def _disable_pool_for_duplicate_id(self, transcript_id: str) -> None:
        """Clear and permanently disable the pool once a transcript id names two rows."""
        if self._pool_poisoned:
            return
        self._pool_poisoned = True
        if self.transcript_pool is not None:
            self.transcript_pool.clear()
        logger.warning(
            "%s has a repeated transcript id %r; disabling the transcript pool for this source",
            self.file_path,
            transcript_id,
        )

    def _records_for_region(
        self,
        seq_region_name: str,
        start: int | None = None,
        end: int | None = None,
        *,
        feature_types: frozenset[str] | None = None,
    ) -> list[_GXFRecord]:
        """Return raw records from an indexed query, without assembling incomplete models."""
        source_name = seq_region_name
        if self.source_region_resolver is not None:
            source_name = self.source_region_resolver.resolve(seq_region_name)
        records = self._query_records(source_name, start, end, feature_types=feature_types)
        source_name_is_declared = (
            self.source_region_resolver is not None
            and source_name in self.source_region_resolver.source_names
        )
        if records or source_name_is_declared:
            return records

        # If no records are found using the source name declared in the index file or header,
        # try common fallback names (chr prefix variations). This can happen for an un-indexed
        # GFF/GTF file that doesn't have the ##sequence-region directive.
        fallback_names = [seq_region_name]
        fallback_names.append(
            seq_region_name[3:] if seq_region_name[:3].lower() == "chr" else f"chr{seq_region_name}"
        )
        for fallback_name in fallback_names:
            if fallback_name == source_name:
                continue
            records = self._query_records(fallback_name, start, end, feature_types=feature_types)
            if records:
                return records
        return []

    def _query_records(
        self,
        seq_region_name: str,
        start: int | None = None,
        end: int | None = None,
        *,
        feature_types: frozenset[str] | None = None,
    ) -> list[_GXFRecord]:
        """Read raw records using one source-specific sequence-region name.

        Args:
            seq_region_name: The source-spelling sequence-region name to query.
            start: Inclusive query start.
            end: Inclusive query end.
            feature_types: When given, only rows whose column 3 matches one of these values
                (case-insensitively) are returned. Applied on the Rust side before a non-matching
                line is parsed into a record.
        """
        region = format_region_string(seq_region_name=seq_region_name, start=start, end=end)
        return self._query_typed_records(region, feature_types=feature_types)

    @abstractmethod
    def _query_typed_records(
        self,
        region: str,
        *,
        feature_types: frozenset[str] | None = None,
    ) -> list[_GXFRecord]:
        """Seek the persistent query parser to a region and read its records.

        Implemented per format because only the GFF3 and GTF parsers accept ``feature_types``;
        ``Parser``, which bounds ``P``, does not. A subclass binds ``P`` to its own parser type, so
        the call type-checks there without narrowing ``self.adaptor`` here.

        Returns:
            list[_GXFRecord]: The records the query returned, empty when it matched nothing.
        """

        raise NotImplementedError

    def close(self) -> None:
        """Close the region-query parser and the attached sequence source."""

        if self._query_parser is not None:
            self._query_parser.close()
            self._query_parser = None
            self._query_parser_pid = None
        super().close()

    def _open_query_parser(self) -> P:
        """Return a query-mode parser kept open for this process's region queries.

        Annotation issues thousands of region queries; re-opening the file and its index for each
        dominated the run time. The parser is separate from ``adaptor``, which other methods open
        and close in other modes, and a process reached through fork opens its own so it never
        shares a file offset with its parent.
        """
        if self._query_parser is None or self._query_parser_pid != os.getpid():
            query_parser = evolve(self.adaptor)
            query_parser.open(self.file_path, mode="query")
            self._query_parser = query_parser
            self._query_parser_pid = os.getpid()
        return self._query_parser

    @staticmethod
    def _merged_spans(spans: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
        """Return overlapping transcript spans as minimal query ranges."""
        ordered = sorted(spans)
        merged: list[tuple[int, int]] = []
        for start, end in ordered:
            if not merged or start > merged[-1][1] + 1:
                merged.append((start, end))
                continue
            previous_start, previous_end = merged[-1]
            merged[-1] = (previous_start, max(previous_end, end))
        return merged

    def _records_to_features(self, records: list) -> list[Transcript]:
        """Convert a batch of records into transcripts, with a fresh gene and slice context.

        The abstract conversion ``FileFeatureSource`` requires, and what its ``get_all_features``
        calls. Fresh dictionaries are passed every time: ``_build_transcripts`` defaults to shared
        mutable dictionaries, which would carry genes and slices over from one call to the next.
        """
        return self._build_transcripts(records, gene_records={}, gene_cache={})

    def _iter_features_from_adaptor(
        self,
        adaptor: P,
        *,
        transcript_predicate: Callable[[_GXFRecord], bool] | None = None,
    ) -> Iterator[Transcript]:
        """Yield transcripts from an open parser, buffering records into blocks.

        Buffers one contiguous gene block from the parser, then converts its grouped records into
        Transcript objects.

        Args:
            adaptor: An open parser positioned at the records to stream.

        Returns:
            Iterator[Transcript]: Transcripts built from buffered GFF3 records.
        """

        gene_cache: dict[str, Gene] = {}
        buffer: list[_GXFRecord] = []
        current_gene_key: str | None = None

        while (record := adaptor.next()) is not None:
            gene_key = self._stream_gene_key(record)
            if gene_key is not None and current_gene_key not in (None, gene_key):
                yield from self._stream_buffer_features(
                    buffer,
                    gene_cache=gene_cache,
                    transcript_predicate=transcript_predicate,
                )
                buffer = []
                gene_cache.clear()
            if gene_key is not None:
                current_gene_key = gene_key
            buffer.append(record)

        if buffer:
            yield from self._stream_buffer_features(
                buffer,
                gene_cache=gene_cache,
                transcript_predicate=transcript_predicate,
            )

    def _stream_buffer_features(
        self,
        records: list[_GXFRecord],
        *,
        gene_cache: dict[str, Gene],
        transcript_predicate: Callable[[_GXFRecord], bool] | None,
    ) -> list[Transcript]:
        """Build requested transcripts from one gene-sized streamed record block."""

        groups = self._group_records(records)
        if transcript_predicate is not None:
            groups = [group for group in groups if transcript_predicate(group.transcript_record)]
        return self._groups_to_features_with_context(
            groups,
            gene_cache=gene_cache,
        )

    def _stream_gene_key(self, record: _GXFRecord) -> str | None:
        """Return the gene-sized streaming block key for one ordered record."""

        explicit_gene_id = self._first_attr(record, "gene_id")
        if explicit_gene_id is not None:
            return explicit_gene_id
        feature_type = record.feature_type.lower()
        if feature_type in self.gene_feature_types:
            return self.get_record_id(record)
        if feature_type in self.transcript_feature_types:
            parent_ids = self.get_parent_ids(record)
            if parent_ids:
                return parent_ids[0]
            return f"transcript:{self.get_record_id(record)}"
        return None

    @abstractmethod
    def get_record_id(self, record: _GXFRecord) -> str | None:
        """Return the stable identifier for a GFF3/GTF record.

        Args:
            record: The record to inspect.

        Returns:
            str | None: The record identifier, or None when the record is not
            identified.
        """

        raise NotImplementedError

    @abstractmethod
    def get_parent_ids(self, record: _GXFRecord) -> tuple[str, ...]:
        """Return the parent identifiers for a GFF3/GTF record, if they exist.

        Args:
            record: The record to inspect.

        Returns:
            tuple[str, ...]: A tuple of parent identifiers, or an empty tuple when the record has no
            parent.
        """

        raise NotImplementedError

    def _build_transcripts(
        self,
        records: list[_GXFRecord],
        *,
        gene_records: dict[str, _GXFRecord] | None = None,
        gene_cache: dict[str, Gene] | None = None,
        candidate_ids: Mapping[str, tuple[int, int]] | None = None,
    ) -> list[Transcript]:
        """Convert buffered records into Transcript objects while reusing caches.

        Args:
            records (list[_GXFRecord]): Buffered records that belong to one streaming chunk.
            gene_records (dict[str, _GXFRecord]): Previously seen gene records keyed by gene id.
            gene_cache (dict[str, Gene]): Previously built Gene objects keyed by gene id.
            candidate_ids (Mapping[str, tuple[int, int]] | None): When given, only transcripts
                whose id is a key of this mapping are built; every other transcript is dropped
                while grouping, before any model is constructed for it.

        Returns:
            list[Transcript]: Transcript objects assembled from the buffered records.
        """

        gene_records = {} if gene_records is None else gene_records
        gene_cache = {} if gene_cache is None else gene_cache
        groups = self._group_records(records, candidate_ids=candidate_ids)
        return self._groups_to_features_with_context(
            groups,
            gene_records=gene_records,
            gene_cache=gene_cache,
        )

    def _groups_to_features_with_context(
        self,
        groups: list[_GXFTranscriptGroup],
        *,
        gene_records: dict[str, _GXFRecord] | None = None,
        gene_cache: dict[str, Gene] | None = None,
    ) -> list[Transcript]:
        """Build already-grouped transcripts while reusing gene and slice models."""

        gene_records = {} if gene_records is None else gene_records
        gene_cache = {} if gene_cache is None else gene_cache
        transcripts = [self._build_transcript_from_record_group(group) for group in groups]

        for group, transcript in zip(groups, transcripts, strict=True):
            gene_record = self._resolve_gene_record(group, gene_records)
            self._resolve_gene(gene_record, group.gene_id, transcript, gene_cache)

        return transcripts

    def _group_records(
        self,
        records: Iterable[_GXFRecord],
        *,
        candidate_ids: Mapping[str, tuple[int, int]] | None = None,
    ) -> list[_GXFTranscriptGroup]:
        """Group parser records into per-transcript blocks.

        ``record_id`` is computed only for gene and transcript rows, and ``parent_id`` only for
        transcript, exon and coding rows -- the ids each routing decision below actually needs. An
        exon or CDS record's own id is still resolved later, per exon, when the transcript itself
        is built.

        Args:
            records (Iterable[_GXFRecord]): The buffered records to group.
            candidate_ids (Mapping[str, tuple[int, int]] | None): When given, only transcripts
                whose id is a key of this mapping are grouped; every other transcript's records are
                still scanned (a record's transcript is not known until the whole pass completes),
                but no group, and therefore no model, is built for it. Relative order among the
                kept transcripts is unchanged.

        Returns:
            list[_GXFTranscriptGroup]: One grouped transcript block per (candidate) transcript
            record, in the order transcripts were first seen.
        """

        gene_map: dict[str, _GXFRecord] = {}
        transcript_map: dict[str, tuple[_GXFRecord, str | None]] = {}
        exon_map: dict[str, list[_GXFRecord]] = defaultdict(list)
        cds_map: dict[str, list[_GXFRecord]] = defaultdict(list)

        for record in records:
            feature_type = record.feature_type.lower()

            if feature_type in self.gene_feature_types:
                record_id = self.get_record_id(record)
                if record_id:
                    gene_map[record_id] = record
            elif feature_type in self.transcript_feature_types:
                record_id = self.get_record_id(record)
                if record_id:
                    previous = transcript_map.get(record_id)
                    if previous is not None and (previous[0].start, previous[0].end) != (
                        record.start,
                        record.end,
                    ):
                        # One row can arrive once per merged span; a different span is a real duplicate.
                        self._disable_pool_for_duplicate_id(record_id)
                    # A later duplicate without a Parent keeps the earlier gene id.
                    parent_ids = self.get_parent_ids(record)
                    parent_id = parent_ids[0] if parent_ids else None
                    previous_gene_id = previous[1] if previous is not None else None
                    transcript_map[record_id] = (record, parent_id or previous_gene_id)
            elif feature_type in self.exon_feature_types:
                parent_ids = self.get_parent_ids(record)
                for parent_id in parent_ids:
                    exon_map[parent_id].append(record)
            elif feature_type in self.coding_feature_types:
                parent_ids = self.get_parent_ids(record)
                for parent_id in parent_ids:
                    cds_map[parent_id].append(record)

        return [
            _GXFTranscriptGroup(
                transcript_id=transcript_id,
                gene_id=gene_id,
                transcript_record=transcript_record,
                gene_record=gene_map.get(gene_id) if gene_id else None,
                exon_records=exon_map.get(transcript_id, []),
                cds_records=cds_map.get(transcript_id, []),
            )
            for transcript_id, (transcript_record, gene_id) in transcript_map.items()
            if candidate_ids is None or transcript_id in candidate_ids
        ]

    def _resolve_gene_record(
        self,
        group: _GXFTranscriptGroup,
        gene_records: dict[str, _GXFRecord] | None = None,
    ) -> _GXFRecord | None:
        """Return the gene record for a grouped transcript, if one exists.

        Args:
            group (_GXFTranscriptGroup): The grouped transcript records.
            gene_records (dict[str, _GXFRecord]): Previously seen gene records keyed by gene id.

        Returns:
            _GXFRecord | None: The gene record for the grouped transcript, or None if not found.
        """

        gene_records = {} if gene_records is None else gene_records
        if group.gene_record is not None:
            return group.gene_record

        if group.gene_id:
            return gene_records.get(group.gene_id)
        return None

    def _resolve_gene(
        self,
        gene_record: _GXFRecord | None,
        gene_id: str | None,
        transcript: Transcript,
        gene_cache: dict[str, Gene] | None = None,
    ) -> None:
        """Attach a Gene to a Transcript, reusing cached Gene objects when possible.

        Args:
            gene_record (_GXFRecord | None): The gene record to convert, if present.
            gene_id (str | None): The gene's own identifier, already resolved by the caller
                (``_GXFTranscriptGroup.gene_id`` for the common in-group case) rather than
                recomputed from ``gene_record``.
            transcript (Transcript): The transcript to attach the gene to.
            gene_cache (dict[str, Gene]): Previously built Gene objects keyed by gene id.
        """

        if gene_record is None:
            return

        gene_cache = {} if gene_cache is None else gene_cache
        gene_key = gene_id or ""
        gene = gene_cache.get(gene_key)
        if gene is None:
            gene = self._build_gene(gene_record)
            gene_cache[gene_key] = gene

        assert gene.transcripts is not None
        gene.transcripts.append(transcript)
        transcript.gene = gene

    def _build_transcript_from_record_group(
        self,
        group: _GXFTranscriptGroup,
    ) -> Transcript:
        """Construct a Transcript from one grouped record block.

        Args:
            group: The grouped records to convert.

        Returns:
            Transcript: The assembled transcript object.
        """

        transcript_record = group.transcript_record
        strand = STRAND_SYMBOL_TO_ENUM.get(transcript_record.strand_symbol, Strand.UNKNOWN)
        ### TO DO: Resolve the seqid to reference seq_region_name to get slice
        genomic_slice = self.slice_cache.get_slice(transcript_record.seqid, self.sequence_source)

        transcript = Transcript(
            stable_id=group.transcript_id,
            start=transcript_record.start,
            end=transcript_record.end,
            strand=strand,
            slice=genomic_slice,
        )
        # Both spellings appear in the wild: GFF3 writes "biotype", GTF "transcript_biotype". It was
        # read for the gene but never for the transcript, so a protein_coding transcript read from a
        # file came back with no biotype at all.
        biotype = self._first_attr(transcript_record, "biotype") or self._first_attr(
            transcript_record, "transcript_biotype"
        )
        exons = self._build_exons(group.exon_records, group.cds_records, strand, genomic_slice)
        # Brand-new transcript: nothing is cached yet, so the on_setattr hook has nothing to
        # invalidate -- skip straight to the fields it would leave set.
        transcript.assign_without_hooks(biotype=biotype, spliced_exons=exons)

        if group.cds_records:
            translation = self._build_translation(group.cds_records, transcript)
            if translation is not None:
                translation.transcript = transcript
                # translation is derived and read-only, so it is seeded rather than assigned.
                # seed_attrs, not seed_and_invalidate_attrs: nothing downstream is resolved yet.
                transcript.seed_attrs(translation=translation)

        return transcript

    def _build_exons(
        self,
        exon_records: list[_GXFRecord],
        cds_records: list[_GXFRecord],
        strand: Strand,
        genomic_slice: Slice,
    ) -> list[Exon]:
        """Construct exon objects from grouped records.

        An exon's phase is read from the first CDS row, in transcript order, that starts inside it,
        as VEP reads it (``BaseGXF::_create_transcript``); an exon no CDS row starts in has an
        unknown phase. A transcript with no CDS rows leaves every phase unset.

        Args:
            exon_records: Exon records attached to one transcript.
            cds_records: Coding records attached to the same transcript.
            strand: The transcript strand.
            genomic_slice: The transcript's own Slice, reused for every exon -- an exon is always
                on the same sequence region as its transcript.

        Returns:
            list[Exon]: The exon objects for the transcript, ordered 5' to 3' along the transcript
            (descending genomic start for a reverse-strand transcript).
        """

        reverse = strand == Strand.REVERSE
        ordered_cds_records = sorted(cds_records, key=lambda record: record.start, reverse=reverse)
        exons: list[Exon] = []
        for exon_record in sorted(exon_records, key=lambda record: record.start, reverse=reverse):
            stable_id = self.get_record_id(exon_record) or ""
            exon = Exon(
                stable_id=stable_id,
                start=exon_record.start,
                end=exon_record.end,
                strand=strand,
                slice=genomic_slice,
            )
            if ordered_cds_records:
                # Brand-new exon: nothing cached yet, same reasoning as spliced_exons above.
                exon.assign_without_hooks(phase=self._exon_phase(exon_record, ordered_cds_records))
            exons.append(exon)
        return exons

    @staticmethod
    def _exon_phase(exon_record: _GXFRecord, ordered_cds_records: list[_GXFRecord]) -> ExonPhase:
        """Return the Ensembl phase of the first CDS row starting inside an exon, else unknown."""
        for cds_record in ordered_cds_records:
            if exon_record.start <= cds_record.start <= exon_record.end:
                if cds_record.phase_value is None:
                    return ExonPhase.UNKNOWN
                return GFF_TO_ENSEMBL_PHASE.get(cds_record.phase_value, ExonPhase.UNKNOWN)
        return ExonPhase.UNKNOWN

    def _build_translation(self, cds_records: list[_GXFRecord], transcript: Transcript) -> Translation | None:
        """Construct a Translation from CDS records.

        ``Translation.start`` and ``end`` are offsets from the start of the transcript, not genomic
        coordinates -- ``Transcript._coding_start`` re-maps them with ``self.start + start - 1``.
        Seeding them with the raw genomic span therefore double-counts the transcript's own start,
        putting the coding boundaries outside the transcript: ``cds`` resolves to None, and every
        coding consequence silently disappears while positional ones survive.

        Args:
            cds_records: Coding records attached to the transcript.
            transcript: The transcript the offsets are relative to, which also gives the strand:
                on the reverse strand the translation starts at the highest genomic coordinate.

        Returns:
            Translation | None: The translation object, or None when no coding
            records are present.
        """

        if not cds_records:
            return None

        cds_start = min(record.start for record in cds_records)
        cds_end = max(record.end for record in cds_records)
        if transcript.strand == Strand.FORWARD:
            start = cds_start - transcript.start + 1
            end = cds_end - transcript.start + 1
        else:
            start = transcript.end - cds_end + 1
            end = transcript.end - cds_start + 1
        stable_id = self._first_attr(cds_records[0], "protein_id") or ""
        return Translation(start=start, end=end, stable_id=stable_id)

    def _build_gene(self, gene_record: _GXFRecord) -> Gene:
        """Construct a Gene object from a grouped gene record.

        Args:
            gene_record: The gene record to convert.

        Returns:
            Gene: The assembled gene object.
        """

        strand = STRAND_SYMBOL_TO_ENUM.get(gene_record.strand_symbol, Strand.UNKNOWN)
        stable_id = self.get_record_id(gene_record) or ""

        gene = Gene(
            stable_id=stable_id,
            start=gene_record.start,
            end=gene_record.end,
            strand=strand,
            slice=self.slice_cache.get_slice(gene_record.seqid, self.sequence_source),
        )
        name = self._first_attr(gene_record, "Name") or self._first_attr(gene_record, "gene_name")
        biotype = self._first_attr(gene_record, "biotype") or self._first_attr(gene_record, "gene_biotype")
        # Brand-new gene: nothing cached yet, same reasoning as the transcript and exon paths above.
        gene.assign_without_hooks(name=name, biotype=biotype, transcripts=[])
        return gene

    @staticmethod
    def _first_attr(record: _GXFRecord, key: str) -> str | None:
        """Return the first attribute value for a key, if present.

        Uses the ``first_attribute`` pymethod -- a linear scan of the stored attribute pairs --
        rather than the ``attributes`` getter, which builds a whole dict per call.

        Args:
            record: The record to inspect.
            key: The attribute key to look up.

        Returns:
            str | None: The first attribute value for the key, if present.
        """

        return record.first_attribute(key)
