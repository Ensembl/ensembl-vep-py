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
* build Transcript, Gene, Exon, Intron, and Translation objects, and
* reuse Slice and Gene instances across records that share the same region.

Only the parser-specific record identity and parent resolution rules are left
to the concrete GFF3 and GTF adapters.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from collections import defaultdict
from typing import ClassVar, Generic, Iterable, Protocol, TypeVar, Iterator

from attrs import define

from ensembl.vep.common import CoordSystemLevel, ExonPhase, PHASE_VALUE_TO_ENUM, STRAND_SYMBOL_TO_ENUM, Strand, format_region_string
from ensembl.vep.models.genic import Exon, Gene, Intron, Transcript, Translation
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.parsers import Parser
from ensembl.vep.sources.file.file_feature_source import FileFeatureSource

logger = logging.getLogger(__name__)

_CHROMOSOME_COORD_SYSTEM = CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME, default=True)
_PLACEHOLDER_SEQ_REGION_END = 2**31 - 1
_ORDERED_FILE_ASSUMPTION_MESSAGE = (
    "GFF3/GTF streaming iteration assumes a sorted file: a gene record precedes its "
    "transcripts, and each transcript's exon/CDS children immediately follow it. "
    "Files violating this ordering may yield incomplete or incorrectly grouped transcripts."
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
    phase_value: str | None
    attributes: dict[str, tuple[str, ...]]


@define(kw_only=True, auto_attribs=True)
class _GXFTranscriptGroup:
    """Grouped parser records describing one transcript.

    Attributes:
        transcript_record (_GXFRecord): Primary transcript record for the group.
        gene_record (_GXFRecord | None): Associated gene record, if one was present in the input.
        exon_records (list[_GXFRecord]): Exon records attached to the transcript.
        cds_records (list[_GXFRecord]): Coding records attached to the transcript.
    """

    transcript_record: _GXFRecord
    gene_record: _GXFRecord | None
    exon_records: list[_GXFRecord]
    cds_records: list[_GXFRecord]


@define(kw_only=True, auto_attribs=True)
class GXFTranscriptSource(FileFeatureSource[P, Transcript], ABC, Generic[P]):
    """Transcript source backed by GFF3 or GTF file.

    Attributes:
        gene_feature_types (frozenset[str]): Lower-cased feature types that identify a gene record.
        transcript_feature_types (frozenset[str]): Lower-cased feature types that identify a transcript record.
        exon_feature_types (frozenset[str]): Lower-cased feature types that identify an exon record.
        coding_feature_types (frozenset[str]): Lower-cased feature types that contribute to the transcript translation.
    """

    gene_feature_types: ClassVar[frozenset[str]] = frozenset({"gene"})
    transcript_feature_types: ClassVar[frozenset[str]] = frozenset({"transcript"})
    exon_feature_types: ClassVar[frozenset[str]] = frozenset({"exon"})
    coding_feature_types: ClassVar[frozenset[str]] = frozenset({"cds"})
    adaptor: P

    def iter_all_features(self) -> Iterator[Transcript]:
        """Yield transcripts by streaming the file sequentially.

        Returns:
            Iterator[Transcript]: Transcripts built from the entire file.
        """

        logger.info(_ORDERED_FILE_ASSUMPTION_MESSAGE)

        with self.adaptor.open(self.file_path, mode="sequential") as adaptor:
            yield from self._iter_features_from_adaptor(adaptor)

    def iter_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> Iterator[Transcript]:
        """Yield transcripts overlapping a genomic region.

        Returns:
            Iterator[Transcript]: Transcripts overlapping the requested region.
        """

        try:
            region = format_region_string(region, seq_region_name=seq_region_name, start=start, end=end)
        except ValueError:
            return

        logger.info(_ORDERED_FILE_ASSUMPTION_MESSAGE)

        with self.adaptor.open(self.file_path, mode="query") as adaptor:
            if not adaptor.seek(region):
                return
            yield from self._iter_features_from_adaptor(adaptor)

    def _iter_features_from_adaptor(self, adaptor: P) -> Iterator[Transcript]:
        """Yield transcripts from an open parser, buffering records into blocks.

        Buffers records from the parser until a new transcript or gene record is seen, then
        converts the buffered records into Transcript objects.

        Args:
            adaptor: An open parser positioned at the records to stream.

        Returns:
            Iterator[Transcript]: Transcripts built from buffered GFF3 records.
        """

        gene_records: dict[str, _GXFRecord] = {}
        gene_cache: dict[str, Gene] = {}
        slice_cache: dict[str, Slice] = {}
        buffer: list[_GXFRecord] = []

        while (record := adaptor.next()) is not None:
            feature_type = record.feature_type.lower()
            if feature_type in self.gene_feature_types:
                raw_id = self.get_record_id(record)
                if raw_id:
                    gene_records[raw_id] = record
                if buffer:
                    yield from self._records_to_features_with_context(
                        buffer,
                        gene_records=gene_records,
                        gene_cache=gene_cache,
                        slice_cache=slice_cache,
                    )
                    buffer = []
                continue

            if feature_type in self.transcript_feature_types and buffer:
                yield from self._records_to_features_with_context(
                    buffer,
                    gene_records=gene_records,
                    gene_cache=gene_cache,
                    slice_cache=slice_cache,
                )
                buffer = []

            buffer.append(record)

        if buffer:
            yield from self._records_to_features_with_context(
                buffer,
                gene_records=gene_records,
                gene_cache=gene_cache,
                slice_cache=slice_cache,
            )

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
    def get_parent_id(self, record: _GXFRecord) -> str | None:
        """Return the parent identifier for a GFF3/GTF record, if one exists.

        Args:
            record: The record to inspect.

        Returns:
            str | None: The parent identifier, or None when the record has no
            parent.
        """

        raise NotImplementedError

    def _records_to_features_with_context(
        self,
        records: list[_GXFRecord],
        *,
        gene_records: dict[str, _GXFRecord] = {},
        gene_cache: dict[str, Gene] = {},
        slice_cache: dict[str, Slice] = {},
    ) -> list[Transcript]:
        """Convert buffered records into Transcript objects while reusing caches.

        Args:
            records (list[_GXFRecord]): Buffered records that belong to one streaming chunk.
            gene_records (dict[str, _GXFRecord]): Previously seen gene records keyed by gene id.
            gene_cache (dict[str, Gene]): Previously built Gene objects keyed by gene id.
            slice_cache (dict[str, Slice]): Previously built Slice objects keyed by sequence-region name.

        Returns:
            list[Transcript]: Transcript objects assembled from the buffered records.
        """

        groups = self._group_records(records)
        transcripts = [self._build_transcript_from_record_group(group, slice_cache) for group in groups]

        for group, transcript in zip(groups, transcripts):
            gene_record = self._resolve_gene_record(group, gene_records)
            self._resolve_gene(gene_record, transcript, gene_cache, slice_cache)

        return transcripts

    def _group_records(self, records: Iterable[_GXFRecord]) -> list[_GXFTranscriptGroup]:
        """Group parser records into per-transcript blocks.

        Args:
            records (Iterable[_GXFRecord]): The buffered records to group.

        Returns:
            list[_GXFTranscriptGroup]: One grouped transcript block per transcript record.
        """

        gene_map: dict[str, _GXFRecord] = {}
        transcript_map: dict[str, _GXFRecord] = {}
        transcript_parent: dict[str, str] = {}
        exon_map: dict[str, list[_GXFRecord]] = defaultdict(list)
        cds_map: dict[str, list[_GXFRecord]] = defaultdict(list)

        for record in records:
            feature_type = record.feature_type.lower()
            record_id = self.get_record_id(record)
            parent_id = self.get_parent_id(record)

            if feature_type in self.gene_feature_types:
                if record_id:
                    gene_map[record_id] = record
            elif feature_type in self.transcript_feature_types:
                if record_id:
                    transcript_map[record_id] = record
                    if parent_id:
                        transcript_parent[record_id] = parent_id
            elif feature_type in self.exon_feature_types:
                if parent_id:
                    exon_map[parent_id].append(record)
            elif feature_type in self.coding_feature_types:
                if parent_id:
                    cds_map[parent_id].append(record)

        return [
            _GXFTranscriptGroup(
                transcript_record=transcript_record,
                gene_record=gene_map.get(transcript_parent.get(transcript_id, ""), None),
                exon_records=exon_map.get(transcript_id, []),
                cds_records=cds_map.get(transcript_id, []),
            )
            for transcript_id, transcript_record in transcript_map.items()
        ]

    def _resolve_gene_record(
        self,
        group: _GXFTranscriptGroup,
        gene_records: dict[str, _GXFRecord]  = {},
    ) -> _GXFRecord | None:
        """Return the gene record for a grouped transcript, if one exists.
        
        Args:
            group (_GXFTranscriptGroup): The grouped transcript records.
            gene_records (dict[str, _GXFRecord]): Previously seen gene records keyed by gene id.

        Returns:
            _GXFRecord | None: The gene record for the grouped transcript, or None if not found.
        """

        gene_record = group.gene_record
        if gene_record is not None:
            return gene_record

        parent = self.get_parent_id(group.transcript_record)
        if parent:
            return gene_records.get(parent)
        return None

    def _resolve_gene(
        self,
        gene_record: _GXFRecord | None,
        transcript: Transcript,
        gene_cache: dict[str, Gene] = {},
        slice_cache: dict[str, Slice] = {},
    ) -> None:
        """Attach a Gene to a Transcript, reusing cached Gene objects when possible.
        
        Args:
            gene_record (_GXFRecord | None): The gene record to convert, if present.
            transcript (Transcript): The transcript to attach the gene to.
            gene_cache (dict[str, Gene]): Previously built Gene objects keyed by gene id.
            slice_cache (dict[str, Slice]): Previously built Slice objects keyed by sequence-region name.
        """

        if gene_record is None:
            return

        gene_id = self.get_record_id(gene_record) or ""
        gene = gene_cache.get(gene_id)
        if gene is None:
            gene = self._build_gene(gene_record, slice_cache)
            gene_cache[gene_id] = gene

        gene.transcripts.append(transcript)
        transcript.gene = gene

    def _build_transcript_from_record_group(
        self,
        group: _GXFTranscriptGroup,
        slice_cache: dict[str, Slice],
    ) -> Transcript:
        """Construct a Transcript from one grouped record block.

        Args:
            group: The grouped records to convert.
            slice_cache: Previously built Slice cache keyed by sequence-region name.

        Returns:
            Transcript: The assembled transcript object.
        """

        transcript_record = group.transcript_record
        strand = STRAND_SYMBOL_TO_ENUM.get(transcript_record.strand_symbol, Strand.UNKNOWN)
        genomic_slice = self._get_or_build_slice(transcript_record.seqid, slice_cache)
        stable_id = self.get_record_id(transcript_record)

        transcript = Transcript(
            stable_id=stable_id,
            start=transcript_record.start,
            end=transcript_record.end,
            strand=strand,
            slice=genomic_slice,
        )
        transcript.spliced_exons = self._build_exons(group.exon_records, strand, slice_cache)
        transcript.introns = self._build_introns(transcript.spliced_exons, strand, genomic_slice)

        if group.cds_records:
            translation = self._build_translation(group.cds_records)
            if translation is not None:
                translation.transcript = transcript
                transcript.translation = translation

        return transcript

    def _build_exons(
        self,
        exon_records: list[_GXFRecord],
        strand: Strand,
        slice_cache: dict[str, Slice],
    ) -> list[Exon]:
        """Construct exon objects from grouped records.

        Args:
            exon_records: Exon records attached to one transcript.
            strand: The transcript strand.
            slice_cache: Previously built Slice cache keyed by sequence-region name.

        Returns:
            list[Exon]: The exon objects for the transcript, ordered 5' to 3' along the transcript
            (descending genomic start for a reverse-strand transcript).
        """

        exons: list[Exon] = []
        for exon_record in sorted(
            exon_records, key=lambda record: record.start, reverse=(strand == Strand.REVERSE)
        ):
            stable_id = self.get_record_id(exon_record)
            genomic_slice = self._get_or_build_slice(exon_record.seqid, slice_cache)
            exon = Exon(
                stable_id=stable_id,
                start=exon_record.start,
                end=exon_record.end,
                strand=strand,
                slice=genomic_slice,
            )
            if exon_record.phase_value is not None:
                exon.phase = PHASE_VALUE_TO_ENUM.get(exon_record.phase_value, ExonPhase.UNKNOWN)
            exons.append(exon)
        return exons

    def _build_introns(self, spliced_exons: list[Exon], strand: Strand, genomic_slice: Slice) -> list[Intron]:
        """Construct introns between consecutive exons.

        Args:
            spliced_exons: Exons ordered 5' to 3' along the transcript.
            strand: The transcript strand.
            genomic_slice: The shared chromosome-level slice.

        Returns:
            list[Intron]: The introns between adjacent exons.
        """

        introns: list[Intron] = []
        for prev_exon, next_exon in zip(spliced_exons, spliced_exons[1:]):
            first, second = (next_exon, prev_exon) if strand == Strand.REVERSE else (prev_exon, next_exon)
            intron = Intron(
                start=first.end + 1,
                end=second.start - 1,
                strand=strand,
                slice=genomic_slice,
            )
            intron.prev_exon = prev_exon
            intron.next_exon = next_exon
            introns.append(intron)
        return introns

    def _build_translation(self, cds_records: list[_GXFRecord]) -> Translation | None:
        """Construct a Translation from CDS records.

        Args:
            cds_records: Coding records attached to the transcript.

        Returns:
            Translation | None: The translation object, or None when no coding
            records are present.
        """

        if not cds_records:
            return None

        cds_start = min(record.start for record in cds_records)
        cds_end = max(record.end for record in cds_records)
        stable_id = self._first_attr(cds_records[0], "protein_id") or ""
        return Translation(start=cds_start, end=cds_end, stable_id=stable_id)

    def _build_gene(self, gene_record: _GXFRecord, slice_cache: dict[str, Slice]) -> Gene:
        """Construct a Gene object from a grouped gene record.

        Args:
            gene_record: The gene record to convert.
            slice_cache: Previously built Slice cache keyed by sequence-region name.

        Returns:
            Gene: The assembled gene object.
        """

        strand = STRAND_SYMBOL_TO_ENUM.get(gene_record.strand_symbol, Strand.UNKNOWN)
        stable_id = self.get_record_id(gene_record)

        gene = Gene(
            stable_id=stable_id,
            start=gene_record.start,
            end=gene_record.end,
            strand=strand,
            slice=self._get_or_build_slice(gene_record.seqid, slice_cache),
        )
        gene.name = self._first_attr(gene_record, "Name") or self._first_attr(gene_record, "gene_name")
        gene.biotype = self._first_attr(gene_record, "biotype") or self._first_attr(gene_record, "gene_biotype")
        gene.transcripts = []
        return gene

    def _get_or_build_slice(self, seq_region_name: str, slice_cache: dict[str, Slice]) -> Slice:
        """Return the shared Slice for a sequence region, creating it when needed.

        Args:
            seq_region_name: The sequence-region name to look up.
            slice_cache: Previously built Slice cache keyed by sequence-region name.

        Returns:
            Slice: The shared Slice object for the sequence region.
        """

        if seq_region_name not in slice_cache:
            seq_region_length = None
            if self.sequence_source is not None:
                seq_region_length = self.sequence_source.get_seq_region_length(seq_region_name)

            slice_cache[seq_region_name] = Slice(
                coord_system=_CHROMOSOME_COORD_SYSTEM,
                seq_region_name=seq_region_name,
                start=1,
                end=seq_region_length if seq_region_length is not None else _PLACEHOLDER_SEQ_REGION_END,
                strand=Strand.FORWARD,
                seq_region_length=seq_region_length,
                sequence_source=self.sequence_source,
            )
        return slice_cache[seq_region_name]

    @staticmethod
    def _first_attr(record: _GXFRecord, key: str) -> str | None:
        """Return the first attribute value for a key, if present.
        
        Args:
            record: The record to inspect.
            key: The attribute key to look up.

        Returns:
            str | None: The first attribute value for the key, if present.
        """

        values = record.attributes.get(key, [None])
        return values[0]
