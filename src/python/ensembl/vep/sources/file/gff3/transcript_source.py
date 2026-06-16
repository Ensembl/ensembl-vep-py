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
"""GFF3-backed transcript feature source."""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Iterable, Iterator

from attrs import define, field

from ensembl.vep.common import ExonPhase, Strand, CoordSystemLevel, STRAND_SYMBOL_TO_ENUM, PHASE_VALUE_TO_ENUM
from ensembl.vep.common.methods import format_region_string
from ensembl.vep.parsers import GFF3Parser, GFF3Record
from ensembl.vep.models.genic import Transcript, Exon, Gene, Translation, Intron
from ensembl.vep.models.sequence import Slice, CoordSystem
from ensembl.vep.sources.file.file_source import FileSource

logger = logging.getLogger(__name__)

_CHROMOSOME_COORD_SYSTEM = CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME, default=True)
_PLACEHOLDER_SEQ_REGION_END = 2**31 - 1  # placeholder until real sequence-region lengths are available

_TRANSCRIPT_FEATURE_TYPES = frozenset(
    {"mRNA", "transcript", "lnc_RNA", "pseudogenic_transcript", "processed_transcript"}
)

_ORDERED_FILE_ASSUMPTION_MESSAGE = (
    "GFF3 streaming iteration assumes a sorted file: a gene record precedes its "
    "transcripts, and each transcript's exon/CDS children immediately follow it. "
    "Files violating this ordering may yield incomplete or incorrectly grouped transcripts."
)


@define(kw_only=True, auto_attribs=True)
class _GFF3TranscriptGroup:
    """Groups raw GFF3 records that together describe one transcript.

    Attributes:
        transcript_record (GFF3Record): The main transcript record.
        gene_record (GFF3Record | None): The associated gene record, if any.
        exon_records (list[GFF3Record]): The exon records for this transcript.
        cds_records (list[GFF3Record]): The CDS records for this transcript.
    """

    transcript_record: GFF3Record
    gene_record: GFF3Record | None
    exon_records: list[GFF3Record]
    cds_records: list[GFF3Record]


@define(kw_only=True, auto_attribs=True)
class GFF3TranscriptSource(FileSource[GFF3Parser, Transcript]):
    """Reads a GFF3 file and produces Transcript domain objects.

    Handles GFF3 parent-child grouping: one mRNA record plus its exon and CDS
    children (linked via ID/Parent attributes) are assembled into one Transcript.
    CDS records are also used to build a Translation, exons are used to build
    Introns, and gene records are used to build a Gene linked back to each of
    its transcripts.

    Attributes:
        adaptor (GFF3Parser): The GFF3 parser adaptor.

    Usage::

        with GFF3TranscriptSource(file_path="/path/to/annotation.gff3.gz") as source:
            transcripts = source.get_features_by_region("1", start=1_000_000, end=2_000_000)
    """

    adaptor: GFF3Parser = field(init=False, factory=GFF3Parser)

    def _records_to_features(self, records: list[GFF3Record]) -> list[Transcript]:
        """Convert a list of GFF3 records into Transcript domain objects.

        Args:
            records (list[GFF3Record]): The GFF3 records to convert.

        Returns:
            list[Transcript]: A list of Transcript domain objects.
        """

        groups = _group_gff3_records(records)
        gene_objects: dict[str, Gene] = {}
        slice_cache: dict[str, Slice] = {}
        transcripts = [_build_transcript_from_record_group(group, slice_cache) for group in groups]

        for group, transcript in zip(groups, transcripts):
            _resolve_gene(group.gene_record, transcript, gene_objects, slice_cache)

        return transcripts

    def iter_all_features(self) -> Iterator[Transcript]:
        """Lazily yield Transcripts, assuming the GFF3 file is sorted.

        See module-level _ORDERED_FILE_ASSUMPTION_MESSAGE for the ordering this relies on.
        """

        logger.info(_ORDERED_FILE_ASSUMPTION_MESSAGE)
        gene_records: dict[str, GFF3Record] = {}
        gene_objects: dict[str, Gene] = {}
        slice_cache: dict[str, Slice] = {}
        buffer: list[GFF3Record] = []

        with self.adaptor.open(self.file_path, mode="sequential") as adaptor:
            while (record := adaptor.next()) is not None:
                if record.feature_type == "gene":
                    raw_id = _first_attr(record, "ID")
                    if raw_id:
                        gene_records[raw_id] = record
                    if buffer:
                        yield from _flush_buffer(buffer, gene_records, gene_objects, slice_cache)
                        buffer = []
                    continue
                if record.feature_type in _TRANSCRIPT_FEATURE_TYPES and buffer:
                    yield from _flush_buffer(buffer, gene_records, gene_objects, slice_cache)
                    buffer = []
                buffer.append(record)
            if buffer:
                yield from _flush_buffer(buffer, gene_records, gene_objects, slice_cache)

    def iter_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> Iterator[Transcript]:
        """Lazily yield Transcripts overlapping a region, assuming the GFF3 file is sorted.

        See module-level _ORDERED_FILE_ASSUMPTION_MESSAGE for the ordering this relies on.
        """

        try:
            region = format_region_string(region, seq_region_name=seq_region_name, start=start, end=end)
        except ValueError:
            return

        logger.info(_ORDERED_FILE_ASSUMPTION_MESSAGE)
        gene_records: dict[str, GFF3Record] = {}
        gene_objects: dict[str, Gene] = {}
        slice_cache: dict[str, Slice] = {}
        buffer: list[GFF3Record] = []

        with self.adaptor.open(self.file_path) as adaptor:
            if not adaptor.seek(region):
                return
            while (record := adaptor.next()) is not None:
                if record.feature_type == "gene":
                    raw_id = _first_attr(record, "ID")
                    if raw_id:
                        gene_records[raw_id] = record
                    if buffer:
                        yield from _flush_buffer(buffer, gene_records, gene_objects, slice_cache)
                        buffer = []
                    continue
                if record.feature_type in _TRANSCRIPT_FEATURE_TYPES and buffer:
                    yield from _flush_buffer(buffer, gene_records, gene_objects, slice_cache)
                    buffer = []
                buffer.append(record)
            if buffer:
                yield from _flush_buffer(buffer, gene_records, gene_objects, slice_cache)


def _flush_buffer(
    buffer: list[GFF3Record],
    gene_records: dict[str, GFF3Record],
    gene_objects: dict[str, Gene],
    slice_cache: dict[str, Slice],
) -> list[Transcript]:
    """Build Transcripts from one buffered block and resolve their genes via persistent caches."""

    groups = _group_gff3_records(buffer)
    transcripts = [_build_transcript_from_record_group(group, slice_cache) for group in groups]

    for group, transcript in zip(groups, transcripts):
        gene_record = group.gene_record
        if gene_record is None and group.transcript_record is not None:
            parent = _first_attr(group.transcript_record, "Parent")
            if parent:
                gene_record = gene_records.get(parent)
        _resolve_gene(gene_record, transcript, gene_objects, slice_cache)

    return transcripts


def _resolve_gene(
    gene_record: GFF3Record | None,
    transcript: Transcript,
    gene_objects: dict[str, Gene],
    slice_cache: dict[str, Slice],
) -> None:
    """Link transcript to its Gene, building or reusing the Gene from the given caches."""

    if gene_record is None:
        return
    gene_id = _first_attr(gene_record, "ID") or ""
    gene = gene_objects.get(gene_id)
    if gene is None:
        gene = _build_gene(gene_record, slice_cache)
        gene_objects[gene_id] = gene
    gene.transcripts.append(transcript)
    transcript.gene = gene


def _group_gff3_records(records: Iterable[GFF3Record]) -> list[_GFF3TranscriptGroup]:
    """Group GFF3 records by parent-child relationship into per-transcript groups.

    Args:
        records (Iterable[GFF3Record]): The GFF3 records to group.

    Returns:
        list[_GFF3TranscriptGroup]: A list of grouped transcript records.
    """

    gene_map: dict[str, GFF3Record] = {}
    transcript_map: dict[str, GFF3Record] = {}
    transcript_parent: dict[str, str] = {}
    exon_map: dict[str, list[GFF3Record]] = defaultdict(list)
    cds_map: dict[str, list[GFF3Record]] = defaultdict(list)

    for record in records:
        feature_type = record.feature_type
        if feature_type == "gene":
            raw_id = _first_attr(record, "ID")
            if raw_id:
                gene_map[raw_id] = record
        elif feature_type in _TRANSCRIPT_FEATURE_TYPES:
            raw_id = _first_attr(record, "ID")
            if raw_id:
                transcript_map[raw_id] = record
                parent = _first_attr(record, "Parent")
                if parent:
                    transcript_parent[raw_id] = parent
        elif feature_type == "exon":
            parent = _first_attr(record, "Parent")
            if parent:
                exon_map[parent].append(record)
        elif feature_type == "CDS":
            parent = _first_attr(record, "Parent")
            if parent:
                cds_map[parent].append(record)

    return [
        _GFF3TranscriptGroup(
            transcript_record=transcript_record,
            gene_record=gene_map.get(transcript_parent.get(transcript_id, ""), None),
            exon_records=exon_map.get(transcript_id, []),
            cds_records=cds_map.get(transcript_id, []),
        )
        for transcript_id, transcript_record in transcript_map.items()
    ]


def _get_or_build_slice(seq_region_name: str, slice_cache: dict[str, Slice]) -> Slice:
    """Return the shared Slice for a seq_region_name, building it on first use.

    The Slice represents the whole sequence region (e.g. chromosome), not any one
    feature's own span, so it is shared across every feature built on that region
    within one call. start=1 and end=_PLACEHOLDER_SEQ_REGION_END are placeholders
    until real sequence-region lengths are available.
    """

    if seq_region_name not in slice_cache:
        slice_cache[seq_region_name] = Slice(
            coord_system=_CHROMOSOME_COORD_SYSTEM,
            seq_region_name=seq_region_name,
            start=1,
            end=_PLACEHOLDER_SEQ_REGION_END,
            strand=Strand.FORWARD,
        )
    return slice_cache[seq_region_name]


def _build_transcript_from_record_group(group: _GFF3TranscriptGroup, slice_cache: dict[str, Slice]) -> Transcript:
    """Construct a Transcript domain object from a grouped set of GFF3 records.

    Args:
        group (_GFF3TranscriptGroup): The grouped GFF3 records for a transcript.
        slice_cache (dict[str, Slice]): Shared Slice cache for this call, keyed by seq_region_name.

    Returns:
        Transcript: The constructed Transcript domain object.
    """

    transcript_record = group.transcript_record
    strand = STRAND_SYMBOL_TO_ENUM.get(transcript_record.strand_symbol, Strand.UNKNOWN)
    genomic_slice = _get_or_build_slice(transcript_record.seqid, slice_cache)
    raw_id = _first_attr(transcript_record, "ID") or ""
    stable_id = raw_id.split(":")[-1] if ":" in raw_id else raw_id

    transcript = Transcript(
        stable_id=stable_id,
        start=transcript_record.start,
        end=transcript_record.end,
        strand=strand,
        slice=genomic_slice,
    )
    transcript.spliced_exons = _build_exons(group.exon_records, strand, slice_cache)
    transcript.introns = _build_introns(transcript.spliced_exons, strand, genomic_slice)

    if group.cds_records:
        translation = _build_translation(group.cds_records)
        if translation:
            translation.transcript = transcript
            transcript.translation = translation

    return transcript


def _build_exons(exon_records: list[GFF3Record], strand: Strand, slice_cache: dict[str, Slice]) -> list[Exon]:
    """Construct Exon domain objects from a transcript's exon records, sorted by genomic start.

    Records are sorted before building Exon objects rather than after, since sorting
    plain GFF3Record.start (int) is cheaper than sorting by an Exon's .start attribute.

    Args:
        exon_records (list[GFF3Record]): The GFF3 exon records for a transcript.
        strand (Strand): The strand of the transcript.
        slice_cache (dict[str, Slice]): Shared Slice cache for this call, keyed by seq_region_name.

    Returns:
        list[Exon]: The constructed Exon domain objects, sorted by genomic start.
    """

    exons = []
    for exon_record in sorted(exon_records, key=lambda record: record.start):
        raw_id = _first_attr(exon_record, "ID") or ""
        stable_id = raw_id.split(":")[-1] if ":" in raw_id else (raw_id or f"exon_{exon_record.start}_{exon_record.end}")
        genomic_slice = _get_or_build_slice(exon_record.seqid, slice_cache)
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


def _build_introns(spliced_exons: list[Exon], strand: Strand, genomic_slice: Slice) -> list[Intron]:
    """Construct Introns between consecutive exons of a transcript.

    spliced_exons is stored genomic-ascending regardless of strand; for reverse-strand
    transcripts the exons are walked in reverse so prev_exon/next_exon reflect the
    transcript's biological reading direction, while each Intron's own start/end stay
    genomic (first <= second).

    Args:
        spliced_exons (list[Exon]): The transcript's exons, sorted genomic-ascending.
        strand (Strand): The strand of the transcript.
        genomic_slice (Slice): The shared Slice for the transcript's seq_region.

    Returns:
        list[Intron]: The introns between consecutive exons, in transcript order.
    """

    ordered = spliced_exons if strand != Strand.REVERSE else list(reversed(spliced_exons))
    introns: list[Intron] = []
    for prev_exon, next_exon in zip(ordered, ordered[1:]):
        first, second = (next_exon, prev_exon) if strand == Strand.REVERSE else (prev_exon, next_exon)
        intron = Intron(
            start=first.end + 1,
            end=second.start - 1,
            strand=strand,
            slice=genomic_slice
        )
        intron.prev_exon = prev_exon
        intron.next_exon = next_exon
        introns.append(intron)
    return introns


def _build_translation(cds_records: list[GFF3Record]) -> Translation | None:
    """Construct a Translation domain object from a transcript's CDS records.

    start/end are genomic coordinates (an approximation of the cDNA-relative
    coordinates Translation conceptually expects); proper coordinate conversion
    is a future improvement.

    Args:
        cds_records (list[GFF3Record]): The CDS records for the transcript.

    Returns:
        Translation | None: The constructed Translation domain object, or None if no CDS records are provided.
    """

    if not cds_records:
        return None
    cds_start = min(record.start for record in cds_records)
    cds_end = max(record.end for record in cds_records)
    stable_id = _first_attr(cds_records[0], "protein_id") or ""
    return Translation(
        start=cds_start,
        end=cds_end,
        stable_id=stable_id
    )


def _build_gene(gene_record: GFF3Record, slice_cache: dict[str, Slice]) -> Gene:
    """Construct a Gene domain object from a gene record.

    Args:
        gene_record (GFF3Record): The GFF3 gene record.
        slice_cache (dict[str, Slice]): Shared Slice cache for this call, keyed by seq_region_name.

    Returns:
        Gene: The constructed Gene domain object, with an empty transcripts list.
    """

    strand = STRAND_SYMBOL_TO_ENUM.get(gene_record.strand_symbol, Strand.UNKNOWN)
    raw_id = _first_attr(gene_record, "ID") or ""
    stable_id = raw_id.split(":")[-1] if ":" in raw_id else raw_id

    gene = Gene(
        stable_id=stable_id,
        start=gene_record.start,
        end=gene_record.end,
        strand=strand,
        slice=_get_or_build_slice(gene_record.seqid, slice_cache),
    )
    gene.name = _first_attr(gene_record, "Name")
    gene.biotype = _first_attr(gene_record, "biotype")
    gene.transcripts = []

    return gene


def _first_attr(record: GFF3Record, key: str) -> str | None:
    """Return the first value of a GFF3 attribute, or None if absent.

    Args:
        record (GFF3Record): The GFF3 record.
        key (str): The attribute key to retrieve.

    Returns:
        str | None: The first value of the attribute, or None if absent.
    """

    values = record.attributes.get(key)
    if values:
        return values[0]
    return None
