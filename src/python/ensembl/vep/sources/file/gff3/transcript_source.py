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

from collections import defaultdict
from typing import Iterable
from attrs import define, field

from ensembl.vep.common import ExonPhase, Strand, STRAND_SYMBOL_TO_ENUM, PHASE_VALUE_TO_ENUM
from ensembl.vep.parsers import GFF3Parser, GFF3Record
from ensembl.vep.models.genic import Transcript, Exon, Gene, Translation
from ensembl.vep.models.sequence.slice import Slice
from ensembl.vep.sources.file.file_source import FileSource

_COORD_SYSTEM = "chromosome"

_TRANSCRIPT_FEATURE_TYPES = frozenset(
    {"mRNA", "transcript", "lnc_RNA", "pseudogenic_transcript", "processed_transcript"}
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
    CDS records are also used to build a Translation, and gene records are used
    to build a Gene linked back to each of its transcripts.

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
        transcripts = [
            _build_transcript_from_record_group(group) for group in groups
        ]

        gene_buckets: dict[str, tuple[GFF3Record, list[Transcript]]] = {}
        for group, transcript in zip(groups, transcripts):
            if group.gene_record:
                gene_id = _first_attr(group.gene_record, "ID") or ""
                gene_buckets.setdefault(gene_id, (group.gene_record, []))[1].append(transcript)

        for gene_rec, gene_transcripts in gene_buckets.values():
            _build_gene(gene_rec, gene_transcripts)

        return transcripts


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


def _build_transcript_from_record_group(group: _GFF3TranscriptGroup) -> Transcript:
    """Construct a Transcript domain object from a grouped set of GFF3 records.

    Args:
        group (_GFF3TranscriptGroup): The grouped GFF3 records for a transcript.

    Returns:
        Transcript: The constructed Transcript domain object.
    """

    transcript_record = group.transcript_record
    strand = STRAND_SYMBOL_TO_ENUM.get(transcript_record.strand_symbol, Strand.UNKNOWN)
    genomic_slice = Slice(
        coord_system=_COORD_SYSTEM,
        seq_region_name=transcript_record.seqid,
        start=transcript_record.start,
        end=transcript_record.end,
        strand=strand,
    )
    raw_id = _first_attr(transcript_record, "ID") or ""
    stable_id = raw_id.split(":")[-1] if ":" in raw_id else raw_id

    transcript = Transcript(
        stable_id=stable_id,
        start=transcript_record.start,
        end=transcript_record.end,
        strand=strand,
        slice=genomic_slice,
    )
    transcript.spliced_exons = sorted(
        [_build_exon(exon_record, strand) for exon_record in group.exon_records],
        key=lambda exon: exon.start,
    )

    if group.cds_records:
        translation = _build_translation(group.cds_records)
        if translation:
            translation.transcript = transcript
            transcript.translation = translation

    return transcript


def _build_exon(exon_record: GFF3Record, strand: Strand) -> Exon:
    """Construct an Exon domain object from a single GFF3 exon record.

    Args:
        exon_record (GFF3Record): The GFF3 exon record.
        strand (Strand): The strand of the exon.

    Returns:
        Exon: The constructed Exon domain object.
    """

    raw_id = _first_attr(exon_record, "ID") or ""
    stable_id = raw_id.split(":")[-1] if ":" in raw_id else (raw_id or f"exon_{exon_record.start}_{exon_record.end}")
    genomic_slice = Slice(
        coord_system=_COORD_SYSTEM,
        seq_region_name=exon_record.seqid,
        start=exon_record.start,
        end=exon_record.end,
        strand=strand,
    )
    exon = Exon(
        stable_id=stable_id,
        start=exon_record.start,
        end=exon_record.end,
        strand=strand,
        slice=genomic_slice,
    )
    if exon_record.phase_value is not None:
        exon.phase = PHASE_VALUE_TO_ENUM.get(exon_record.phase_value, ExonPhase.UNKNOWN)
    return exon


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
    return Translation(start=cds_start, end=cds_end, stable_id=stable_id)


def _build_gene(gene_record: GFF3Record, transcripts: list[Transcript]) -> Gene:
    """Construct a Gene domain object from a gene record and link it to its transcripts.

    Args:
        gene_record (GFF3Record): The GFF3 gene record.
        transcripts (list[Transcript]): The transcripts associated with the gene.

    Returns:
        Gene: The constructed Gene domain object.
    """

    strand = STRAND_SYMBOL_TO_ENUM.get(gene_record.strand_symbol, Strand.UNKNOWN)
    raw_id = _first_attr(gene_record, "ID") or ""
    stable_id = raw_id.split(":")[-1] if ":" in raw_id else raw_id

    gene = Gene(
        stable_id=stable_id,
        start=gene_record.start,
        end=gene_record.end,
        strand=strand,
        slice=Slice(
            coord_system=_COORD_SYSTEM,
            seq_region_name=gene_record.seqid,
            start=gene_record.start,
            end=gene_record.end,
            strand=strand,
        ),
    )
    gene.name = _first_attr(gene_record, "Name")
    gene.biotype = _first_attr(gene_record, "biotype")
    gene.transcripts = transcripts

    for t in transcripts:
        t.gene = gene

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
