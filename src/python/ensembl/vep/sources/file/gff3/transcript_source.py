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
from ensembl.vep.ensembl_vep_core import GFF3Parser, GFF3Record
from ensembl.vep.models.genic import Transcript, Exon, Gene, Translation
from ensembl.vep.models.sequence.slice import Slice
from ensembl.vep.sources.file.file_source import FileSource

_COORD_SYSTEM = "chromosome"

_TRANSCRIPT_FEATURE_TYPES = frozenset(
    {"mRNA", "transcript", "lnc_RNA", "pseudogenic_transcript", "processed_transcript"}
)


@define(kw_only=True, auto_attribs=True)
class _GFF3TranscriptGroup:
    """Groups raw GFF3 records that together describe one transcript."""

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

    Usage::

        with GFF3TranscriptSource(file_path="/path/to/annotation.gff3.gz") as source:
            transcripts = source.get_features_by_region("1", start=1_000_000, end=2_000_000)
    """

    adaptor: GFF3Parser = field(init=False, factory=GFF3Parser)

    def _records_to_features(self, records: list[GFF3Record]) -> list[Transcript]:
        groups = _group_gff3_records(records)
        transcripts = [_build_transcript(g) for g in groups]

        gene_buckets: dict[str, tuple[GFF3Record, list[Transcript]]] = {}
        for group, transcript in zip(groups, transcripts):
            if group.gene_record:
                gene_id = _first_attr(group.gene_record, "ID") or ""
                gene_buckets.setdefault(gene_id, (group.gene_record, []))[1].append(transcript)

        for gene_rec, gene_transcripts in gene_buckets.values():
            _build_gene(gene_rec, gene_transcripts)

        return transcripts


def _group_gff3_records(records: Iterable[GFF3Record]) -> list[_GFF3TranscriptGroup]:
    """Group GFF3 records by parent-child relationship into per-transcript groups."""
    gene_map: dict[str, GFF3Record] = {}
    transcript_map: dict[str, GFF3Record] = {}
    transcript_parent: dict[str, str] = {}
    exon_map: dict[str, list[GFF3Record]] = defaultdict(list)
    cds_map: dict[str, list[GFF3Record]] = defaultdict(list)

    for rec in records:
        ft = rec.feature_type
        if ft == "gene":
            raw_id = _first_attr(rec, "ID")
            if raw_id:
                gene_map[raw_id] = rec
        elif ft in _TRANSCRIPT_FEATURE_TYPES:
            raw_id = _first_attr(rec, "ID")
            if raw_id:
                transcript_map[raw_id] = rec
                parent = _first_attr(rec, "Parent")
                if parent:
                    transcript_parent[raw_id] = parent
        elif ft == "exon":
            parent = _first_attr(rec, "Parent")
            if parent:
                exon_map[parent].append(rec)
        elif ft == "CDS":
            parent = _first_attr(rec, "Parent")
            if parent:
                cds_map[parent].append(rec)

    return [
        _GFF3TranscriptGroup(
            transcript_record=t_rec,
            gene_record=gene_map.get(transcript_parent.get(t_id, ""), None),
            exon_records=exon_map.get(t_id, []),
            cds_records=cds_map.get(t_id, []),
        )
        for t_id, t_rec in transcript_map.items()
    ]


def _build_transcript(group: _GFF3TranscriptGroup) -> Transcript:
    """Construct a Transcript domain object from a grouped set of GFF3 records."""
    t = group.transcript_record
    strand = STRAND_SYMBOL_TO_ENUM.get(t.strand_symbol, Strand.UNKNOWN)
    genomic_slice = Slice(
        coord_system=_COORD_SYSTEM,
        seq_region_name=t.seqid,
        start=t.start,
        end=t.end,
        strand=strand,
    )
    raw_id = _first_attr(t, "ID") or ""
    stable_id = raw_id.split(":")[-1] if ":" in raw_id else raw_id

    transcript = Transcript(
        stable_id=stable_id,
        start=t.start,
        end=t.end,
        strand=strand,
        slice=genomic_slice,
    )
    transcript.spliced_exons = sorted(
        [_build_exon(e, strand) for e in group.exon_records],
        key=lambda e: e.start,
    )

    if group.cds_records:
        translation = _build_translation(group.cds_records)
        if translation:
            translation.transcript = transcript
            transcript.translation = translation

    return transcript


def _build_exon(rec: GFF3Record, strand: Strand) -> Exon:
    """Construct an Exon domain object from a single GFF3 exon record."""
    raw_id = _first_attr(rec, "ID") or ""
    stable_id = raw_id.split(":")[-1] if ":" in raw_id else (raw_id or f"exon_{rec.start}_{rec.end}")
    genomic_slice = Slice(
        coord_system=_COORD_SYSTEM,
        seq_region_name=rec.seqid,
        start=rec.start,
        end=rec.end,
        strand=strand,
    )
    exon = Exon(
        stable_id=stable_id,
        start=rec.start,
        end=rec.end,
        strand=strand,
        slice=genomic_slice,
    )
    if rec.phase_value is not None:
        exon.phase = PHASE_VALUE_TO_ENUM.get(rec.phase_value, ExonPhase.UNKNOWN)
    return exon


def _build_translation(cds_records: list[GFF3Record]) -> Translation | None:
    """Construct a Translation domain object from a transcript's CDS records.

    start/end are genomic coordinates (an approximation of the cDNA-relative
    coordinates Translation conceptually expects); proper coordinate conversion
    is a future improvement.
    """
    if not cds_records:
        return None
    cds_start = min(rec.start for rec in cds_records)
    cds_end = max(rec.end for rec in cds_records)
    stable_id = _first_attr(cds_records[0], "protein_id") or ""
    return Translation(start=cds_start, end=cds_end, stable_id=stable_id)


def _build_gene(gene_rec: GFF3Record, transcripts: list[Transcript]) -> Gene:
    """Construct a Gene domain object from a gene record and link it to its transcripts."""
    strand = STRAND_SYMBOL_TO_ENUM.get(gene_rec.strand_symbol, Strand.UNKNOWN)
    raw_id = _first_attr(gene_rec, "ID") or ""
    stable_id = raw_id.split(":")[-1] if ":" in raw_id else raw_id

    gene = Gene(
        stable_id=stable_id,
        start=gene_rec.start,
        end=gene_rec.end,
        strand=strand,
        slice=Slice(
            coord_system=_COORD_SYSTEM,
            seq_region_name=gene_rec.seqid,
            start=gene_rec.start,
            end=gene_rec.end,
            strand=strand,
        ),
    )
    gene.name = _first_attr(gene_rec, "Name")
    gene.biotype = _first_attr(gene_rec, "biotype")
    gene.transcripts = transcripts

    for t in transcripts:
        t.gene = gene

    return gene


def _first_attr(rec: GFF3Record, key: str) -> str | None:
    """Return the first value of a GFF3 attribute, or None if absent."""
    values = rec.attributes.get(key)
    if values:
        return values[0]
    return None
