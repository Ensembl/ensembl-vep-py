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
"""GFF3 rendering for corpora."""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from urllib.parse import quote

from ensembl.vep.common import Strand
from ensembl.vep.corpus.base import Corpus
from ensembl.vep.models.genic import Gene, Transcript

__all__ = ["GFF3_NAME", "escape_attribute_value", "gff3_lines"]

GFF3_NAME = "corpus.gff3"

_SOURCE = "ensembl_vep_synthetic"
_ATTRIBUTE_SAFE_CHARACTERS = ".:^*$@!+_?-|"


def escape_attribute_value(value: str) -> str:
    """Percent-encode one scalar GFF3 attribute value."""
    return quote(value, safe=_ATTRIBUTE_SAFE_CHARACTERS)


def _typed_id(kind: str, stable_id: str) -> str:
    """Return a GFF3 identifier without duplicating an existing type prefix."""
    prefix = f"{kind}:"
    return stable_id if stable_id.startswith(prefix) else f"{prefix}{stable_id}"


def _display_name(kind: str, stable_id: str) -> str:
    """Return the untyped portion of a conventionally prefixed GFF3 identifier."""
    prefix = f"{kind}:"
    return stable_id.removeprefix(prefix)


def _translation_id(transcript: Transcript) -> str:
    """Return a stable, non-empty GFF3 identifier for a transcript's translation."""
    translation = transcript.translation
    assert translation is not None
    stable_id = translation.stable_id or f"{_display_name('transcript', transcript.stable_id)}_translation"
    return _typed_id("CDS", stable_id)


def _strand_symbol(strand: Strand) -> str:
    """Return the GFF3 strand column for a strand."""
    return "+" if strand == Strand.FORWARD else "-"


def _row(
    seq_region_name: str,
    feature_type: str,
    start: int,
    end: int,
    strand: Strand,
    attributes: dict[str, str],
    phase: int | None = None,
) -> str:
    """Return one nine-column GFF3 line."""
    attribute_text = ";".join(f"{key}={escape_attribute_value(value)}" for key, value in attributes.items())
    phase_text = "." if phase is None else str(phase)
    return "\t".join(
        (
            seq_region_name,
            _SOURCE,
            feature_type,
            str(start),
            str(end),
            ".",
            _strand_symbol(strand),
            phase_text,
            attribute_text,
        )
    )


def _cds_rows(transcript: Transcript, seq_region_name: str) -> list[tuple[int, str]]:
    """Return the CDS rows, each keyed by start, with the phase each block begins in."""
    spans = transcript.cds_genomic_spans
    translation = transcript.translation
    if not spans or translation is None:
        return []
    rows = []
    coding_bases_before = 0
    for block_start, block_end in spans:
        phase = (3 - coding_bases_before % 3) % 3
        rows.append(
            (
                block_start,
                _row(
                    seq_region_name,
                    "CDS",
                    block_start,
                    block_end,
                    transcript.strand,
                    {
                        "ID": _translation_id(transcript),
                        "Parent": _typed_id("transcript", transcript.stable_id),
                        "protein_id": _translation_id(transcript),
                    },
                    phase=phase,
                ),
            )
        )
        coding_bases_before += block_end - block_start + 1
    return rows


def _utr_rows(transcript: Transcript, seq_region_name: str) -> list[tuple[int, str]]:
    """Return the UTR rows, each keyed by start."""
    rows = []
    for feature_type, spans in (
        ("five_prime_UTR", transcript.five_prime_utr_genomic_spans),
        ("three_prime_UTR", transcript.three_prime_utr_genomic_spans),
    ):
        for block_start, block_end in spans or []:
            rows.append(
                (
                    block_start,
                    _row(
                        seq_region_name,
                        feature_type,
                        block_start,
                        block_end,
                        transcript.strand,
                        {"Parent": _typed_id("transcript", transcript.stable_id)},
                    ),
                )
            )
    return rows


def _transcript_rows(transcript: Transcript, seq_region_name: str) -> list[str]:
    """Return a transcript's row followed by its exon, CDS and UTR rows in genomic order."""
    attributes = {
        "ID": _typed_id("transcript", transcript.stable_id),
        "Name": _display_name("transcript", transcript.stable_id),
    }
    if transcript.gene is not None:
        attributes["Parent"] = _typed_id("gene", transcript.gene.stable_id)
    biotype = transcript.biotype
    if biotype is None and transcript.gene is not None:
        biotype = transcript.gene.biotype
    if biotype is not None:
        attributes["biotype"] = biotype
    rows = [
        _row(seq_region_name, "transcript", transcript.start, transcript.end, transcript.strand, attributes)
    ]

    children: list[tuple[int, str]] = []
    for exon in transcript.spliced_exons:
        children.append(
            (
                exon.start,
                _row(
                    seq_region_name,
                    "exon",
                    exon.start,
                    exon.end,
                    exon.strand,
                    {
                        "ID": _typed_id("exon", exon.stable_id),
                        "Parent": _typed_id("transcript", transcript.stable_id),
                    },
                ),
            )
        )
    children.extend(_cds_rows(transcript, seq_region_name))
    children.extend(_utr_rows(transcript, seq_region_name))
    children.sort(key=lambda entry: entry[0])
    rows.extend(row for _start, row in children)
    return rows


def _gene_row(gene: Gene, seq_region_name: str) -> str:
    """Return a gene's row."""
    attributes = {"ID": _typed_id("gene", gene.stable_id), "Name": _display_name("gene", gene.stable_id)}
    if gene.biotype is not None:
        attributes["biotype"] = gene.biotype
    return _row(seq_region_name, "gene", gene.start, gene.end, gene.strand, attributes)


def gff3_lines(corpora: Sequence[Corpus]) -> Iterator[str]:
    """Yield the GFF3 records for corpora."""
    yield "##gff-version 3"
    for corpus in corpora:
        for sequence_region in corpus.sequence_regions:
            yield f"##sequence-region {sequence_region.seq_region_name} 1 {sequence_region.end}"
    for corpus in corpora:
        written_genes: set[str] = set()
        for transcript in corpus.transcripts:
            if transcript.slice is None:
                continue
            seq_region_name = transcript.slice.seq_region_name
            gene = transcript.gene
            if gene is not None and gene.stable_id not in written_genes:
                written_genes.add(gene.stable_id)
                yield _gene_row(gene, seq_region_name)
            yield from _transcript_rows(transcript, seq_region_name)
