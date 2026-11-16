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
"""Build synthetic sequences and transcripts."""

from __future__ import annotations

import random
from collections.abc import Sequence

from ensembl.vep.common import Strand, reverse_complement
from ensembl.vep.models.genic import Exon, Gene, Intron, Transcript, Translation
from ensembl.vep.models.sequence import Slice

__all__ = [
    "build_transcript",
]

_STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
_BASES = "ACGT"


def _random_sequence(length: int, random_generator: random.Random) -> str:
    """Returns a random nucleotide sequence.

    Args:
        length (int): Number of bases to generate.
        random_generator (random.Random): Seeded random generator, for reproducibility.

    Returns:
        str: A random sequence of ``ACGT`` bases.
    """
    return "".join(random_generator.choice(_BASES) for _ in range(length))


def _make_clean_open_reading_frame(codon_count: int, random_generator: random.Random) -> str:
    """Build a clean ORF: ``ATG``, no internal stop codons, and a terminal stop codon.

    Args:
        codon_count (int): Total number of codons, including the start and stop codon.
        random_generator (random.Random): Seeded random generator, for reproducibility.

    Returns:
        str: The coding sequence (transcript-sense), ``codon_count * 3`` bases long.
    """
    codons = ["ATG"]
    while len(codons) < codon_count - 1:
        candidate_codon = "".join(random_generator.choice(_BASES) for _ in range(3))
        if candidate_codon not in _STOP_CODONS:
            codons.append(candidate_codon)
    codons.append("TAA")
    return "".join(codons)


def _translation_offsets_for_genomic_cds_bounds(
    exons_genomic_ascending: Sequence[tuple[int, int]],
    strand: Strand,
    cds_genomic_low: int,
    cds_genomic_high: int,
) -> tuple[int, int]:
    """Derive Translation's ``exon_coding_start``/``exon_coding_end`` from CDS bounds.

    Inverts Transcript._cds_genomic_spans's clipping arithmetic: given the desired
    outermost genomic bounds of the CDS and the transcript's exon layout, returns the
    exon-relative offsets that will make ``Transcript.cds_genomic_spans`` reproduce exactly
    ``(cds_genomic_low, cds_genomic_high)`` as its envelope.

    Args:
        exons_genomic_ascending (Sequence[tuple[int, int]]): Every exon of the transcript, in
            ascending genomic order (regardless of strand).
        strand (Strand): The transcript's strand.
        cds_genomic_low (int): Lowest genomic position of the CDS.
        cds_genomic_high (int): Highest genomic position of the CDS.

    Returns:
        tuple[int, int]: ``(exon_coding_start, exon_coding_end)``.
    """
    if strand == Strand.FORWARD:
        transcript_first_exon = exons_genomic_ascending[0]
        transcript_last_exon = exons_genomic_ascending[-1]
        exon_coding_start = cds_genomic_low - transcript_first_exon[0] + 1
        exon_coding_end = cds_genomic_high - transcript_last_exon[0] + 1
    else:
        transcript_first_exon = exons_genomic_ascending[-1]
        transcript_last_exon = exons_genomic_ascending[0]
        exon_coding_start = transcript_first_exon[1] - cds_genomic_high + 1
        exon_coding_end = transcript_last_exon[1] - cds_genomic_low + 1
    return exon_coding_start, exon_coding_end


def build_transcript(
    *,
    stable_id: str,
    synthetic_slice: Slice,
    strand: Strand,
    exons_genomic_ascending: Sequence[tuple[int, int]],
    cds_genomic_bounds: tuple[int, int] | None,
    gene: Gene | None,
    biotype: str | None = None,
) -> Transcript:
    """Build one Transcript (with exons, introns, and optionally a translation) over a slice.

    Args:
        stable_id (str): Stable identifier for the transcript.
        synthetic_slice (Slice): The shared synthetic-chromosome slice.
        strand (Strand): The transcript's strand.
        exons_genomic_ascending (Sequence[tuple[int, int]]): Exon (start, end) genomic spans, in
            ascending genomic order.
        cds_genomic_bounds (tuple[int, int] | None): (low, high) genomic bounds of the CDS, or
            None for a non-coding transcript.
        gene (Gene | None): The owning gene, or None.
        biotype (str | None): The transcript's own biotype, or None to leave it undeclared so
            that ``FactContext.biotype`` falls back to the gene.

    Returns:
        Transcript: The constructed transcript, with ``spliced_exons`` in transcript (5' to 3')
            order and ``introns`` in ascending genomic order.
    """
    exons_transcript_order = (
        list(exons_genomic_ascending) if strand == Strand.FORWARD else list(reversed(exons_genomic_ascending))
    )
    spliced_exons = [
        Exon(
            start=exon_start,
            end=exon_end,
            strand=strand,
            stable_id=f"{stable_id}_exon{exon_index}",
            slice=synthetic_slice,
        )
        for exon_index, (exon_start, exon_end) in enumerate(exons_transcript_order)
    ]

    introns = [
        Intron(
            start=exons_genomic_ascending[index][1] + 1,
            end=exons_genomic_ascending[index + 1][0] - 1,
            strand=strand,
            slice=synthetic_slice,
        )
        for index in range(len(exons_genomic_ascending) - 1)
        if exons_genomic_ascending[index + 1][0] - exons_genomic_ascending[index][1] > 1
    ]

    transcript = Transcript(
        start=exons_genomic_ascending[0][0],
        end=exons_genomic_ascending[-1][1],
        strand=strand,
        stable_id=stable_id,
        slice=synthetic_slice,
    )
    transcript.spliced_exons = spliced_exons
    # Seeded after spliced_exons: that assignment's on_setattr hook invalidates _introns, which
    # would otherwise discard this seed immediately. introns is seeded rather than derived because
    # the derivation records a zero-gap adjacent-exon pair as an intron; this generator's own gap
    # list skips those.
    transcript.seed_attrs(introns=introns)
    transcript.gene = gene
    transcript.biotype = biotype

    if cds_genomic_bounds is not None:
        cds_genomic_low, cds_genomic_high = cds_genomic_bounds
        exon_coding_start, exon_coding_end = _translation_offsets_for_genomic_cds_bounds(
            exons_genomic_ascending, strand, cds_genomic_low, cds_genomic_high
        )
        # Transcript-relative and **1-based**, measured from this transcript's own 5' end --
        # self.start on the forward strand, self.end on the reverse strand -- the exact convention
        # Transcript._coding_start/_coding_end's translation-fallback branches expect, and the
        # inverse of Transcript._translation's own (derived-translation) formula. 1-based means
        # translation_start == 1 is the transcript's first base, matching Feature.start/end and
        # exon_coding_start; the +1 is what makes that round-trip.
        if strand == Strand.FORWARD:
            translation_start = cds_genomic_low - transcript.start + 1
            translation_end = cds_genomic_high - transcript.start + 1
        else:
            translation_start = transcript.end - cds_genomic_high + 1
            translation_end = transcript.end - cds_genomic_low + 1
        translation = Translation(
            start=translation_start,
            end=translation_end,
            stable_id=f"{stable_id}_translation",
            exon_coding_start=exon_coding_start,
            exon_coding_end=exon_coding_end,
        )
        translation.transcript = transcript
        transcript.seed_attrs(translation=translation)

    return transcript


def _write_transcript_cds(chromosome_chars: list[str], transcript: Transcript, coding_sequence: str) -> None:
    """Write a transcript-sense coding sequence into the raw chromosome character list.

    ``cds_genomic_spans`` is transcript order (not necessarily genomic-ascending -- see its own
    docstring), so this writes per block: take the next chunk of ``coding_sequence`` (transcript
    order matches block order), reverse-complement *that chunk alone* if the transcript is on the
    reverse strand, then write it ascending into that block's own genomic range -- the inverse of
    ``Transcript._coding_sequence``'s per-block extraction, not a whole-string reverse-complement
    (which only works when blocks are walked in genomic-ascending order).

    Args:
        chromosome_chars (list[str]): The whole synthetic chromosome, as a mutable list of
            single-character strings, modified in place.
        transcript (Transcript): The (already-linked) transcript to write the CDS under.
        coding_sequence (str): Transcript-sense (5' to 3') coding sequence, exactly as long as
            ``transcript.cds_length``.
    """
    genomic_spans = transcript.cds_genomic_spans
    assert genomic_spans is not None
    read_offset = 0
    for block_start, block_end in genomic_spans:
        block_length = block_end - block_start + 1
        chunk = coding_sequence[read_offset : read_offset + block_length]
        read_offset += block_length
        forward_sense_chunk = chunk if transcript.strand == Strand.FORWARD else reverse_complement(chunk)
        for position, base in zip(range(block_start, block_end + 1), forward_sense_chunk):
            chromosome_chars[position - 1] = base
    assert read_offset == len(coding_sequence)
