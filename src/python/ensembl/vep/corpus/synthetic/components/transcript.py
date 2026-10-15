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
"""Build synthetic transcript model objects."""

from __future__ import annotations

from collections.abc import Sequence

from ensembl.vep.common import Strand, reverse_complement
from ensembl.vep.models.genic import Exon, Gene, Intron, Transcript, Translation
from ensembl.vep.models.sequence import Slice

__all__ = [
    "build_cds_genomic_bounds",
    "build_transcript",
    "build_translation",
    "calculate_cds_length",
    "patch_terminal_codons",
]


def calculate_cds_length(
    exons_genomic_ascending: Sequence[tuple[int, int]],
    cds_genomic_low: int,
    cds_genomic_high: int,
) -> int:
    """Return coding bases within inclusive bounds across spliced exons.

    Args:
        exons_genomic_ascending (Sequence[tuple[int, int]]): Inclusive exon bounds in genomic order.
        cds_genomic_low (int): Inclusive lower genomic CDS boundary.
        cds_genomic_high (int): Inclusive upper genomic CDS boundary.

    Returns:
        int: Exonic bases intersecting the CDS bounds.
    """
    return sum(
        min(exon_end, cds_genomic_high) - max(exon_start, cds_genomic_low) + 1
        for exon_start, exon_end in exons_genomic_ascending
        if max(exon_start, cds_genomic_low) <= min(exon_end, cds_genomic_high)
    )


def _genomic_position_at_transcript_offset(
    exons_genomic_ascending: Sequence[tuple[int, int]],
    strand: Strand,
    transcript_offset: int,
) -> int:
    """Return the genomic position at a zero-based spliced-transcript offset."""
    exons_in_transcript_order = (
        exons_genomic_ascending if strand == Strand.FORWARD else tuple(reversed(exons_genomic_ascending))
    )
    remaining = transcript_offset
    for exon_start, exon_end in exons_in_transcript_order[:-1]:
        exon_length = exon_end - exon_start + 1
        if remaining < exon_length:
            return exon_start + remaining if strand == Strand.FORWARD else exon_end - remaining
        remaining -= exon_length
    exon_start, exon_end = exons_in_transcript_order[-1]
    return exon_start + remaining if strand == Strand.FORWARD else exon_end - remaining


def build_cds_genomic_bounds(
    exons_genomic_ascending: Sequence[tuple[int, int]],
    *,
    strand: Strand,
    frame_index: int,
    start_offset: int = 50,
    stop_offset: int = 50,
) -> tuple[int, int]:
    """Return CDS genomic bounds from spliced-transcript terminal offsets.

    ``start_offset`` is measured from the transcript's 5' end and carries into later exons when
    the first exon is shorter. ``stop_offset`` is measured back from the 3' end and is clamped
    towards that end when the requested position would cross the start codon. The terminal bound
    is then moved towards the 5' end until the CDS length is a multiple of three.

    Args:
        exons_genomic_ascending (Sequence[tuple[int, int]]): Exons in ascending genomic order.
        strand (Strand): Transcript strand.
        frame_index (int): Reading-frame index used to spread alternative starts across frames.
        start_offset (int): Preferred zero-based start offset from the transcript's 5' end.
        stop_offset (int): Preferred zero-based terminal offset back from the transcript's 3' end.

    Returns:
        tuple[int, int]: Inclusive genomic ``(low, high)`` CDS bounds.

    Raises:
        ValueError: If the exon layout or requested offsets cannot hold a codon.
    """
    if not exons_genomic_ascending:
        raise ValueError("a coding transcript needs at least one exon")
    if frame_index < 0 or start_offset < 0 or stop_offset < 0:
        raise ValueError("CDS frame and terminal offsets must not be negative")

    five_prime_exon = exons_genomic_ascending[0] if strand == Strand.FORWARD else exons_genomic_ascending[-1]
    five_prime_exon_length = five_prime_exon[1] - five_prime_exon[0] + 1
    adjusted_start_offset = start_offset + ((five_prime_exon_length - start_offset - frame_index) % 3)
    spliced_length = sum(exon_end - exon_start + 1 for exon_start, exon_end in exons_genomic_ascending)
    available_after_start = spliced_length - adjusted_start_offset
    if available_after_start < 3:
        raise ValueError("transcript is too short to place a complete CDS start codon")

    maximum_stop_offset = available_after_start - 3
    clamped_stop_offset = min(stop_offset, maximum_stop_offset)
    terminal_offset = spliced_length - clamped_stop_offset - 1
    coding_length = terminal_offset - adjusted_start_offset + 1
    terminal_offset -= coding_length % 3

    start_position = _genomic_position_at_transcript_offset(
        exons_genomic_ascending,
        strand,
        adjusted_start_offset,
    )
    terminal_position = _genomic_position_at_transcript_offset(
        exons_genomic_ascending,
        strand,
        terminal_offset,
    )
    cds_genomic_low = min(start_position, terminal_position)
    cds_genomic_high = max(start_position, terminal_position)
    return cds_genomic_low, cds_genomic_high


def patch_terminal_codons(chromosome_chars: list[str], transcripts: Sequence[Transcript]) -> None:
    """Write ATG start and TAA stop codons for every coding transcript.

    Args:
        chromosome_chars (list[str]): Mutable forward-strand chromosome sequence.
        transcripts (Sequence[Transcript]): Transcripts whose CDS terminals are patched.

    Raises:
        ValueError: If a terminal codon cannot be mapped inside its transcript CDS.
    """
    for transcript in transcripts:
        if transcript.translation is None:
            continue
        transcript_cds_length = transcript.cds_length
        assert transcript_cds_length is not None
        start_positions = [transcript.cds_offset_to_genomic(offset) for offset in range(3)]
        stop_positions = [
            transcript.cds_offset_to_genomic(transcript_cds_length - 3 + offset) for offset in range(3)
        ]
        for codon_bases, genomic_positions in (("ATG", start_positions), ("TAA", stop_positions)):
            if any(genomic_position is None for genomic_position in genomic_positions):
                raise ValueError(f"{transcript.stable_id} terminal codon falls outside its CDS")
            for base, genomic_position in zip(codon_bases, genomic_positions, strict=True):
                assert genomic_position is not None
                chromosome_chars[genomic_position - 1] = (
                    base if transcript.strand == Strand.FORWARD else reverse_complement(base)
                )


def build_translation(transcript: Transcript, genomic_bounds: tuple[int, int]) -> Translation:
    """Build a transcript translation from inclusive genomic CDS bounds.

    Args:
        transcript (Transcript): Transcript whose spliced exons contain both CDS boundaries.
        genomic_bounds (tuple[int, int]): Inclusive ``(low, high)`` genomic CDS bounds.

    Returns:
        Translation: Translation linked to ``transcript`` with exon-relative coding offsets.

    Raises:
        ValueError: If either CDS boundary falls outside the transcript's exons.
    """
    cds_genomic_low, cds_genomic_high = genomic_bounds
    coding_start_position = cds_genomic_low if transcript.strand == Strand.FORWARD else cds_genomic_high
    coding_end_position = cds_genomic_high if transcript.strand == Strand.FORWARD else cds_genomic_low
    coding_start_exon = next(
        (exon for exon in transcript.spliced_exons if exon.start <= coding_start_position <= exon.end),
        None,
    )
    coding_end_exon = next(
        (exon for exon in transcript.spliced_exons if exon.start <= coding_end_position <= exon.end),
        None,
    )
    if coding_start_exon is None or coding_end_exon is None:
        raise ValueError("CDS bounds must fall within transcript exons")

    if transcript.strand == Strand.FORWARD:
        exon_coding_start = coding_start_position - coding_start_exon.start + 1
        exon_coding_end = coding_end_position - coding_end_exon.start + 1
        translation_start = cds_genomic_low - transcript.start + 1
        translation_end = cds_genomic_high - transcript.start + 1
    else:
        exon_coding_start = coding_start_exon.end - coding_start_position + 1
        exon_coding_end = coding_end_exon.end - coding_end_position + 1
        translation_start = transcript.end - cds_genomic_high + 1
        translation_end = transcript.end - cds_genomic_low + 1

    translation = Translation(
        start=translation_start,
        end=translation_end,
        stable_id=f"{transcript.stable_id}_translation",
        exon_coding_start=exon_coding_start,
        exon_coding_end=exon_coding_end,
    )
    translation.transcript = transcript
    return translation


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
    """Build one transcript with exons, introns, and an optional translation.

    Args:
        stable_id (str): Transcript stable identifier.
        synthetic_slice (Slice): Sequence region shared by the transcript and its exons.
        strand (Strand): Transcript orientation.
        exons_genomic_ascending (Sequence[tuple[int, int]]): Inclusive exon bounds in genomic order.
        cds_genomic_bounds (tuple[int, int] | None): Inclusive CDS bounds, or ``None`` for a
            non-coding transcript.
        gene (Gene | None): Optional parent gene.
        biotype (str | None): Optional transcript biotype.

    Returns:
        Transcript: Constructed transcript with seeded introns and translation.

    Raises:
        ValueError: If a CDS boundary does not fall inside a transcript exon.
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
    introns: list[Intron] = []
    for prev_exon, next_exon in zip(spliced_exons, spliced_exons[1:]):
        first, second = (next_exon, prev_exon) if strand == Strand.REVERSE else (prev_exon, next_exon)
        intron = Intron(start=first.end + 1, end=second.start - 1, strand=strand, slice=synthetic_slice)
        intron.assign_without_hooks(prev_exon=prev_exon, next_exon=next_exon)
        introns.append(intron)
    transcript = Transcript(
        start=exons_genomic_ascending[0][0],
        end=exons_genomic_ascending[-1][1],
        strand=strand,
        stable_id=stable_id,
        slice=synthetic_slice,
    )
    transcript.spliced_exons = spliced_exons
    transcript.seed_attrs(introns=introns)
    transcript.gene = gene
    transcript.biotype = biotype

    if cds_genomic_bounds is not None:
        transcript.seed_attrs(translation=build_translation(transcript, cds_genomic_bounds))
    return transcript
