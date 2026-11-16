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
"""Build the dense heavy-locus synthetic corpus."""

from __future__ import annotations

import random
from collections.abc import Sequence

from ensembl.vep.common import CoordSystemLevel, Strand, invalidate, reverse_complement
from ensembl.vep.models.genic import Gene, Transcript
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.models.variation import VariationAllele

from ensembl.vep.corpus.registry import register_corpus
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.corpus.synthetic.base import SyntheticCorpus
from ensembl.vep.corpus.synthetic.components.sequence import _random_sequence, build_transcript
from ensembl.vep.corpus.synthetic.components.variation import (
    _deletion,
    _insertion,
    _snv,
    build_transcript_variation_alleles,
)

__all__ = [
    "build_heavy_locus_genome",
    "build_heavy_locus_transcripts",
    "heavy_locus_corpus",
    "heavy_locus_variant_alleles",
]

# --------------------------------------------------------------------------
# A heavy transcript locus.
#
# The four transcripts above are disjoint, which makes them useless for measuring anything that
# shares work *between* transcripts. A real dense locus has many transcripts over one gene
# sharing most of their internal exons and differing at the edges, so this builds that shape (2
# strands x 18 exon-sharing transcripts per strand), porting ``socore.data.build_heavy_locus``:
#
# * three alternative first exons (alternative TSS), of different lengths, so transcripts enter
#   the shared exons in different reading frames;
# * five constitutive internal exons present in every transcript;
# * one cassette exon of 61bp -- deliberately not a multiple of three -- so including it shifts
#   the frame of everything downstream of it;
# * three alternative last exons.
#
# 3 x 2 x 3 = 18 transcripts per locus, one locus per strand. This is the only corpus that
# exercises multi-variant-per-transcript zone reuse (the zone map cache, caches/zone_map.py)
# and, in a later phase, cross-transcript codon sharing.
# --------------------------------------------------------------------------

_HEAVY_LOCUS_FIRST_EXONS: tuple[tuple[int, int], ...] = ((0, 200), (1000, 1200), (2000, 2200))
_HEAVY_LOCUS_CONSTITUTIVE_EXONS_UPSTREAM_OF_CASSETTE: tuple[tuple[int, int], ...] = (
    (3000, 3150),
    (4000, 4120),
)
_HEAVY_LOCUS_CASSETTE_EXON: tuple[int, int] = (4500, 4560)
_HEAVY_LOCUS_CONSTITUTIVE_EXONS_DOWNSTREAM_OF_CASSETTE: tuple[tuple[int, int], ...] = (
    (5000, 5180),
    (6000, 6090),
    (7000, 7210),
)
_HEAVY_LOCUS_LAST_EXONS: tuple[tuple[int, int], ...] = ((8000, 8300), (9000, 9400), (10000, 10300))
_HEAVY_LOCUS_SPAN = 10300
_HEAVY_LOCUS_CHROMOSOME_LENGTH = 130000
_HEAVY_LOCUS_CHROMOSOME_NAME = "chrHeavyLocus"
_HEAVY_LOCUS_FORWARD_ORIGIN = 10000
_HEAVY_LOCUS_REVERSE_ORIGIN = 70000


def _heavy_locus_cds_length(
    exons_genomic_ascending: Sequence[tuple[int, int]], cds_genomic_low: int, cds_genomic_high: int
) -> int:
    """Returns the number of coding bases across ``exons_genomic_ascending`` clipped to CDS bounds.

    Args:
        exons_genomic_ascending (Sequence[tuple[int, int]]): Exon (start, end) genomic spans, in
            ascending genomic order.
        cds_genomic_low (int): Candidate lowest genomic bound of the CDS.
        cds_genomic_high (int): Candidate highest genomic bound of the CDS.

    Returns:
        int: The number of coding bases, i.e. what ``Transcript.cds_length`` would report for
            these bounds.
    """
    return sum(
        min(block_end, cds_genomic_high) - max(block_start, cds_genomic_low) + 1
        for block_start, block_end in exons_genomic_ascending
        if max(block_start, cds_genomic_low) <= min(block_end, cds_genomic_high)
    )


def _heavy_locus_biotype(transcript_index: int) -> tuple[str, bool]:
    """Returns ``(biotype, is_coding)`` for a heavy-locus transcript, by its index in the locus.

    Args:
        transcript_index (int): The transcript's 0-based index within its 18-transcript locus.

    Returns:
        tuple[str, bool]: The transcript's biotype and whether it is coding.
    """
    if transcript_index % 11 == 5:
        return "lincRNA", False
    if transcript_index % 7 == 3:
        return "nonsense_mediated_decay", True
    return "protein_coding", True


def _heavy_locus_cds_genomic_bounds(
    exons_genomic_ascending: Sequence[tuple[int, int]],
    *,
    strand: Strand,
    frame_index: int,
    transcript_index: int,
) -> tuple[int, int]:
    """Returns the CDS genomic bounds for one heavy-locus transcript.

    Args:
        exons_genomic_ascending (Sequence[tuple[int, int]]): The transcript's exon blocks,
            ascending
            genomic order.
        strand (Strand): The transcript's strand.
        frame_index (int): Index of the transcript's 5' exon among its alternatives (see
            build_heavy_locus_transcripts), used to spread alternative TSSs across the
            three reading frames.
        transcript_index (int): The transcript's 0-based index within its locus, used to space
            distinct stop codons for transcripts sharing a 3' exon.

    Returns:
        tuple[int, int]: ``(cds_genomic_low, cds_genomic_high)``.
    """
    # The start codon sits 50 + frame_index bases into the 5' exon. The +frame_index is what
    # makes the three alternative TSSs enter the shared exons in three *different* reading
    # frames -- the interesting case for a zone map, since the same genomic footprint is
    # "overlap_cds" for every transcript sharing it, but a frame-based (rather than
    # genomic-footprint-based) fact would get the codon boundaries wrong for two of the three.
    five_prime_exon = exons_genomic_ascending[0] if strand == Strand.FORWARD else exons_genomic_ascending[-1]
    five_prime_exon_length = five_prime_exon[1] - five_prime_exon[0] + 1
    start_offset = 50 + ((five_prime_exon_length - 50 - frame_index) % 3)
    # The stop is spaced 6bp per transcript so transcripts sharing a 3' exon get distinct,
    # non-overlapping stop codons instead of overwriting each other's.
    stop_offset = 50 + 6 * transcript_index
    if strand == Strand.FORWARD:
        cds_genomic_low = exons_genomic_ascending[0][0] + start_offset
        cds_genomic_high = exons_genomic_ascending[-1][1] - stop_offset
    else:
        # On the reverse strand the low genomic coordinate is the 3' end.
        cds_genomic_low = exons_genomic_ascending[0][0] + stop_offset
        cds_genomic_high = exons_genomic_ascending[-1][1] - start_offset
    while _heavy_locus_cds_length(exons_genomic_ascending, cds_genomic_low, cds_genomic_high) % 3:
        if strand == Strand.FORWARD:
            cds_genomic_high -= 1
        else:
            cds_genomic_low += 1
    return cds_genomic_low, cds_genomic_high


def _build_one_heavy_locus_transcript(
    synthetic_slice: Slice,
    *,
    strand: Strand,
    stable_id_tag: str,
    transcript_index: int,
    first_exon_index: int,
    last_exon_index: int,
    exons_genomic_ascending: Sequence[tuple[int, int]],
) -> Transcript:
    """Build a single heavy-locus transcript over a fixed exon layout.

    Args:
        synthetic_slice (Slice): The shared synthetic-chromosome slice.
        strand (Strand): The transcript's strand.
        stable_id_tag (str): Prefix for the transcript's stable_id.
        transcript_index (int): The transcript's 0-based index within its locus.
        first_exon_index (int): Index of the alternative first exon this transcript uses.
        last_exon_index (int): Index of the alternative last exon this transcript uses.
        exons_genomic_ascending (Sequence[tuple[int, int]]): This transcript's exon blocks,
            ascending
            genomic order.

    Returns:
        Transcript: The constructed transcript, not yet carrying real start/stop codon sequence.
    """
    biotype, is_coding = _heavy_locus_biotype(transcript_index)

    cds_genomic_bounds: tuple[int, int] | None = None
    if is_coding:
        frame_index = first_exon_index if strand == Strand.FORWARD else last_exon_index
        cds_genomic_bounds = _heavy_locus_cds_genomic_bounds(
            exons_genomic_ascending, strand=strand, frame_index=frame_index, transcript_index=transcript_index
        )

    stable_id = f"{stable_id_tag}{transcript_index:02d}"
    gene: Gene | None = None
    if biotype == "nonsense_mediated_decay":
        gene = Gene(
            start=exons_genomic_ascending[0][0],
            end=exons_genomic_ascending[-1][1],
            strand=strand,
            stable_id=f"{stable_id}_gene",
            slice=synthetic_slice,
        )
        gene.biotype = biotype

    return build_transcript(
        stable_id=stable_id,
        synthetic_slice=synthetic_slice,
        strand=strand,
        exons_genomic_ascending=exons_genomic_ascending,
        cds_genomic_bounds=cds_genomic_bounds,
        gene=gene,
    )


def build_heavy_locus_transcripts(
    synthetic_slice: Slice, *, origin: int, strand: Strand, stable_id_tag: str
) -> list[Transcript]:
    """Build 18 transcripts sharing constitutive exons over a single locus.

    Mirrors ``socore.data.build_heavy_locus``.

    Args:
        synthetic_slice (Slice): The shared synthetic-chromosome slice.
        origin (int): Genomic (slice-relative) offset the locus's exon coordinates are relative
            to.
        strand (Strand): The strand every transcript in this locus is built on.
        stable_id_tag (str): Prefix for each transcript's stable_id (e.g. "FWD", "REV").

    Returns:
        list[Transcript]: The 18 transcripts, not yet carrying real start/stop codon sequence
            (see build_heavy_locus_genome, which patches that in afterwards).
    """
    transcripts: list[Transcript] = []
    transcript_index = 0
    for first_exon_index, first_exon in enumerate(_HEAVY_LOCUS_FIRST_EXONS):
        for include_cassette_exon in (False, True):
            for last_exon_index, last_exon in enumerate(_HEAVY_LOCUS_LAST_EXONS):
                exon_blocks = (
                    [first_exon]
                    + list(_HEAVY_LOCUS_CONSTITUTIVE_EXONS_UPSTREAM_OF_CASSETTE)
                    + ([_HEAVY_LOCUS_CASSETTE_EXON] if include_cassette_exon else [])
                    + list(_HEAVY_LOCUS_CONSTITUTIVE_EXONS_DOWNSTREAM_OF_CASSETTE)
                    + [last_exon]
                )
                exons_genomic_ascending = sorted(
                    (origin + block_start, origin + block_end) for block_start, block_end in exon_blocks
                )
                transcripts.append(
                    _build_one_heavy_locus_transcript(
                        synthetic_slice,
                        strand=strand,
                        stable_id_tag=stable_id_tag,
                        transcript_index=transcript_index,
                        first_exon_index=first_exon_index,
                        last_exon_index=last_exon_index,
                        exons_genomic_ascending=exons_genomic_ascending,
                    )
                )
                transcript_index += 1
    return transcripts


def _patch_heavy_locus_start_stop_codons(
    chromosome_chars: list[str], transcripts: Sequence[Transcript]
) -> None:
    """Write real ATG/TAA start/stop codons for every coding transcript in ``transcripts``.

    Args:
        chromosome_chars (list[str]): The whole synthetic chromosome, as a mutable list of
            single-character strings, modified in place.
        transcripts (Sequence[Transcript]): Transcripts to patch codons into (only coding ones are
            touched).
    """
    for transcript in transcripts:
        if transcript.translation is None:
            continue
        cds_length = transcript.cds_length
        assert cds_length is not None
        start_codon_positions = [transcript.cds_offset_to_genomic(offset) for offset in range(3)]
        stop_codon_positions = [
            transcript.cds_offset_to_genomic(cds_length - 3 + offset) for offset in range(3)
        ]
        for codon_bases, genomic_positions in (("ATG", start_codon_positions), ("TAA", stop_codon_positions)):
            if any(genomic_position is None for genomic_position in genomic_positions):
                continue
            for base, genomic_position in zip(codon_bases, genomic_positions):
                assert genomic_position is not None
                chromosome_chars[genomic_position - 1] = (
                    base if transcript.strand == Strand.FORWARD else reverse_complement(base)
                )


def _fill_heavy_locus_background(
    chromosome_chars: list[str], origin: int, random_generator: random.Random
) -> None:
    """Fill one locus's span with a stop-codon-free ``{A, C, G}`` background.

    A sequence drawn only from ``{A, C, G}`` cannot contain ``TAA``, ``TAG`` or ``TGA`` in any
    reading frame, nor their reverse complements ``TTA``, ``CTA``, ``TCA`` -- so no transcript
    sharing this background picks up a spurious internal stop regardless of where its CDS starts
    or which of the three reading frames it enters the shared exons in. Start and stop codons are
    then patched in explicitly per transcript by _patch_heavy_locus_start_stop_codons.
    Mirrors ``socore.data._write_locus_sequence``.

    Args:
        chromosome_chars (list[str]): The whole synthetic chromosome, as a mutable list of
            single-character strings, modified in place.
        origin (int): Genomic (slice-relative) start of the locus span to fill.
        random_generator (random.Random): Seeded random generator, for reproducibility.
    """
    for position in range(origin, origin + _HEAVY_LOCUS_SPAN + 1):
        chromosome_chars[position - 1] = random_generator.choice("ACG")


def build_heavy_locus_genome(seed: int = 31) -> tuple[Slice, list[Transcript]]:
    """Build the two-locus (one per strand), 36-transcript heavy-locus synthetic chromosome.

    Mirrors ``socore.data.build_locus_genome``.

    Args:
        seed (int): Random seed, for reproducibility.

    Returns:
        tuple[Slice, list[Transcript]]: The shared chromosome slice (with real, codon-patched
            sequence) and all 36 transcripts (18 forward-strand, 18 reverse-strand).
    """
    random_generator = random.Random(seed)
    chromosome_chars = list(_random_sequence(_HEAVY_LOCUS_CHROMOSOME_LENGTH, random_generator))

    synthetic_slice = Slice(
        coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
        seq_region_name=_HEAVY_LOCUS_CHROMOSOME_NAME,
        start=1,
        end=_HEAVY_LOCUS_CHROMOSOME_LENGTH,
        strand=Strand.FORWARD,
    )
    with synthetic_slice.unsafe() as unsafe_slice:
        setattr(unsafe_slice, "_seq", "".join(chromosome_chars))

    forward_transcripts = build_heavy_locus_transcripts(
        synthetic_slice, origin=_HEAVY_LOCUS_FORWARD_ORIGIN, strand=Strand.FORWARD, stable_id_tag="FWD"
    )
    reverse_transcripts = build_heavy_locus_transcripts(
        synthetic_slice, origin=_HEAVY_LOCUS_REVERSE_ORIGIN, strand=Strand.REVERSE, stable_id_tag="REV"
    )
    transcripts = forward_transcripts + reverse_transcripts

    _fill_heavy_locus_background(chromosome_chars, _HEAVY_LOCUS_FORWARD_ORIGIN, random_generator)
    _fill_heavy_locus_background(chromosome_chars, _HEAVY_LOCUS_REVERSE_ORIGIN, random_generator)
    _patch_heavy_locus_start_stop_codons(chromosome_chars, transcripts)

    with synthetic_slice.unsafe() as unsafe_slice:
        setattr(unsafe_slice, "_seq", "".join(chromosome_chars))

    # Unlike build_synthetic_genome's four disjoint transcripts, transcripts here legitimately
    # share exons in different reading frames (that is the whole point of this corpus -- see the
    # module-level comment on frame diversity above), so one transcript's real start/stop codon
    # can land as a premature or missing stop in a sibling's frame. That is expected, not a data
    # bug -- mirrors ``socore.data.build_locus_genome``, which asserts nothing about translated
    # protein shape for this dataset (unlike ``build_genome``, the four-transcript corpus, which
    # does assert a clean ORF). Only the coding-sequence cache needs invalidating here.
    for transcript in transcripts:
        if transcript.translation is not None:
            invalidate(transcript, "_coding_sequence")

    return synthetic_slice, transcripts


def heavy_locus_variant_alleles(
    synthetic_slice: Slice, step: int = 13, seed: int = 41
) -> list[tuple[str, VariationAllele]]:
    """Build a dense sweep of variants across each locus, so shared exons are hit repeatedly.

    Mirrors ``socore.data.locus_variants``.

    Args:
        synthetic_slice (Slice): The chromosome slice built by build_heavy_locus_genome.
        step (int): Genomic spacing between consecutive sweep positions.
        seed (int): Random seed, for reproducibility.

    Returns:
        list[tuple[str, VariationAllele]]: Named heavy-locus sweep variants.

    Raises:
        ValueError: If ``step`` is not positive.
    """
    if step <= 0:
        raise ValueError(f"step must be positive, got {step}")
    random_generator = random.Random(seed)
    named_variants: list[tuple[str, VariationAllele]] = []
    for origin, tag in ((_HEAVY_LOCUS_FORWARD_ORIGIN, "FWD"), (_HEAVY_LOCUS_REVERSE_ORIGIN, "REV")):
        for sweep_index, position in enumerate(range(origin, origin + _HEAVY_LOCUS_SPAN, step)):
            named_variants.append(
                _snv(synthetic_slice, position, f"{tag}_snv{sweep_index}", random_generator)
            )
            if sweep_index % 3 == 0:
                named_variants.append(_deletion(synthetic_slice, position, 1, f"{tag}_del1_{sweep_index}"))
            if sweep_index % 5 == 0:
                named_variants.append(_deletion(synthetic_slice, position, 3, f"{tag}_del3_{sweep_index}"))
            if sweep_index % 7 == 0:
                named_variants.append(
                    _insertion(synthetic_slice, position, "ACG", f"{tag}_ins3_{sweep_index}")
                )
            if sweep_index % 11 == 0:
                named_variants.append(_insertion(synthetic_slice, position, "A", f"{tag}_ins1_{sweep_index}"))
    return named_variants


def heavy_locus_corpus(
    step: int = 13,
    seed: int = 31,
    *,
    sweep_seed: int = 41,
    spec: CorpusSpec | None = None,
) -> SyntheticCorpus:
    """Build the heavy-locus corpus: 2 strands x 18 exon-sharing transcripts, a dense variant sweep.

    Mirrors ``socore.data.heavy_locus_dataset``. Unlike default_corpus's four disjoint
    transcripts, this is the only corpus that exercises multi-variant-per-transcript zone reuse
    (the zone map cache, caches/zone_map.py) and, in a later phase, cross-transcript codon
    sharing -- many transcripts here share the same genomic exon footprint, so a transcript's
    zone map (or codon cache entry) built for one variant is reused by every other variant that
    lands in the same shared exon, across every transcript sharing it.

    Args:
        step (int): Genomic spacing between consecutive sweep positions.
        seed (int): Random seed for the underlying chromosome, for reproducibility.
        sweep_seed (int): Random seed for the variant sweep's alternate bases.

    Returns:
        SyntheticCorpus: The synthetic slice, its 36 transcripts, variants, and overlaps.
    """
    synthetic_slice, transcripts = build_heavy_locus_genome(seed=seed)
    named_variants = heavy_locus_variant_alleles(synthetic_slice, step=step, seed=sweep_seed)
    named_overlaps = build_transcript_variation_alleles(named_variants, transcripts)
    return SyntheticCorpus(
        synthetic_slice=synthetic_slice,
        transcripts=tuple(transcripts),
        named_variants=tuple(named_variants),
        named_overlaps=tuple(named_overlaps),
        spec=spec or CorpusSpec(kind="synthetic.heavy"),
    )


@register_corpus("synthetic.heavy", aliases=("heavy",))
def _registered_heavy_corpus(specification: CorpusSpec) -> SyntheticCorpus:
    """Build the registered heavy-locus synthetic corpus."""
    return heavy_locus_corpus(
        step=specification.step,
        seed=specification.heavy_genome_seed,
        sweep_seed=specification.sweep_seed,
        spec=specification,
    )
