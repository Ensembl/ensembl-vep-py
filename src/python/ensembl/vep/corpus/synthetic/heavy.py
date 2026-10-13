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

from ensembl.vep.common import CoordSystemLevel, Strand, invalidate
from ensembl.vep.corpus.overlap import build_transcript_overlaps
from ensembl.vep.corpus.registry import register_corpus
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.corpus.synthetic.base import SyntheticCorpus
from ensembl.vep.corpus.synthetic.components.sequence import fill_locus_sequence
from ensembl.vep.corpus.synthetic.components.transcript import (
    build_cds_genomic_bounds,
    build_transcript,
    patch_terminal_codons,
)
from ensembl.vep.corpus.synthetic.components.variation import (
    _deletion,
    _insertion,
    _snv,
)
from ensembl.vep.models.genic import Gene, Transcript
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.models.variation import VariationAllele

__all__ = [
    "build_genome",
    "build_transcripts",
    "build_variation_alleles",
    "heavy_corpus",
]

# --------------------------------------------------------------------------
# A heavy transcript locus.
#
# The four transcripts in default corpus are disjoint, which makes them useless for measuring
# anything that shares work *between* transcripts. A real dense locus has many transcripts over
# one gene sharing most of their internal exons and differing at the edges, so this builds that
# shape (2 strands x 18 exon-sharing transcripts per strand):
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

_CHROMOSOME_NAME = "chrHeavyLocus"
_CHROMOSOME_LENGTH = 130000

_FIRST_EXONS: tuple[tuple[int, int], ...] = ((0, 200), (1000, 1200), (2000, 2200))
_CONSTITUTIVE_EXONS_UPSTREAM_OF_CASSETTE: tuple[tuple[int, int], ...] = (
    (3000, 3150),
    (4000, 4120),
)
_CASSETTE_EXON: tuple[int, int] = (4500, 4560)
_CONSTITUTIVE_EXONS_DOWNSTREAM_OF_CASSETTE: tuple[tuple[int, int], ...] = (
    (5000, 5180),
    (6000, 6090),
    (7000, 7210),
)
_LAST_EXONS: tuple[tuple[int, int], ...] = ((8000, 8300), (9000, 9400), (10000, 10300))

_LOCUS_SPAN = 10300
_FORWARD_ORIGIN = 10000
_REVERSE_ORIGIN = 70000
_CDS_START_OFFSET = 50
_CDS_STOP_OFFSET = 50
_CDS_STOP_SPACING = 6


def _biotype(transcript_index: int) -> tuple[str, bool]:
    """Return ``(biotype, is_coding)`` for a heavy-locus transcript, by its index in the locus.

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


def _build_transcript(
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
    biotype, is_coding = _biotype(transcript_index)

    translation_bounds: tuple[int, int] | None = None
    if is_coding:
        frame_index = first_exon_index if strand == Strand.FORWARD else last_exon_index
        translation_bounds = build_cds_genomic_bounds(
            exons_genomic_ascending,
            strand=strand,
            frame_index=frame_index,
            start_offset=_CDS_START_OFFSET,
            stop_offset=_CDS_STOP_OFFSET + _CDS_STOP_SPACING * transcript_index,
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
        cds_genomic_bounds=translation_bounds,
        gene=gene,
    )


def build_transcripts(
    synthetic_slice: Slice, *, origin: int, strand: Strand, stable_id_tag: str
) -> list[Transcript]:
    """Build 18 transcripts sharing constitutive exons over a single locus.

    Args:
        synthetic_slice (Slice): The shared synthetic-chromosome slice.
        origin (int): Genomic (slice-relative) offset the locus's exon coordinates are relative
            to.
        strand (Strand): The strand every transcript in this locus is built on.
        stable_id_tag (str): Prefix for each transcript's stable_id (e.g. "FWD", "REV").

    Returns:
        list[Transcript]: The 18 transcripts, not yet carrying real start/stop codon sequence
            (see build_genome, which patches that in afterwards).
    """
    transcripts: list[Transcript] = []
    transcript_index = 0
    for first_exon_index, first_exon in enumerate(_FIRST_EXONS):
        for include_cassette_exon in (False, True):
            for last_exon_index, last_exon in enumerate(_LAST_EXONS):
                exon_blocks = (
                    [first_exon]
                    + list(_CONSTITUTIVE_EXONS_UPSTREAM_OF_CASSETTE)
                    + ([_CASSETTE_EXON] if include_cassette_exon else [])
                    + list(_CONSTITUTIVE_EXONS_DOWNSTREAM_OF_CASSETTE)
                    + [last_exon]
                )
                exons_genomic_ascending = sorted(
                    (origin + block_start, origin + block_end) for block_start, block_end in exon_blocks
                )
                transcripts.append(
                    _build_transcript(
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


def build_genome(seed: int = 31) -> tuple[Slice, tuple[Transcript, ...]]:
    """Build the two-locus (one per strand), 36-transcript heavy-locus synthetic chromosome.

    Mirrors ``socore.data.build_locus_genome``.

    Args:
        seed (int): Random seed, for reproducibility.

    Returns:
        tuple[Slice, tuple[Transcript, ...]]: The shared chromosome slice (with real, codon-patched
            sequence) and all 36 transcripts (18 forward-strand, 18 reverse-strand).
    """
    random_generator = random.Random(seed)
    chromosome_chars = [""] * _CHROMOSOME_LENGTH
    fill_locus_sequence(
        chromosome_chars,
        1,
        _CHROMOSOME_LENGTH,
        random_generator,
        bases="ACGT",
    )

    synthetic_slice = Slice(
        coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
        seq_region_name=_CHROMOSOME_NAME,
        start=1,
        end=_CHROMOSOME_LENGTH,
        strand=Strand.FORWARD,
    )
    with synthetic_slice.unsafe() as unsafe_slice:
        setattr(unsafe_slice, "_seq", "".join(chromosome_chars))

    forward_transcripts = build_transcripts(
        synthetic_slice, origin=_FORWARD_ORIGIN, strand=Strand.FORWARD, stable_id_tag="FWD"
    )
    reverse_transcripts = build_transcripts(
        synthetic_slice, origin=_REVERSE_ORIGIN, strand=Strand.REVERSE, stable_id_tag="REV"
    )
    transcripts = tuple(forward_transcripts + reverse_transcripts)

    fill_locus_sequence(
        chromosome_chars,
        _FORWARD_ORIGIN,
        _FORWARD_ORIGIN + _LOCUS_SPAN,
        random_generator,
    )
    fill_locus_sequence(
        chromosome_chars,
        _REVERSE_ORIGIN,
        _REVERSE_ORIGIN + _LOCUS_SPAN,
        random_generator,
    )
    patch_terminal_codons(chromosome_chars, transcripts)

    with synthetic_slice.unsafe() as unsafe_slice:
        setattr(unsafe_slice, "_seq", "".join(chromosome_chars))

    # Unlike the default corpus's four disjoint transcripts, transcripts here legitimately
    # share exons in different reading frames, so one transcript's real start/stop codon
    # can land as a premature or missing stop in a sibling's frame. That is expected, not a data
    # bug. Only the coding-sequence cache needs invalidating here.
    for transcript in transcripts:
        if transcript.translation is not None:
            invalidate(transcript, "_coding_sequence")

    return synthetic_slice, transcripts


def build_variation_alleles(
    synthetic_slice: Slice, step: int = 13, seed: int = 41
) -> tuple[tuple[str, VariationAllele], ...]:
    """Build a dense sweep of variants across each locus, so shared exons are hit repeatedly.

    Args:
        synthetic_slice (Slice): The chromosome slice built by build_genome.
        step (int): Genomic spacing between consecutive sweep positions.
        seed (int): Random seed, for reproducibility.

    Returns:
        tuple[tuple[str, VariationAllele], ...]: Named heavy-locus sweep variants.

    Raises:
        ValueError: If ``step`` is not positive.
    """
    if step <= 0:
        raise ValueError(f"step must be positive, got {step}")
    random_generator = random.Random(seed)
    named_variants: list[tuple[str, VariationAllele]] = []
    for origin, tag in ((_FORWARD_ORIGIN, "FWD"), (_REVERSE_ORIGIN, "REV")):
        for sweep_index, position in enumerate(range(origin, origin + _LOCUS_SPAN, step)):
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
    return tuple(named_variants)


def heavy_corpus(
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
        spec (CorpusSpec | None): Construction specification retained by the corpus.

    Returns:
        SyntheticCorpus: The synthetic slice, its 36 transcripts, variants, and overlaps.
    """
    synthetic_slice, transcripts = build_genome(seed=seed)
    named_variants = build_variation_alleles(synthetic_slice, step=step, seed=sweep_seed)
    named_overlaps = build_transcript_overlaps(named_variants, transcripts)
    return SyntheticCorpus(
        synthetic_slice=synthetic_slice,
        transcripts=transcripts,
        named_variants=named_variants,
        named_overlaps=named_overlaps,
        spec=spec or CorpusSpec(kind="synthetic.heavy"),
    )


@register_corpus("synthetic.heavy", aliases=("heavy",))
def _registered_heavy_corpus(specification: CorpusSpec) -> SyntheticCorpus:
    """Build the registered heavy-locus synthetic corpus."""
    return heavy_corpus(
        step=specification.step,
        seed=specification.heavy_genome_seed,
        sweep_seed=specification.sweep_seed,
        spec=specification,
    )
