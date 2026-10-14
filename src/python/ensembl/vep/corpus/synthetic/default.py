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
"""Build the default synthetic corpus."""

from __future__ import annotations

import random
from collections.abc import Sequence

from ensembl.vep.common import CoordSystemLevel, Strand
from ensembl.vep.common.codons import translate_dna
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
    boundary_variant_alleles,
    fuzz_variant_alleles,
)
from ensembl.vep.models.genic import Gene, Transcript
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.models.variation import VariationAllele

__all__ = [
    "build_genome",
    "build_variation_alleles",
    "default_corpus",
]

_CHROMOSOME_NAME = "chrSynthetic"
_CHROMOSOME_LENGTH = 60000

_ENST0001_EXONS: tuple[tuple[int, int], ...] = ((1000, 1200), (2000, 2200), (3000, 3300))
_ENST0002_EXONS: tuple[tuple[int, int], ...] = ((10000, 10300), (11000, 11200), (12000, 12150))
_ENST0003_EXONS: tuple[tuple[int, int], ...] = ((20000, 20200), (21000, 21300))
_ENST0004_EXONS: tuple[tuple[int, int], ...] = ((30000, 30200), (31000, 31250))

_CDS_START_OFFSET = 50
_ENST0001_CDS_FRAME_INDEX = 1
_ENST0001_CDS_STOP_OFFSET = 101
_ENST0002_CDS_FRAME_INDEX = 2
_ENST0002_CDS_STOP_OFFSET = 102
_ENST0004_CDS_FRAME_INDEX = 1
_ENST0004_CDS_STOP_OFFSET = 51


def build_genome(seed: int = 7) -> tuple[Slice, tuple[Transcript, ...]]:
    """Build the synthetic chromosome and its four transcripts.

    Mirrors ``socore.data.build_genome``: one forward-strand protein-coding transcript, one
    reverse-strand protein-coding transcript, one non-coding (lincRNA-like) transcript, and one
    nonsense_mediated_decay transcript.

    The biotypes are declared on purposely different features so that both branches of
    ``FactContext.biotype`` (see consequence/context/fact.py) are covered: ENST0001 and ENST0003
    carry their own, ENST0002 and ENST0004 leave theirs None and inherit the gene's.

    Args:
        seed (int): Random seed, for reproducibility.

    Returns:
        tuple[Slice, tuple[Transcript, ...]]: The shared chromosome slice (with real, CDS-patched
            sequence) and the four transcripts.
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

    # ENST0002 deliberately declares no biotype of its own, so that it exercises
    # FactContext.biotype's fall back to the gene; ENST0001 declares its own and carries no gene.
    coding_gene = Gene(
        start=_ENST0002_EXONS[0][0],
        end=_ENST0002_EXONS[-1][1],
        strand=Strand.REVERSE,
        stable_id="ENSG0002",
        slice=synthetic_slice,
    )
    coding_gene.biotype = "protein_coding"

    nmd_gene = Gene(
        start=_ENST0004_EXONS[0][0],
        end=_ENST0004_EXONS[-1][1],
        strand=Strand.FORWARD,
        stable_id="ENSG0004",
        slice=synthetic_slice,
    )
    nmd_gene.biotype = "nonsense_mediated_decay"

    transcripts = (
        build_transcript(
            stable_id="ENST0001",
            synthetic_slice=synthetic_slice,
            strand=Strand.FORWARD,
            exons_genomic_ascending=_ENST0001_EXONS,
            cds_genomic_bounds=build_cds_genomic_bounds(
                _ENST0001_EXONS,
                strand=Strand.FORWARD,
                frame_index=_ENST0001_CDS_FRAME_INDEX,
                start_offset=_CDS_START_OFFSET,
                stop_offset=_ENST0001_CDS_STOP_OFFSET,
            ),
            gene=None,
            biotype="protein_coding",
        ),
        build_transcript(
            stable_id="ENST0002",
            synthetic_slice=synthetic_slice,
            strand=Strand.REVERSE,
            exons_genomic_ascending=_ENST0002_EXONS,
            cds_genomic_bounds=build_cds_genomic_bounds(
                _ENST0002_EXONS,
                strand=Strand.REVERSE,
                frame_index=_ENST0002_CDS_FRAME_INDEX,
                start_offset=_CDS_START_OFFSET,
                stop_offset=_ENST0002_CDS_STOP_OFFSET,
            ),
            gene=coding_gene,
        ),
        build_transcript(
            stable_id="ENST0003",
            synthetic_slice=synthetic_slice,
            strand=Strand.FORWARD,
            exons_genomic_ascending=_ENST0003_EXONS,
            cds_genomic_bounds=None,
            gene=None,
            biotype="lincRNA",
        ),
        build_transcript(
            stable_id="ENST0004",
            synthetic_slice=synthetic_slice,
            strand=Strand.FORWARD,
            exons_genomic_ascending=_ENST0004_EXONS,
            cds_genomic_bounds=build_cds_genomic_bounds(
                _ENST0004_EXONS,
                strand=Strand.FORWARD,
                frame_index=_ENST0004_CDS_FRAME_INDEX,
                start_offset=_CDS_START_OFFSET,
                stop_offset=_ENST0004_CDS_STOP_OFFSET,
            ),
            gene=nmd_gene,
        ),
    )

    for transcript in transcripts:
        fill_locus_sequence(
            chromosome_chars,
            transcript.start,
            transcript.end,
            random_generator,
        )
    patch_terminal_codons(chromosome_chars, transcripts)

    with synthetic_slice.unsafe() as unsafe_slice:
        setattr(unsafe_slice, "_seq", "".join(chromosome_chars))

    for transcript in transcripts:
        if transcript.translation is None:
            continue
        written_coding_sequence = transcript.coding_sequence
        assert written_coding_sequence is not None
        translated = translate_dna(written_coding_sequence)
        assert translated.startswith("M"), f"{transcript.stable_id} does not start with M: {translated[:5]}"
        assert translated.endswith("*"), f"{transcript.stable_id} does not end with a stop"
        assert "*" not in translated[:-1], f"{transcript.stable_id} has an internal stop"

    return synthetic_slice, transcripts


def build_variation_alleles(
    synthetic_slice: Slice,
    transcripts: Sequence[Transcript],
    *,
    fuzz_count: int = 4000,
    boundary_seed: int = 11,
    fuzz_seed: int = 23,
) -> tuple[tuple[str, VariationAllele], ...]:
    """Build named boundary and random alleles for the default corpus.

    Args:
        synthetic_slice (Slice): Sequence region holding the variants.
        transcripts (Sequence[Transcript]): Transcripts used to derive boundary variants.
        fuzz_count (int): Number of additional random variants.
        boundary_seed (int): Boundary-variant random seed.
        fuzz_seed (int): Random-variant seed.

    Returns:
        tuple[tuple[str, VariationAllele], ...]: Boundary variants followed by random variants.
    """
    return tuple(
        boundary_variant_alleles(synthetic_slice, transcripts, seed=boundary_seed)
        + fuzz_variant_alleles(synthetic_slice, fuzz_count, seed=fuzz_seed)
    )


def default_corpus(
    fuzz_count: int = 4000,
    seed: int = 7,
    *,
    boundary_seed: int = 11,
    fuzz_seed: int = 23,
    spec: CorpusSpec | None = None,
) -> SyntheticCorpus:
    """Build the default synthetic corpus: four transcripts, boundary variants, and fuzz variants.

    Mirrors ``socore.data.default_dataset``.

    Args:
        fuzz_count (int): Number of additional random ("fuzz") variants to add.
        seed (int): Random seed for the underlying chromosome, for reproducibility.
        boundary_seed (int): Random seed for the boundary variants' alternate bases.
        fuzz_seed (int): Random seed for the fuzz variants.
        spec (CorpusSpec | None): Construction specification retained by the corpus.

    Returns:
        SyntheticCorpus: The synthetic slice, its transcripts, variants, and overlaps.
    """
    synthetic_slice, transcripts = build_genome(seed=seed)
    named_variants = build_variation_alleles(
        synthetic_slice,
        transcripts,
        fuzz_count=fuzz_count,
        boundary_seed=boundary_seed,
        fuzz_seed=fuzz_seed,
    )
    named_overlaps = build_transcript_overlaps(named_variants, transcripts)
    return SyntheticCorpus(
        synthetic_slice=synthetic_slice,
        transcripts=transcripts,
        named_variants=named_variants,
        named_overlaps=named_overlaps,
        spec=spec or CorpusSpec(kind="synthetic.default"),
    )


@register_corpus("synthetic.default", aliases=("default",))
def _registered_default_corpus(specification: CorpusSpec) -> SyntheticCorpus:
    """Build the registered default synthetic corpus."""
    return default_corpus(
        fuzz_count=specification.fuzz,
        seed=specification.genome_seed,
        boundary_seed=specification.boundary_seed,
        fuzz_seed=specification.fuzz_seed,
        spec=specification,
    )
