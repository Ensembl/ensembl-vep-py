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

from ensembl.vep.common import CoordSystemLevel, Strand, invalidate
from ensembl.vep.common.codons import translate_dna
from ensembl.vep.models.genic import Gene, Transcript
from ensembl.vep.models.sequence import CoordSystem, Slice

from ensembl.vep.corpus.registry import register_corpus
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.corpus.synthetic.base import SyntheticCorpus
from ensembl.vep.corpus.synthetic.components.sequence import (
    _make_clean_open_reading_frame,
    _random_sequence,
    _write_transcript_cds,
    build_transcript,
)
from ensembl.vep.corpus.synthetic.components.variation import (
    boundary_variant_alleles,
    build_transcript_variation_alleles,
    fuzz_variant_alleles,
)

__all__ = [
    "build_synthetic_genome",
    "default_corpus",
]

_CHROMOSOME_LENGTH = 60000
_CHROMOSOME_NAME = "chrSynthetic"


def build_synthetic_genome(seed: int = 7) -> tuple[Slice, list[Transcript]]:
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
        tuple[Slice, list[Transcript]]: The shared chromosome slice (with real, CDS-patched
            sequence) and the four transcripts.
    """
    random_generator = random.Random(seed)
    chromosome_chars = list(_random_sequence(_CHROMOSOME_LENGTH, random_generator))

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
        start=10000, end=12150, strand=Strand.REVERSE, stable_id="ENSG0002", slice=synthetic_slice
    )
    coding_gene.biotype = "protein_coding"

    nmd_gene = Gene(
        start=30000, end=31250, strand=Strand.FORWARD, stable_id="ENSG0004", slice=synthetic_slice
    )
    nmd_gene.biotype = "nonsense_mediated_decay"

    transcripts = [
        build_transcript(
            stable_id="ENST0001",
            synthetic_slice=synthetic_slice,
            strand=Strand.FORWARD,
            exons_genomic_ascending=[(1000, 1200), (2000, 2200), (3000, 3300)],
            cds_genomic_bounds=(1050, 3199),
            gene=None,
            biotype="protein_coding",
        ),
        build_transcript(
            stable_id="ENST0002",
            synthetic_slice=synthetic_slice,
            strand=Strand.REVERSE,
            exons_genomic_ascending=[(10000, 10300), (11000, 11200), (12000, 12150)],
            cds_genomic_bounds=(10102, 12100),
            gene=coding_gene,
        ),
        build_transcript(
            stable_id="ENST0003",
            synthetic_slice=synthetic_slice,
            strand=Strand.FORWARD,
            exons_genomic_ascending=[(20000, 20200), (21000, 21300)],
            cds_genomic_bounds=None,
            gene=None,
            biotype="lincRNA",
        ),
        build_transcript(
            stable_id="ENST0004",
            synthetic_slice=synthetic_slice,
            strand=Strand.FORWARD,
            exons_genomic_ascending=[(30000, 30200), (31000, 31250)],
            cds_genomic_bounds=(30050, 31199),
            gene=nmd_gene,
        ),
    ]

    for transcript in transcripts:
        if transcript.translation is None:
            continue
        cds_length = transcript.cds_length
        assert (
            cds_length is not None and cds_length % 3 == 0
        ), f"{transcript.stable_id} CDS length {cds_length} is not a multiple of 3"
        coding_sequence = _make_clean_open_reading_frame(cds_length // 3, random_generator)
        _write_transcript_cds(chromosome_chars, transcript, coding_sequence)

    with synthetic_slice.unsafe() as unsafe_slice:
        setattr(unsafe_slice, "_seq", "".join(chromosome_chars))

    for transcript in transcripts:
        if transcript.translation is None:
            continue
        invalidate(transcript, "_coding_sequence")
        written_coding_sequence = transcript.coding_sequence
        assert written_coding_sequence is not None
        translated = translate_dna(written_coding_sequence)
        assert translated.startswith("M"), f"{transcript.stable_id} does not start with M: {translated[:5]}"
        assert translated.endswith("*"), f"{transcript.stable_id} does not end with a stop"
        assert "*" not in translated[:-1], f"{transcript.stable_id} has an internal stop"

    return synthetic_slice, transcripts


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

    Returns:
        SyntheticCorpus: The synthetic slice, its transcripts, variants, and overlaps.
    """
    synthetic_slice, transcripts = build_synthetic_genome(seed=seed)
    named_variants = boundary_variant_alleles(
        synthetic_slice, transcripts, seed=boundary_seed
    ) + fuzz_variant_alleles(synthetic_slice, fuzz_count, seed=fuzz_seed)
    named_overlaps = build_transcript_variation_alleles(named_variants, transcripts)
    return SyntheticCorpus(
        synthetic_slice=synthetic_slice,
        transcripts=tuple(transcripts),
        named_variants=tuple(named_variants),
        named_overlaps=tuple(named_overlaps),
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
