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
"""Generate synthetic variation alleles."""

from __future__ import annotations

import random
from collections.abc import Sequence

from ensembl.vep.common import Strand
from ensembl.vep.corpus.synthetic.components.sequence import random_sequence
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.sequence import Slice
from ensembl.vep.models.variation import VariationAllele

__all__ = [
    "boundary_variant_alleles",
    "fuzz_variant_alleles",
]

_BASES = "ACGT"


def _snv(
    synthetic_slice: Slice, position: int, name: str, random_generator: random.Random
) -> tuple[str, VariationAllele]:
    """Build a named single-nucleotide variant at a genomic position.

    Args:
        synthetic_slice (Slice): The chromosome slice to read the reference base from.
        position (int): Genomic (slice-relative) position.
        name (str): A descriptive name for test IDs.
        random_generator (random.Random): Seeded random generator, for reproducibility.

    Returns:
        tuple[str, VariationAllele]: The variant's name and allele.
    """
    assert synthetic_slice.seq is not None
    reference_base = synthetic_slice.seq[position - 1]
    alternate_base = random_generator.choice([base for base in _BASES if base != reference_base])
    return name, VariationAllele(
        start=position,
        end=position,
        strand=Strand.FORWARD,
        slice=synthetic_slice,
        reference_allele=reference_base,
        alternate_allele=alternate_base,
    )


def _deletion(synthetic_slice: Slice, position: int, length: int, name: str) -> tuple[str, VariationAllele]:
    """Build a named deletion starting at a genomic position.

    Args:
        synthetic_slice (Slice): The chromosome slice to read the deleted bases from.
        position (int): Genomic (slice-relative) start position of the deletion.
        length (int): Number of bases deleted.
        name (str): A descriptive name for test IDs.

    Returns:
        tuple[str, VariationAllele]: The variant's name and allele.
    """
    assert synthetic_slice.seq is not None
    reference_bases = synthetic_slice.seq[position - 1 : position - 1 + length]
    return name, VariationAllele(
        start=position,
        end=position + length - 1,
        strand=Strand.FORWARD,
        slice=synthetic_slice,
        reference_allele=reference_bases,
        alternate_allele="",
    )


def _insertion(
    synthetic_slice: Slice, position: int, inserted_sequence: str, name: str
) -> tuple[str, VariationAllele]:
    """Build a named insertion immediately before a genomic position (``end == start - 1``).

    Args:
        synthetic_slice (Slice): The chromosome slice the insertion is placed on.
        position (int): Genomic (slice-relative) position the insertion is placed before.
        inserted_sequence (str): The inserted bases.
        name (str): A descriptive name for test IDs.

    Returns:
        tuple[str, VariationAllele]: The variant's name and allele.
    """
    return name, VariationAllele(
        start=position,
        end=position - 1,
        strand=Strand.FORWARD,
        slice=synthetic_slice,
        reference_allele="",
        alternate_allele=inserted_sequence,
    )


def _splice_boundary_variants(
    synthetic_slice: Slice,
    transcript: Transcript,
    random_generator: random.Random,
) -> list[tuple[str, VariationAllele]]:
    """Build the splice-site and intronic boundary variants for one transcript."""
    named_variants: list[tuple[str, VariationAllele]] = []
    tag = transcript.stable_id
    strand_step = 1 if transcript.strand == Strand.FORWARD else -1

    for intron_index, intron in enumerate(transcript.introns):
        five_prime_pos, three_prime_pos = (
            (intron.start, intron.end) if transcript.strand == Strand.FORWARD else (intron.end, intron.start)
        )
        for offset in (1, 2, 3, 4, 5, 6, 8):
            named_variants.append(
                _snv(
                    synthetic_slice,
                    five_prime_pos + strand_step * (offset - 1),
                    f"{tag}_intron{intron_index}_donor+{offset}",
                    random_generator,
                )
            )
        for offset in (1, 2, 3, 5, 8, 17):
            named_variants.append(
                _snv(
                    synthetic_slice,
                    three_prime_pos - strand_step * (offset - 1),
                    f"{tag}_intron{intron_index}_acceptor-{offset}",
                    random_generator,
                )
            )
        for offset in (1, 2, 3, 4):
            named_variants.append(
                _snv(
                    synthetic_slice,
                    five_prime_pos - strand_step * offset,
                    f"{tag}_intron{intron_index}_exon_before_donor{offset}",
                    random_generator,
                )
            )
            named_variants.append(
                _snv(
                    synthetic_slice,
                    three_prime_pos + strand_step * offset,
                    f"{tag}_intron{intron_index}_exon_after_acceptor{offset}",
                    random_generator,
                )
            )
        deep_intron_midpoint = (intron.start + intron.end) // 2
        named_variants.append(
            _snv(
                synthetic_slice,
                deep_intron_midpoint,
                f"{tag}_intron{intron_index}_deep",
                random_generator,
            )
        )
    return named_variants


def _coding_boundary_variants(
    synthetic_slice: Slice,
    transcript: Transcript,
    random_generator: random.Random,
) -> list[tuple[str, VariationAllele]]:
    """Build codon, CDS, and UTR boundary variants for one coding transcript."""
    named_variants: list[tuple[str, VariationAllele]] = []
    tag = transcript.stable_id
    for codon_name, codon_spans in (
        ("start", transcript.start_codon_genomic_spans()),
        ("stop", transcript.stop_codon_genomic_spans()),
    ):
        if codon_spans is None:
            continue
        for span_start, span_end in codon_spans:
            for position in range(span_start, span_end + 1):
                named_variants.append(
                    _snv(
                        synthetic_slice,
                        position,
                        f"{tag}_{codon_name}_codon_pos{position}",
                        random_generator,
                    )
                )
        # A 3bp genomic deletion/insertion at codon_spans[0][0] only actually deletes/edits
        # the codon when it is one contiguous genomic run. A codon split across introns cannot
        # be represented by a single genomic edit, so these cases are single-span only.
        if len(codon_spans) == 1:
            codon_start = codon_spans[0][0]
            named_variants.append(
                _deletion(synthetic_slice, codon_start, 3, f"{tag}_{codon_name}_codon_del3")
            )
            named_variants.append(
                _deletion(synthetic_slice, codon_start, 1, f"{tag}_{codon_name}_codon_del1")
            )
            named_variants.append(
                _insertion(synthetic_slice, codon_start, "GGG", f"{tag}_{codon_name}_codon_ins3")
            )

    cds_length = transcript.cds_length or 0
    for cds_offset in range(3, cds_length - 3, 17):
        genomic_position = transcript.cds_offset_to_genomic(cds_offset)
        if genomic_position is None:
            continue
        named_variants.append(
            _snv(synthetic_slice, genomic_position, f"{tag}_cds{cds_offset}", random_generator)
        )
        named_variants.append(_deletion(synthetic_slice, genomic_position, 1, f"{tag}_cds{cds_offset}_del1"))
        named_variants.append(_deletion(synthetic_slice, genomic_position, 3, f"{tag}_cds{cds_offset}_del3"))
        named_variants.append(
            _insertion(synthetic_slice, genomic_position, "A", f"{tag}_cds{cds_offset}_ins1")
        )
        named_variants.append(
            _insertion(synthetic_slice, genomic_position, "ACG", f"{tag}_cds{cds_offset}_ins3")
        )

    utr5_position, utr3_position = (
        (transcript.start + 5, transcript.end - 5)
        if transcript.strand == Strand.FORWARD
        else (transcript.end - 5, transcript.start + 5)
    )
    named_variants.append(_snv(synthetic_slice, utr5_position, f"{tag}_utr5_edge", random_generator))
    named_variants.append(_snv(synthetic_slice, utr3_position, f"{tag}_utr3_edge", random_generator))
    return named_variants


def _non_coding_boundary_variants(
    synthetic_slice: Slice,
    transcript: Transcript,
    random_generator: random.Random,
) -> list[tuple[str, VariationAllele]]:
    """Build the representative exonic variant for one non-coding transcript."""
    return [
        _snv(
            synthetic_slice,
            transcript.spliced_exons[0].start + 10,
            f"{transcript.stable_id}_noncoding_exon",
            random_generator,
        )
    ]


def _transcript_flank_and_span_variants(
    synthetic_slice: Slice,
    transcript: Transcript,
    chromosome_length: int,
    random_generator: random.Random,
) -> list[tuple[str, VariationAllele]]:
    """Build flank and whole-transcript structural variants for one transcript."""
    named_variants: list[tuple[str, VariationAllele]] = []
    tag = transcript.stable_id
    for flank_distance in (1, 100, 4999, 5000, 5001, 8000):
        named_variants.append(
            _snv(
                synthetic_slice,
                max(1, transcript.start - flank_distance),
                f"{tag}_upstream{flank_distance}",
                random_generator,
            )
        )
        named_variants.append(
            _snv(
                synthetic_slice,
                min(chromosome_length, transcript.end + flank_distance),
                f"{tag}_downstream{flank_distance}",
                random_generator,
            )
        )
    named_variants.append(
        _deletion(
            synthetic_slice,
            transcript.start - 10,
            (transcript.end - transcript.start) + 20,
            f"{tag}_ablation",
        )
    )
    named_variants.append(_insertion(synthetic_slice, transcript.start - 10, "A" * 50, f"{tag}_elongation"))
    # Zero-length "insertion at a boundary" case (end == start - 1) exercised directly.
    named_variants.append(_insertion(synthetic_slice, transcript.start, "T", f"{tag}_insertion_at_start"))
    return named_variants


def boundary_variant_alleles(
    synthetic_slice: Slice, transcripts: Sequence[Transcript], seed: int = 11
) -> list[tuple[str, VariationAllele]]:
    """Build one variant per interesting structural offset, for every transcript.

    Mirrors ``socore.data.boundary_variants``: every splice donor/acceptor offset, the exonic
    side of every splice boundary, a deep-intronic control, every start/stop codon base (plus
    whole-codon deletions/insertions), a sweep through the CDS, UTR edges, upstream/downstream
    flanks (including the 4999/5000/5001 boundary), and whole-transcript ablation/elongation.

    Args:
        synthetic_slice (Slice): The chromosome slice built by the default corpus generator.
        transcripts (Sequence[Transcript]): The transcripts built by the default corpus generator.
        seed (int): Random seed, for reproducibility.

    Returns:
        list[tuple[str, VariationAllele]]: Named boundary variants.
    """
    random_generator = random.Random(seed)
    named_variants: list[tuple[str, VariationAllele]] = []
    chromosome_length = synthetic_slice.end

    for transcript in transcripts:
        named_variants.extend(_splice_boundary_variants(synthetic_slice, transcript, random_generator))
        if transcript.translation is not None:
            named_variants.extend(_coding_boundary_variants(synthetic_slice, transcript, random_generator))
        else:
            named_variants.extend(
                _non_coding_boundary_variants(synthetic_slice, transcript, random_generator)
            )
        named_variants.extend(
            _transcript_flank_and_span_variants(
                synthetic_slice, transcript, chromosome_length, random_generator
            )
        )

    return named_variants


def fuzz_variant_alleles(
    synthetic_slice: Slice, variant_count: int = 4000, seed: int = 23
) -> list[tuple[str, VariationAllele]]:
    """Build random variants across the whole synthetic chromosome.

    Mirrors ``socore.data.fuzz_variants``.

    Args:
        synthetic_slice (Slice): The chromosome slice built by the default corpus generator.
        variant_count (int): Number of fuzz variants to generate.
        seed (int): Random seed, for reproducibility.

    Returns:
        list[tuple[str, VariationAllele]]: Named fuzz variants.

    Raises:
        ValueError: If ``variant_count`` is negative.
    """
    if variant_count < 0:
        raise ValueError(f"variant_count must be non-negative, got {variant_count}")
    random_generator = random.Random(seed)
    named_variants: list[tuple[str, VariationAllele]] = []
    chromosome_length = synthetic_slice.end
    for index in range(variant_count):
        position = random_generator.randint(20, chromosome_length - 60)
        variant_kind_roll = random_generator.random()
        if variant_kind_roll < 0.6:
            named_variants.append(_snv(synthetic_slice, position, f"fuzz{index}", random_generator))
        elif variant_kind_roll < 0.8:
            deletion_length = random_generator.choice([1, 2, 3, 4, 6, 12])
            named_variants.append(_deletion(synthetic_slice, position, deletion_length, f"fuzz{index}"))
        else:
            insertion_length = random_generator.choice([1, 2, 3, 5, 9])
            named_variants.append(
                _insertion(
                    synthetic_slice,
                    position,
                    random_sequence(insertion_length, random_generator),
                    f"fuzz{index}",
                )
            )
    return named_variants
