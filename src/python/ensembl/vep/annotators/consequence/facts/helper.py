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
"""Unregistered helper functions used by consequence facts."""

from __future__ import annotations

from typing import cast

from ensembl.vep.annotators.consequence.facts.base import get_fact
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.common import DEFAULT_DISTANCE, ClassSOTerm, Strand
from ensembl.vep.common.constants import EMPTY_ALLELE
from ensembl.vep.models.genic import Intron
from ensembl.vep.models.variation import VariationAllele
from ensembl.vep.models.overlap import FeatureOverlap

__all__ = [
    "DOWNSTREAM_DISTANCE",
    "UPSTREAM_DISTANCE",
    "allele_length_difference",
    "alternate_allele",
    "alternate_allele_length",
    "downstream_distance",
    "genomic_overlap",
    "intron_offsets_to_genomic_span",
    "reference_allele",
    "reference_allele_length",
    "upstream_distance",
    "variant_genomic_span",
]

UPSTREAM_DISTANCE = DOWNSTREAM_DISTANCE = DEFAULT_DISTANCE
"""The distances a run uses when it configures none, matching Perl VEP's ``--distance`` default."""


def upstream_distance(fact_context: FactContext) -> int:
    """Return the run's upstream pairing distance, or the default outside a run.

    Read from the annotator context because ``--distance`` is per run, and a fact evaluated
    standalone (a unit test, say) has no context to read.
    """

    context = fact_context.annotator_context
    return UPSTREAM_DISTANCE if context is None else context.upstream_distance


def downstream_distance(fact_context: FactContext) -> int:
    """Return the run's downstream pairing distance, or the default outside a run."""

    context = fact_context.annotator_context
    return DOWNSTREAM_DISTANCE if context is None else context.downstream_distance


# --------------------------------------------------------------------------
# Span and overlaps
# --------------------------------------------------------------------------


def variant_genomic_span(fact_context: FactContext[FeatureOverlap]) -> tuple[int, int]:
    """Return the variant's genomic ``(minimum, maximum)`` span.

    Thin wrapper over ``FactContext.variant_span``, which memoises the result for the whole
    overlap -- kept as a free function because every position/sequence fact calls it by this name.
    """
    return fact_context.variant_span


def genomic_overlap(range_start: int, range_end: int, other_start: int, other_end: int) -> bool:
    """Return whether two closed genomic intervals share a position."""
    return range_start <= other_end and range_end >= other_start


# --------------------------------------------------------------------------
# Intron offsets
# --------------------------------------------------------------------------


def intron_offsets_to_genomic_span(
    intron: Intron, strand: Strand, low: int, high: int, from_donor: bool
) -> tuple[int, int]:
    """Return the genomic span of transcript-relative intron offsets.

    As usual the genomic span is ascending.

    Offsets are 1-based, transcript sense, measured from the intron's 5' end (donor) when
    ``from_donor`` is True, else from its 3' end (acceptor). The donor is ``intron.start`` on the
    forward strand and ``intron.end`` on the reverse strand (the intron's 5' end in transcript
    sense is the high genomic coordinate there).
    """
    if strand == Strand.FORWARD:
        five_prime_pos, three_prime_pos, step = intron.start, intron.end, 1
    else:
        five_prime_pos, three_prime_pos, step = intron.end, intron.start, -1
    anchor = five_prime_pos if from_donor else three_prime_pos
    direction = step if from_donor else -step
    position_a = anchor + direction * (low - 1)
    position_b = anchor + direction * (high - 1)
    return (position_a, position_b) if position_a <= position_b else (position_b, position_a)


def intron_offset_hit(
    fact_context: FactContext[FeatureOverlap], low: int, high: int, from_donor: bool
) -> bool:
    """Return whether the variant hits the given offsets in any nearby intron.

    Reads ``_nearby_introns`` via ``get_fact`` rather than calling it directly: the ``@fact``
    decorator returns the original, unguarded, unmemoised function unchanged, so a bare call
    would silently skip both its `requires` guard and the memo table.
    """
    variant_min, variant_max = variant_genomic_span(fact_context)
    strand = fact_context.overlap.feature.strand
    nearby_introns = cast("tuple[Intron, ...]", get_fact(fact_context, "_nearby_introns"))
    for intron in nearby_introns:
        span_start, span_end = intron_offsets_to_genomic_span(intron, strand, low, high, from_donor)
        if genomic_overlap(variant_min, variant_max, span_start, span_end):
            return True
    return False


# --------------------------------------------------------------------------
# Alleles: reference and alternate allele (only relevant for short variants not structural)
# --------------------------------------------------------------------------


def _variant(fact_context: FactContext[FeatureOverlap]) -> VariationAllele:
    """The overlapping variation allele, narrowed from ``Feature``.

    A cast solely because ``FeatureOverlap.overlapping_feature`` is typed with wider ``Feature`` rather than
    ``VariationAllele``.

    Args:
        fact_context (FactContext[FeatureOverlap]): The overlap being evaluated.

    Returns:
        VariationAllele: The overlapping feature, as the variation allele it always is.
    """
    return cast("VariationAllele", fact_context.overlap.overlapping_feature)


def reference_allele(fact_context: FactContext[FeatureOverlap]) -> str:
    """Return the reference allele, treating ``None`` as empty."""
    return _variant(fact_context).reference_allele or ""


def alternate_allele(fact_context: FactContext[FeatureOverlap]) -> str:
    """Return the alternate allele, treating ``None`` as empty."""
    return _variant(fact_context).alternate_allele or ""


def reference_allele_length(fact_context: FactContext[FeatureOverlap]) -> int:
    """Return the reference-allele length, treating empty spellings as zero."""
    allele = _variant(fact_context).reference_allele or ""
    return len(allele.replace(EMPTY_ALLELE, ""))


def alternate_allele_length(fact_context: FactContext[FeatureOverlap]) -> int:
    """Return the alternate-allele length, treating empty spellings as zero."""
    allele = _variant(fact_context).alternate_allele or ""
    return len(allele.replace(EMPTY_ALLELE, ""))


def allele_length_difference(fact_context: FactContext[FeatureOverlap]) -> int:
    """Return the length change from reference to alternate allele."""
    return alternate_allele_length(fact_context) - reference_allele_length(fact_context)


# --------------------------------------------------------------------------
# Alleles: codons, peptides and full CDS
# --------------------------------------------------------------------------


def get_codon_alleles(fact_context: FactContext[FeatureOverlap]) -> tuple[str, str] | None:
    """Return the internal codon-allele fact with its concrete type."""
    return cast("tuple[str, str] | None", get_fact(fact_context, "_codon_alleles"))


def get_peptide_alleles(fact_context: FactContext[FeatureOverlap]) -> tuple[str, str] | None:
    """Return the internal peptide-allele fact with its concrete type."""
    return cast("tuple[str, str] | None", get_fact(fact_context, "_peptide_alleles"))


def get_full_cds_peptide_alleles(fact_context: FactContext[FeatureOverlap]) -> tuple[str, str] | None:
    """Return the internal full-CDS peptide fact with its concrete type."""
    return cast("tuple[str, str] | None", get_fact(fact_context, "_full_cds_peptide_alleles"))


# --------------------------------------------------------------------------
# SO term
# --------------------------------------------------------------------------


def variant_class_so_term(fact_context: FactContext[FeatureOverlap]) -> ClassSOTerm | None:
    """Return the class SO term from the variation that owns the overlapping allele."""
    variation = _variant(fact_context).variation
    return variation.class_so_term if variation is not None else None
