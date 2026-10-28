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
"""Facts describing variant and feature kinds."""

from __future__ import annotations

from ensembl.vep.annotators.consequence.facts.base import get_fact, register_fact
from ensembl.vep.annotators.consequence.facts.helper import (
    reference_allele,
    alternate_allele,
    alternate_allele_length,
    reference_allele_length,
    variant_class_so_term,
)
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.common import ClassSOTerm, CostTier
from ensembl.vep.common.constants import EMPTY_ALLELE
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.regulatory import Motif, Regulatory
from ensembl.vep.models.variation import StructuralVariationAllele
from ensembl.vep.models.overlap import FeatureOverlap, TranscriptVariationAllele

__all__ = [
    "biotype_miRNA",
    "biotype_nmd",
    "biotype_non_coding",
    "biotype_protein_coding",
    "has_translation",
    "is_chromosomal_breakpoint",
    "is_copy_number_gain",
    "is_copy_number_loss",
    "is_deletion",
    "is_duplication",
    "is_indel",
    "is_insertion",
    "is_motif",
    "is_regulatory",
    "is_short",
    "is_short_deletion",
    "is_short_indel",
    "is_short_insertion",
    "is_snv",
    "is_structural",
    "is_structural_deletion",
    "is_structural_indel",
    "is_structural_insertion",
    "is_transcript",
]

# --------------------------------------------------------------------------
# Feature class: what kind of feature is overlapping.
# --------------------------------------------------------------------------


@register_fact(CostTier.FREE)
def is_transcript(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the overlapped feature is a transcript."""
    return isinstance(fact_context.overlap.feature, Transcript)


@register_fact(CostTier.FREE)
def is_regulatory(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the overlapped feature is regulatory."""
    return isinstance(fact_context.overlap.feature, Regulatory)


@register_fact(CostTier.FREE)
def is_motif(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the overlapped feature is a motif."""
    return isinstance(fact_context.overlap.feature, Motif)


@register_fact(CostTier.FREE)
def is_structural(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the overlapping allele is structural."""
    return isinstance(fact_context.overlap.overlapping_feature, StructuralVariationAllele)


@register_fact(CostTier.FREE)
def is_short(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the overlapping allele is a short variant."""
    return not isinstance(fact_context.overlap.overlapping_feature, StructuralVariationAllele)


# --------------------------------------------------------------------------
# Biotype: what kind of transcript this is.
# --------------------------------------------------------------------------


@register_fact(CostTier.FREE, requires=("is_transcript",))
def has_translation(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the transcript has a translation."""
    return fact_context.overlap.feature.translation is not None


@register_fact(CostTier.FREE, requires=("is_transcript",))
def biotype_protein_coding(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the transcript has a protein-coding biotype."""
    return fact_context.biotype == "protein_coding"


@register_fact(CostTier.FREE, requires=("is_transcript",))
def biotype_nmd(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the transcript has a nonsense-mediated-decay biotype."""
    return fact_context.biotype == "nonsense_mediated_decay"


@register_fact(CostTier.FREE, requires=("is_transcript",))
def biotype_non_coding(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the transcript is non-coding."""
    return fact_context.biotype != "protein_coding" and fact_context.overlap.feature.translation is None


@register_fact(CostTier.FREE, requires=("is_transcript",))
def biotype_miRNA(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the transcript has a microRNA biotype."""
    return fact_context.biotype == "miRNA"


# --------------------------------------------------------------------------
# Variant shape: what kind of change the allele describes.
# --------------------------------------------------------------------------


_DUPLICATION_TERMS = frozenset({ClassSOTerm.DUPLICATION, ClassSOTerm.TANDEM_DUPLICATION})
"""Class SO terms that describe a duplicated segment."""

_INSERTION_TERMS = frozenset(
    {
        ClassSOTerm.INSERTION,
        ClassSOTerm.MOBILE_ELEMENT_INSERTION,
        ClassSOTerm.ALU_INSERTION,
        ClassSOTerm.HERV_INSERTION,
        ClassSOTerm.LINE1_INSERTION,
        ClassSOTerm.SVA_INSERTION,
    }
)
"""Class SO terms that add sequence without being a duplication or a copy number gain."""

_DELETION_TERMS = frozenset(
    {
        ClassSOTerm.DELETION,
        ClassSOTerm.MOBILE_ELEMENT_DELETION,
        ClassSOTerm.ALU_DELETION,
        ClassSOTerm.HERV_DELETION,
        ClassSOTerm.LINE1_DELETION,
        ClassSOTerm.SVA_DELETION,
    }
)
"""Class SO terms that remove sequence without being a copy number loss."""

_CHROMOSOMAL_BREAKPOINT_TERMS = frozenset(
    {
        ClassSOTerm.CHROMOSOME_BREAKPOINT,
        ClassSOTerm.INTERCHROMOSOMAL_BREAKPOINT,
        ClassSOTerm.INTRACHROMOSOMAL_BREAKPOINT,
    }
)
"""Class SO terms that describe a chromosomal breakpoint, undirected and both directed forms."""


@register_fact(CostTier.FREE, requires=("is_short",))
def is_short_insertion(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the variant is a short insertion."""
    return reference_allele(fact_context) == EMPTY_ALLELE or alternate_allele_length(
        fact_context
    ) > reference_allele_length(fact_context)


@register_fact(CostTier.FREE, requires=("is_short",))
def is_short_deletion(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the variant is a short deletion."""
    return alternate_allele(fact_context) == EMPTY_ALLELE or reference_allele_length(
        fact_context
    ) > alternate_allele_length(fact_context)


@register_fact(CostTier.FREE, deps=("is_short_insertion", "is_short_deletion"), requires=("is_short",))
def is_short_indel(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the short variant changes allele length."""
    return bool(get_fact(fact_context, "is_short_insertion")) or bool(
        get_fact(fact_context, "is_short_deletion")
    )


@register_fact(CostTier.FREE, requires=("is_structural",))
def is_structural_insertion(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the structural variant adds sequence."""
    return variant_class_so_term(fact_context) in (
        _INSERTION_TERMS | _DUPLICATION_TERMS | {ClassSOTerm.COPY_NUMBER_GAIN}
    )


@register_fact(CostTier.FREE, requires=("is_structural",))
def is_duplication(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the variant represents a duplicated segment."""
    return variant_class_so_term(fact_context) in (_DUPLICATION_TERMS | {ClassSOTerm.COPY_NUMBER_GAIN})


@register_fact(CostTier.FREE, requires=("is_structural",))
def is_copy_number_gain(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the variant represents a copy-number gain."""
    return variant_class_so_term(fact_context) == ClassSOTerm.COPY_NUMBER_GAIN


@register_fact(CostTier.FREE, requires=("is_structural",))
def is_structural_deletion(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the structural variant removes sequence."""
    return variant_class_so_term(fact_context) in (_DELETION_TERMS | {ClassSOTerm.COPY_NUMBER_LOSS})


@register_fact(CostTier.FREE, requires=("is_structural",))
def is_copy_number_loss(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the variant represents a copy-number loss."""
    return variant_class_so_term(fact_context) == ClassSOTerm.COPY_NUMBER_LOSS


@register_fact(CostTier.FREE, requires=("is_structural",))
def is_chromosomal_breakpoint(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the variant represents a chromosomal breakpoint."""
    return variant_class_so_term(fact_context) in _CHROMOSOMAL_BREAKPOINT_TERMS


@register_fact(
    CostTier.FREE,
    deps=("is_structural_insertion", "is_structural_deletion"),
    requires=("is_structural",),
)
def is_structural_indel(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the structural variant changes sequence length."""
    return bool(get_fact(fact_context, "is_structural_insertion")) or bool(
        get_fact(fact_context, "is_structural_deletion")
    )


@register_fact(CostTier.FREE, requires=("is_short",))
def is_snv(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the variant is a single-nucleotide substitution."""
    reference_length = reference_allele_length(fact_context)
    alternate_length = alternate_allele_length(fact_context)
    return reference_length == alternate_length == 1


@register_fact(CostTier.FREE, deps=("is_short_insertion", "is_structural_insertion"))
def is_insertion(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the variant adds sequence."""
    return bool(get_fact(fact_context, "is_short_insertion")) or bool(
        get_fact(fact_context, "is_structural_insertion")
    )


@register_fact(CostTier.FREE, deps=("is_short_deletion", "is_structural_deletion"))
def is_deletion(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the variant removes sequence."""
    return bool(get_fact(fact_context, "is_short_deletion")) or bool(
        get_fact(fact_context, "is_structural_deletion")
    )


@register_fact(CostTier.FREE, deps=("is_short_indel", "is_structural_indel"))
def is_indel(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the variant changes sequence length."""
    return bool(get_fact(fact_context, "is_short_indel")) or bool(
        get_fact(fact_context, "is_structural_indel")
    )
