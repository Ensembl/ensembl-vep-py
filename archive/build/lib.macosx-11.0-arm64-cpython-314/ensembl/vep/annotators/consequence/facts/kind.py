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
"""Kind facts: describing what kind of features the context has."""

from __future__ import annotations

from ensembl.vep.annotators.consequence.facts.base import fact, get_fact
from ensembl.vep.annotators.consequence.facts.helper import (
    reference_allele,
    alternate_allele,
    alternate_allele_length,
    reference_allele_length,
    variant_class_so_term,
)
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.common import ClassSOTerm, CostTier
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


@fact(CostTier.FREE)
def is_transcript(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The overlapped feature is a Transcript."""
    return isinstance(fact_context.overlap.feature, Transcript)


@fact(CostTier.FREE)
def is_regulatory(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The overlapped feature is a RegulatoryFeature"""
    return isinstance(fact_context.overlap.feature, Regulatory)


@fact(CostTier.FREE)
def is_motif(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The overlapped feature is a MotifFeature."""
    return isinstance(fact_context.overlap.feature, Motif)


@fact(CostTier.FREE)
def is_structural(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The overlapping variant allele is a structural variant."""
    return isinstance(fact_context.overlap.overlapping_feature, StructuralVariationAllele)


@fact(CostTier.FREE)
def is_short(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The overlapping variant allele is a short variant, i.e. not a structural one."""
    return not isinstance(fact_context.overlap.overlapping_feature, StructuralVariationAllele)


# --------------------------------------------------------------------------
# Biotype: what kind of transcript this is.
# --------------------------------------------------------------------------


@fact(CostTier.FREE, requires=("is_transcript",))
def has_translation(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The transcript has a translation (produces protein)."""
    return fact_context.overlap.feature.translation is not None


@fact(CostTier.FREE, requires=("is_transcript",))
def biotype_protein_coding(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The transcript has a biotype of protein coding."""
    return fact_context.biotype == "protein_coding"


@fact(CostTier.FREE, requires=("is_transcript",))
def biotype_nmd(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The transcript has a biotype of nonsense mediated decay."""
    return fact_context.biotype == "nonsense_mediated_decay"


@fact(CostTier.FREE, requires=("is_transcript",))
def biotype_non_coding(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The transcript has no translation and does not have biotype protein coding."""
    return fact_context.biotype != "protein_coding" and fact_context.overlap.feature.translation is None


@fact(CostTier.FREE, requires=("is_transcript",))
def biotype_miRNA(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The transcript has biotype of miRNA."""
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


@fact(CostTier.FREE, requires=("is_short",))
def is_short_insertion(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The variant is a short insertion."""
    return reference_allele(fact_context) == "-" or alternate_allele_length(
        fact_context
    ) > reference_allele_length(fact_context)


@fact(CostTier.FREE, requires=("is_short",))
def is_short_deletion(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The variant is a short deletion."""
    return alternate_allele(fact_context) == "-" or reference_allele_length(
        fact_context
    ) > alternate_allele_length(fact_context)


@fact(CostTier.FREE, deps=("is_short_insertion", "is_short_deletion"), requires=("is_short",))
def is_short_indel(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The short variant's net length change is non-zero: an insertion or a deletion or both at once."""
    return bool(get_fact(fact_context, "is_short_insertion")) or bool(
        get_fact(fact_context, "is_short_deletion")
    )


@fact(CostTier.FREE, requires=("is_structural",))
def is_structural_insertion(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The structural variant adds sequence: an insertion, a duplication, or a copy number gain."""
    return variant_class_so_term(fact_context) in (
        _INSERTION_TERMS | _DUPLICATION_TERMS | {ClassSOTerm.COPY_NUMBER_GAIN}
    )


@fact(CostTier.FREE, requires=("is_structural",))
def is_duplication(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The variant's class SO term describes a duplicated segment."""
    return variant_class_so_term(fact_context) in (_DUPLICATION_TERMS | {ClassSOTerm.COPY_NUMBER_GAIN})


@fact(CostTier.FREE, requires=("is_structural",))
def is_copy_number_gain(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The variant's class SO term describes added copies."""
    return variant_class_so_term(fact_context) == ClassSOTerm.COPY_NUMBER_GAIN


@fact(CostTier.FREE, requires=("is_structural",))
def is_structural_deletion(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The structural variant removes sequence: a deletion or a copy number loss."""
    return variant_class_so_term(fact_context) in (_DELETION_TERMS | {ClassSOTerm.COPY_NUMBER_LOSS})


@fact(CostTier.FREE, requires=("is_structural",))
def is_copy_number_loss(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The variant's class SO term describes removed copies."""
    return variant_class_so_term(fact_context) == ClassSOTerm.COPY_NUMBER_LOSS


@fact(CostTier.FREE, requires=("is_structural",))
def is_chromosomal_breakpoint(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The variant's class SO term describes a chromosomal breakpoint (undirected or directed)."""
    return variant_class_so_term(fact_context) in _CHROMOSOMAL_BREAKPOINT_TERMS


@fact(
    CostTier.FREE,
    deps=("is_structural_insertion", "is_structural_deletion"),
    requires=("is_structural",),
)
def is_structural_indel(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The structural variant's class SO term is a length-changing shape: insertion or deletion."""
    return bool(get_fact(fact_context, "is_structural_insertion")) or bool(
        get_fact(fact_context, "is_structural_deletion")
    )


@fact(CostTier.FREE, requires=("is_short",))
def is_snv(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Reference and alternate are the same non-zero length."""
    reference_length = reference_allele_length(fact_context)
    alternate_length = alternate_allele_length(fact_context)
    return reference_length == alternate_length == 1


@fact(CostTier.FREE, deps=("is_short_insertion", "is_structural_insertion"))
def is_insertion(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The variant adds sequence, whether it is a short variant or a structural one."""
    return bool(get_fact(fact_context, "is_short_insertion")) or bool(
        get_fact(fact_context, "is_structural_insertion")
    )


@fact(CostTier.FREE, deps=("is_short_deletion", "is_structural_deletion"))
def is_deletion(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The variant removes sequence, whether it is a short variant or a structural one."""
    return bool(get_fact(fact_context, "is_short_deletion")) or bool(
        get_fact(fact_context, "is_structural_deletion")
    )


@fact(CostTier.FREE, deps=("is_short_indel", "is_structural_indel"))
def is_indel(fact_context: FactContext[FeatureOverlap]) -> bool:
    """The variant's shape changes length, whether it is a short variant or a structural one."""
    return bool(get_fact(fact_context, "is_short_indel")) or bool(
        get_fact(fact_context, "is_structural_indel")
    )
