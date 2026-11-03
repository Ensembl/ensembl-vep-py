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
"""Facts describing a variant's position relative to its feature."""

from __future__ import annotations

from typing import cast

from ensembl.vep.annotators.consequence.facts.base import get_fact, register_fact
from ensembl.vep.annotators.consequence.facts.helper import (
    downstream_distance,
    upstream_distance,
    genomic_overlap,
    intron_offset_hit,
    variant_genomic_span,
)
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.common import CostTier, Strand
from ensembl.vep.models.genic import Intron
from ensembl.vep.models.overlap import FeatureOverlap, TranscriptVariationAllele

__all__ = [
    "contain_feature",
    "is_downstream",
    "is_upstream",
    "overlap_acceptor_site",
    "overlap_cds",
    "overlap_donor_5th_base",
    "overlap_donor_region",
    "overlap_donor_site",
    "overlap_exon",
    "overlap_feature",
    "overlap_five_prime_utr",
    "overlap_intron",
    "overlap_intron_boundary",
    "overlap_intron_interior",
    "overlap_polypyrimidine_tract",
    "overlap_splice_region",
    "overlap_start_codon",
    "overlap_stop_codon",
    "overlap_three_prime_utr",
    "partial_last_codon",
    "within_any_exon",
    "within_feature",
]

# --------------------------------------------------------------------------
# Placement relative to the whole feature.
# --------------------------------------------------------------------------


@register_fact(CostTier.BOUNDS)
def overlap_feature(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the variant overlaps the feature's genomic span.

    The variant can
        - just intersect the feature's span
        - be fully contained in the feature's span
        - or, fully contain the feature's span
    """
    variant_min, variant_max = variant_genomic_span(fact_context)
    transcript = fact_context.overlap.feature
    return genomic_overlap(variant_min, variant_max, transcript.start, transcript.end)


@register_fact(CostTier.BOUNDS)
def contain_feature(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the variant contains the feature's genomic span.

    Not zone-mappable: a zone-map union answers "does the span touch this interval", not
    "is the span fully inside one specific interval".
    """
    variant_min, variant_max = variant_genomic_span(fact_context)
    transcript = fact_context.overlap.feature
    return variant_min <= transcript.start and variant_max >= transcript.end


@register_fact(CostTier.BOUNDS)
def within_feature(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the variant lies within the feature's genomic span.

    Not zone-mappable.
    """
    variant_min, variant_max = variant_genomic_span(fact_context)
    transcript = fact_context.overlap.feature
    return transcript.start <= variant_min and variant_max <= transcript.end


@register_fact(CostTier.BOUNDS)
def is_upstream(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the variant lies within the upstream pairing distance."""
    variant_min, variant_max = variant_genomic_span(fact_context)
    transcript = fact_context.overlap.feature
    distance = upstream_distance(fact_context)
    if transcript.strand == Strand.FORWARD:
        return 0 < transcript.start - variant_max <= distance
    return 0 < variant_min - transcript.end <= distance


@register_fact(CostTier.BOUNDS)
def is_downstream(fact_context: FactContext[FeatureOverlap]) -> bool:
    """Return whether the variant lies within the downstream pairing distance."""
    variant_min, variant_max = variant_genomic_span(fact_context)
    transcript = fact_context.overlap.feature
    distance = downstream_distance(fact_context)
    if transcript.strand == Strand.FORWARD:
        return 0 < variant_min - transcript.end <= distance
    return 0 < transcript.start - variant_max <= distance


# --------------------------------------------------------------------------
# Placement within a transcript
# --------------------------------------------------------------------------


@register_fact(CostTier.INTERVAL, requires=("is_transcript", "overlap_feature"))
def overlap_exon(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant overlaps an exon."""
    variant_min, variant_max = variant_genomic_span(fact_context)
    return any(
        genomic_overlap(variant_min, variant_max, exon.start, exon.end)
        for exon in fact_context.overlap.feature.spliced_exons
    )


@register_fact(CostTier.INTERVAL, requires=("is_transcript", "overlap_feature"))
def within_any_exon(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant lies within one exon.

    Not zone-mappable.
    """
    variant_min, variant_max = variant_genomic_span(fact_context)
    return any(
        exon.start <= variant_min and variant_max <= exon.end
        for exon in fact_context.overlap.feature.spliced_exons
    )


@register_fact(CostTier.INTERVAL, requires=("is_transcript", "overlap_feature"))
def overlap_intron(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant overlaps an intron."""
    variant_min, variant_max = variant_genomic_span(fact_context)
    return any(
        genomic_overlap(variant_min, variant_max, intron.start, intron.end)
        for intron in fact_context.overlap.feature.introns
    )


@register_fact(CostTier.INTERVAL, requires=("is_transcript", "overlap_feature"))
def overlap_intron_interior(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant overlaps an intron away from its splice sites.

    VEP's ``intronic`` test (``BaseTranscriptVariationAllele::_intron_effects``): overlap with
    ``intron.start + 2 .. intron.end - 2``. An insertion's span is its two flanking bases, so one
    at either end of that window counts, as VEP's own insertion special case has it.
    """
    variant_min, variant_max = variant_genomic_span(fact_context)
    return any(
        genomic_overlap(variant_min, variant_max, intron.start + 2, intron.end - 2)
        for intron in fact_context.overlap.feature.introns
    )


@register_fact(CostTier.INTERVAL, requires=("is_transcript", "overlap_feature"))
def overlap_intron_boundary(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant overlaps an intron boundary window."""
    variant_min, variant_max = variant_genomic_span(fact_context)
    return any(
        genomic_overlap(variant_min, variant_max, intron.start - 3, intron.start + 7)
        or genomic_overlap(variant_min, variant_max, intron.end - 7, intron.end + 3)
        for intron in fact_context.overlap.feature.introns
    )


@register_fact(CostTier.INTERVAL)
def _nearby_introns(fact_context: FactContext[TranscriptVariationAllele]) -> tuple[Intron, ...]:
    """Return every intron within 20 bases of the variant."""
    variant_min, variant_max = variant_genomic_span(fact_context)
    transcript = fact_context.overlap.feature
    return tuple(
        intron
        for intron in transcript.introns
        if genomic_overlap(variant_min - 20, variant_max + 20, intron.start, intron.end)
    )


@register_fact(
    CostTier.INTERVAL,
    deps=("_nearby_introns",),
    requires=("is_transcript", "overlap_feature", "overlap_intron_boundary"),
)
def overlap_donor_site(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant overlaps donor-site bases +1 or +2."""
    return intron_offset_hit(fact_context, 1, 2, from_donor=True)


@register_fact(
    CostTier.INTERVAL,
    deps=("_nearby_introns",),
    requires=("is_transcript", "overlap_feature", "overlap_intron_boundary"),
)
def overlap_acceptor_site(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant overlaps acceptor-site bases -1 or -2."""
    return intron_offset_hit(fact_context, 1, 2, from_donor=False)


@register_fact(
    CostTier.INTERVAL,
    deps=("_nearby_introns",),
    requires=("is_transcript", "overlap_feature", "overlap_intron_boundary"),
)
def overlap_donor_5th_base(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant overlaps donor base +5."""
    return intron_offset_hit(fact_context, 5, 5, from_donor=True)


@register_fact(
    CostTier.INTERVAL,
    deps=("_nearby_introns",),
    requires=("is_transcript", "overlap_feature", "overlap_intron_boundary"),
)
def overlap_donor_region(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant overlaps donor-region bases +3 to +6."""
    return intron_offset_hit(fact_context, 3, 6, from_donor=True)


@register_fact(
    CostTier.INTERVAL,
    deps=("_nearby_introns",),
    requires=("is_transcript", "overlap_feature", "overlap_intron"),
)
def overlap_polypyrimidine_tract(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant overlaps polypyrimidine-tract bases -3 to -17."""
    return intron_offset_hit(fact_context, 3, 17, from_donor=False)


@register_fact(
    CostTier.INTERVAL,
    deps=("_nearby_introns",),
    requires=("is_transcript", "overlap_feature", "overlap_intron_boundary"),
)
def overlap_splice_region(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant overlaps an exonic or intronic splice region."""
    if intron_offset_hit(fact_context, 3, 8, from_donor=True):
        return True
    if intron_offset_hit(fact_context, 3, 8, from_donor=False):
        return True
    variant_min, variant_max = variant_genomic_span(fact_context)
    for intron in fact_context.overlap.feature.introns:
        if genomic_overlap(variant_min, variant_max, intron.start - 3, intron.start - 1):
            return True
        if genomic_overlap(variant_min, variant_max, intron.end + 1, intron.end + 3):
            return True
    return False


# --------------------------------------------------------------------------
# Placement within the coding sequence (CDS) - requires mapping genomic-to-CDS coordinates.
# --------------------------------------------------------------------------


@register_fact(CostTier.MAPPING)
def _cds_footprint(fact_context: FactContext[TranscriptVariationAllele]) -> tuple[int, int] | None:
    """Return the zero-based CDS offset span touched by the variant.

    Delegates to Mapper.genomic_span_to_cds_offsets, which intersects the variant span against the
    CDS's genomic spans directly (O(exon count)) rather than mapping every individual genomic
    position in the variant's span (O(span length)) which can be very large for a large structural variant.

    No accelerator needed here: the block-overlap approach is already O(exon count), which the O(1)
    per-position MapperIndex cannot improve on further.

    NOTE: It is important to check overlap with every coding span, not just the two endpoints of the variant's
    genomic span. A variant can have both endpoints outside of coding exons but still overlap the CDS.
    """
    transcript = fact_context.overlap.feature
    variant_min, variant_max = variant_genomic_span(fact_context)
    mapper = transcript.cds_mapper
    if mapper is None:
        return None
    return mapper.genomic_span_to_offsets(variant_min, variant_max)


@register_fact(
    CostTier.MAPPING,
    deps=("_cds_footprint",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_exon",
        "has_translation",
    ),
)
def overlap_cds(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant overlaps the coding sequence."""
    if fact_context.overlap.feature.translation is None:
        return False
    return get_fact(fact_context, "_cds_footprint") is not None


@register_fact(
    CostTier.INTERVAL,
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_exon",
    ),
)
def overlap_five_prime_utr(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant overlaps the transcript's 5' UTR.

    Tested against genomic spans of the five_prime_utr instead of just the span(start, end) as a UTR can
    contain introns (non-UTR region) and therefore, can give inaccurate results.
    """
    spans = fact_context.overlap.feature.five_prime_utr_genomic_spans
    if not spans:
        return False
    variant_min, variant_max = variant_genomic_span(fact_context)
    return any(genomic_overlap(variant_min, variant_max, start, end) for start, end in spans)


@register_fact(
    CostTier.INTERVAL,
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_exon",
    ),
)
def overlap_three_prime_utr(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant overlaps the transcript's 3' UTR."""
    spans = fact_context.overlap.feature.three_prime_utr_genomic_spans
    if not spans:
        return False
    variant_min, variant_max = variant_genomic_span(fact_context)
    return any(genomic_overlap(variant_min, variant_max, start, end) for start, end in spans)


@register_fact(
    CostTier.MAPPING,
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_exon",
        "has_translation",
    ),
)
def overlap_start_codon(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant overlaps the start codon.

    Similar to UTRs a start codon can split across an exon-exon junction and therefore should be tested
    against all its genomic spans.
    """
    codon_spans = fact_context.overlap.feature.start_codon_genomic_spans()
    if codon_spans is None:
        return False
    variant_min, variant_max = variant_genomic_span(fact_context)
    return any(
        genomic_overlap(variant_min, variant_max, span_start, span_end)
        for span_start, span_end in codon_spans
    )


@register_fact(
    CostTier.MAPPING,
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_exon",
        "overlap_cds",
    ),
)
def overlap_stop_codon(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant overlaps the stop codon."""
    codon_spans = fact_context.overlap.feature.stop_codon_genomic_spans()
    if codon_spans is None:
        return False
    variant_min, variant_max = variant_genomic_span(fact_context)
    return any(
        genomic_overlap(variant_min, variant_max, span_start, span_end)
        for span_start, span_end in codon_spans
    )


@register_fact(
    CostTier.MAPPING,
    deps=("_cds_footprint",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_exon",
    ),
)
def partial_last_codon(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Return whether the variant falls in an incomplete terminal codon."""
    footprint = cast("tuple[int, int] | None", get_fact(fact_context, "_cds_footprint"))
    if footprint is None:
        return False
    _, high_offset = footprint

    cds_length = fact_context.overlap.feature.cds_length
    if cds_length is None:
        return False

    # The last codon that variant overlaps (high_offset) extends beyond the end of the CDS, meaning that
    # the last codon is partial.
    return (high_offset // 3) * 3 + 3 > cds_length
