from __future__ import annotations

from attr import define, field

from ensembl.vep.models.overlap import FeatureOverlap


class TranscriptStructuralVariationAllele(FeatureOverlap):
    """
    Represents the overlap between a transcript and a structural variation allele.
    """
    ...
    exon_numbers: list[int] = field(default=[], init=False)
    intron_numbers: list[int] = field(default=[], init=False)
