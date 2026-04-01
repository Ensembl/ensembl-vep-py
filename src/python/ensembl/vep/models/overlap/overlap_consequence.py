from __future__ import annotations

from attrs import define, field

from ensembl.vep.models import Feature

@define(kw_only=True, auto_attribs=True)
class OverlapConsequence:
    """
    A predicted consequence of an overlap between a feature and a base feature. 
    
    This is typically used to represent the consequence of a variant on a transcript, but it can be used 
    for any type of feature overlap.

    Attributes:
        term: The SO term for the consequence, e.g. "missense_variant".
        accession: The SO accession for the consequence, e.g. "SO:0001583".
        is_default: Whether this consequence is the default consequence for the overlapping feature.
            Only one consequence can be default, for transcript variation allele it is generally intergenic_variant.
        rank: The rank of the consequence, with 1 being the most severe.
        tier: The tier of the consequence. If tier 1 consequence are assigned, then no tier 2 consequences will be
            checked or assigned.
        impact: The predicted impact of the consequence, e.g. "HIGH", "MODERATE", "LOW", "MODIFIER".
        feature_type: The type of the overlapped feature that is affected by this consequence type, e.g. Transcript.
        overlapping_feature_type: The type of the overlapping feature that causes this consequence type, e.g. VariationAllele.
    """
    term: str
    accession: str
    is_default: bool = field(default=False, init=False)
    rank: int = field(default=0, init=False)
    tier: int = field(default=0, init=False)
    impact: str | None = field(default=None, init=False)
    feature_type: type[Feature] | None = field(default=None, init=False)
    overlapping_feature_type: type[Feature] | None = field(default=None, init=False)