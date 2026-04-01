from __future__ import annotations

from attr import define, field
from typing import TYPE_CHECKING
from functools import cached_property
from threading import Lock

from ensembl.vep.models import Feature
if TYPE_CHECKING:
    from ensembl.vep.models.overlap import OverlapConsequence


@define(kw_only=True, auto_attribs=True)
class MultipleFeatureOverlap:
    """
    Represents an overlap between a feature and multiple overlapping features, along with the consequences of that overlap.
    
    Attributes:
        features: The features that are getting overlapped.
        overlapping_feature: The feature that is overlapping.
        overlap_consequences: A list of consequences the overlapping feature has on the feature.
    """
    feature: Feature
    overlapping_features: list[Feature]
    overlap_consequences: list[OverlapConsequence] = field(default=[], init=False)
    pass