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
"""MultipleFeatureOverlap class."""

from __future__ import annotations

from typing import TYPE_CHECKING
from attr import define, field

from ensembl.vep.models import Feature

if TYPE_CHECKING:
    from ensembl.vep.models.overlap import OverlapConsequence


@define(kw_only=True, auto_attribs=True)

__all__ = ["MultipleFeatureOverlap"]


class MultipleFeatureOverlap:
    """
    Represents an overlap between a feature and multiple overlapping features,
    along with the consequences of that overlap.

    Attributes:
        feature (Feature):
            The feature that is getting overlapped.
        overlapping_features (list[Feature]):
            The features that are overlapping.
        overlap_consequences (list[OverlapConsequence]):
            A list of consequences the overlapping features have on the feature.
    """

    feature: Feature
    overlapping_features: list[Feature]
    overlap_consequences: list[OverlapConsequence] = field(factory=list, init=False)
    # more ...
