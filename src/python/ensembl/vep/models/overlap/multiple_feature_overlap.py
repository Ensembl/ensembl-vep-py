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