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
"""FeatureOverlap class."""

from __future__ import annotations

from typing import TYPE_CHECKING
from functools import cached_property
from threading import Lock
from attr import define, field

from ensembl.vep.models import Feature
from ensembl.vep.common import invalidates

if TYPE_CHECKING:
    from ensembl.vep.models.overlap import OverlapConsequence


@define(kw_only=True, auto_attribs=True)
class FeatureOverlap:
    """
    Represents an overlap between a feature and a overlapping feature,
    along with the consequences of that overlap.

    Attributes:
        feature (Feature):
            The feature that is getting overlapped.
        overlapping_feature (Feature):
            The feature that is overlapping.
        overlap_consequences (list[OverlapConsequence]):
            A list of consequences the overlapping feature has on the feature.
    """

    feature: Feature
    overlapping_feature: Feature
    overlap_consequences: list[OverlapConsequence] = field(
        factory=list, init=False, on_setattr=invalidates("__most_severe_consequence", "__consequence_terms")
    )
    _lock: Lock = field(factory=Lock, init=False, repr=False)

    @property
    def most_severe_consequence(self) -> OverlapConsequence | None:
        """
        Public accessor for __most_severe_consequence.

        Returns:
            OverlapConsequence | None: The most severe consequence, or None if there are no consequences.
        """
        return self.__most_severe_consequence

    @cached_property
    def __most_severe_consequence(self) -> OverlapConsequence | None:
        """
        Returns the most severe consequence, or None if there are no consequences.

        Returns:
            OverlapConsequence | None: The most severe consequence, or None if there are no consequences.
        """
        with self._lock:
            most_severe_consequence = None
            current_rank = -1
            for consequence in self.overlap_consequences:
                if consequence.rank > current_rank:
                    most_severe_consequence = consequence
                    current_rank = consequence.rank
            return most_severe_consequence

    @property
    def consequence_terms(self) -> set[str]:
        """
        Public accessor for __consequence_terms.

        Returns:
            set[str]: A set of all consequence terms.
        """
        return self.__consequence_terms

    @cached_property
    def __consequence_terms(self) -> set[str]:
        """
        Returns a set of all consequence terms.

        Returns:
            set[str]: A set of all consequence terms.
        """
        with self._lock:
            return set(consequence.term for consequence in self.overlap_consequences)
