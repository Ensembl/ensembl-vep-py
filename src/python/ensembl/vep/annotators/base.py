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
"""Annotator protocol shared by all annotation domains.

``Annotator`` is the contract callers depend on. ``registry.get_annotator`` returns one without
requiring inheritance, allowing Python and native implementations to share the registry. The
protocol exposes single-overlap and batch annotation; implementations may derive either operation
from the other.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Protocol, TypeVar, runtime_checkable

from ensembl.vep.models.overlap import FeatureOverlap

__all__ = ["Annotator", "ResultT"]

ResultT = TypeVar("ResultT")
"""What an annotator produces per overlap."""


@runtime_checkable
class Annotator(Protocol[ResultT]):
    """Define the common contract for annotators producing ``ResultT``.

    Attributes:
        name (str): The registry name this annotator is registered under, e.g. "consequence.linear".
    """

    name: str

    def annotate_feature_overlap(self, overlap: FeatureOverlap) -> list[ResultT]:
        """Annotate a feature overlap.

        Args:
            overlap (FeatureOverlap): The overlap to annotate.

        Returns:
            list[ResultT]: What this annotator produces for the overlap. May be empty -- for a
                consequence annotator that means callers apply the default/fallback term (e.g.
                ``intergenic_variant``) themselves rather than the annotator doing so.
        """
        raise NotImplementedError

    def annotate_feature_overlaps(self, overlaps: Sequence[FeatureOverlap]) -> list[list[ResultT]]:
        """Annotate a batch of feature overlaps at once.

        Args:
            overlaps (Sequence[FeatureOverlap]): The overlaps to annotate, in order.

        Returns:
            list[list[ResultT]]: One result list per input overlap, in the same order as
                ``overlaps``.
        """
        raise NotImplementedError
