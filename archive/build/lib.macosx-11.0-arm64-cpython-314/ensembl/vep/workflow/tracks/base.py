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
"""Annotation-track contracts and value objects."""

from __future__ import annotations

from collections.abc import Mapping
from itertools import islice
from typing import Generic, Protocol

from attrs import define, field

from ensembl.vep.annotators.base import Annotator, ResultT
from ensembl.vep.models.overlap import FeatureOverlap
from ensembl.vep.overlaps import OverlapFinder
from ensembl.vep.overlaps.base import VariationAlleleBuffer

__all__ = ["AnnotationTrack", "OverlapAnnotationTrack", "TrackResult", "TrackSpec"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class TrackSpec:
    """Describe one configured annotation track.

    Attributes:
        identifier (str): Unique name for this configured track instance.
        kind (str): Registered factory name for the track implementation.
        options (Mapping[str, object]): Track-specific construction values.
    """

    identifier: str
    kind: str
    options: Mapping[str, object] = field(factory=dict)


@define(frozen=True, kw_only=True, auto_attribs=True)
class TrackResult(Generic[ResultT]):
    """Hold annotation results for one overlap.

    Attributes:
        track_identifier (str): Identifier of the producing track.
        overlap (FeatureOverlap): The overlap passed to the annotator.
        annotations (tuple[ResultT, ...]): Values emitted for the overlap.
    """

    track_identifier: str
    overlap: FeatureOverlap
    annotations: tuple[ResultT, ...]


class AnnotationTrack(Protocol[ResultT]):
    """Annotate one input buffer for a configured concern."""

    @property
    def identifier(self) -> str:
        """Return the unique configured track identifier."""
        raise NotImplementedError

    def annotate(self, input_buffer: VariationAlleleBuffer) -> list[TrackResult[ResultT]]:
        """Return annotations for the input buffer."""
        raise NotImplementedError


@define(kw_only=True, auto_attribs=True)
class OverlapAnnotationTrack(Generic[ResultT]):
    """Run an overlap finder and annotator for one configured concern.

    Attributes:
        identifier (str): Unique name for this configured track instance.
        overlap_finder (OverlapFinder[FeatureOverlap]): Produce overlaps for an input buffer.
        annotator (Annotator[ResultT]): Interpret the overlaps found by this track.
        batch_size (int): Maximum overlaps materialised for one annotator call.
    """

    identifier: str
    overlap_finder: OverlapFinder[FeatureOverlap]
    annotator: Annotator[ResultT]
    batch_size: int = field(default=256)

    def __attrs_post_init__(self) -> None:
        if self.batch_size < 1:
            raise ValueError("Annotation track batch_size must be positive.")

    def annotate(self, input_buffer: VariationAlleleBuffer) -> list[TrackResult[ResultT]]:
        """Return annotations for every overlap found in ``input_buffer``."""
        results: list[TrackResult[ResultT]] = []
        overlap_iterator = self.overlap_finder.find(input_buffer)
        while overlaps := list(islice(overlap_iterator, self.batch_size)):
            annotation_groups = self.annotator.annotate_feature_overlaps(overlaps)
            if len(annotation_groups) != len(overlaps):
                raise ValueError(
                    f"Annotator '{self.annotator.name}' returned {len(annotation_groups)} result groups "
                    f"for {len(overlaps)} overlaps."
                )

            for overlap, annotations in zip(overlaps, annotation_groups, strict=True):
                results.append(
                    TrackResult(
                        track_identifier=self.identifier,
                        overlap=overlap,
                        annotations=tuple(annotations),
                    )
                )
        return results
