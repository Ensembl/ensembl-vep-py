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
"""Annotation-track contracts and value objects.

A track owns everything one kind of annotation needs: the source it reads, whatever it caches
between payloads, the rule pairing features with alleles, the annotator interpreting each pair,
and where the results are attached for output. The worker running a track supplies only what
every track shares -- the reference sequence and the synonym map -- so adding a kind of
annotation adds a track and changes nothing above it.
"""

from __future__ import annotations

from collections.abc import Sequence
from itertools import islice
from typing import Generic, Protocol

from attrs import define, field

from ensembl.vep.annotators.base import Annotator, ResultT
from ensembl.vep.config.sections.common.source import SourceConfig
from ensembl.vep.models.overlap import FeatureOverlap
from ensembl.vep.models.variation import Variation
from ensembl.vep.overlaps import OverlapFinder
from ensembl.vep.overlaps.base import VariationAlleleBuffer
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource
from ensembl.vep.sources.synonym import ChromosomeSynonyms

__all__ = ["AnnotationTrack", "OverlapAnnotationTrack", "TrackResources", "TrackResult"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class TrackResources:
    """Resources shared by every track in one workflow.

    Built once per worker and passed to each track builder: a source reads reference sequence and
    resolves sequence-region names through the same objects whichever track built it.

    Attributes:
        source_config (SourceConfig): The configured sources, for the files and distances a track reads.
        sequence_source (FastaSequenceSource | None): The shared reference source, or None for a
            caller that can annotate without reference sequence.
        synonyms (ChromosomeSynonyms): The run's sequence-region synonyms.
    """

    source_config: SourceConfig
    sequence_source: FastaSequenceSource | None
    synonyms: ChromosomeSynonyms


@define(frozen=True, kw_only=True, auto_attribs=True)
class TrackResult(Generic[ResultT]):
    """Annotation results for one overlap.

    Attributes:
        track_identifier (str): Identifier of the producing track.
        overlap (FeatureOverlap): The overlap passed to the annotator.
        annotations (tuple[ResultT, ...]): Values emitted for the overlap.
    """

    track_identifier: str
    overlap: FeatureOverlap
    annotations: tuple[ResultT, ...]


class AnnotationTrack(Protocol[ResultT]):
    """Configured annotation concern for one input buffer."""

    @property
    def identifier(self) -> str:
        """Return the unique configured track identifier."""
        raise NotImplementedError

    def begin_payload(self, input_buffer: VariationAlleleBuffer) -> None:
        """Release the previous payload's state and prime this one's."""
        raise NotImplementedError

    def annotate(self, input_buffer: VariationAlleleBuffer) -> list[TrackResult[ResultT]]:
        """Return annotations for the input buffer."""
        raise NotImplementedError

    def attach(
        self,
        variations: Sequence[Variation],
        results: Sequence[TrackResult[ResultT]],
    ) -> None:
        """Attach this track's results to the alleles output reads them from.

        Every allele of every variation is assigned, so attaching no results releases what a
        rendered window retained.
        """
        raise NotImplementedError

    def close(self) -> None:
        """Release the source and caches this track owns."""
        raise NotImplementedError


@define(kw_only=True, auto_attribs=True)
class OverlapAnnotationTrack(Generic[ResultT]):
    """Overlap finder and annotator for one configured concern.

    Specialised per kind of annotation: a subclass holds the source and caches its finder reads
    through, and attaches the results where that kind of overlap belongs.

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

    def begin_payload(self, input_buffer: VariationAlleleBuffer) -> None:
        """Prime nothing: a track caching between payloads overrides this."""

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

    def attach(
        self,
        variations: Sequence[Variation],
        results: Sequence[TrackResult[ResultT]],
    ) -> None:
        """Reject results no subclass claimed, rather than dropping them silently.

        Raises:
            NotImplementedError: Always. A track that produces results owns where they go.
        """
        raise NotImplementedError(f"Track '{self.identifier}' does not attach its results.")

    def close(self) -> None:
        """Release nothing: a track owning a source overrides this."""
