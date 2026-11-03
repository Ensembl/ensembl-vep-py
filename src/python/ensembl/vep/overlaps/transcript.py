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
"""Transcript-overlap finding."""

from __future__ import annotations

from attrs import define, field

from ensembl.vep.caches.interval_tree import Interval
from ensembl.vep.common import DEFAULT_DISTANCE, Readthrough, Strand
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptStructuralVariationAllele, TranscriptVariationAllele
from ensembl.vep.models.variation import StructuralVariationAllele, VariationAllele
from ensembl.vep.overlaps.base import CandidateLookup, FeatureMatcher
from ensembl.vep.overlaps.feature import FeatureOverlapFinder

__all__ = ["TranscriptOverlapFinder", "transcript_overlap"]


def transcript_overlap(
    transcript: Transcript,
    variation_allele: VariationAllele,
    *,
    readthrough: Readthrough = Readthrough.NONE,
) -> TranscriptVariationAllele:
    """Construct the appropriate transcript overlap model."""
    overlap_class = (
        TranscriptStructuralVariationAllele
        if isinstance(variation_allele, StructuralVariationAllele)
        else TranscriptVariationAllele
    )
    return overlap_class(
        feature=transcript,
        overlapping_feature=variation_allele,
        readthrough=readthrough,
    )


@define(kw_only=True, auto_attribs=True)
class TranscriptOverlapFinder(FeatureOverlapFinder[Transcript, TranscriptVariationAllele]):
    """Transcript-overlap finder for a variation-allele buffer.

    Attributes:
        upstream_distance (int): Genomic bases considered before a transcript's 5' end.
        downstream_distance (int): Genomic bases considered past a transcript's 3' end.
    """

    candidates: CandidateLookup[Transcript]
    matcher: FeatureMatcher
    readthrough: Readthrough = field(default=Readthrough.NONE)
    upstream_distance: int = field(default=DEFAULT_DISTANCE)
    downstream_distance: int = field(default=DEFAULT_DISTANCE)

    def __attrs_post_init__(self) -> None:
        for name in ("upstream_distance", "downstream_distance"):
            distance = getattr(self, name)
            if isinstance(distance, bool) or not isinstance(distance, int):
                raise TypeError(f"Transcript overlap {name} must be an integer.")
            if distance < 0:
                raise ValueError(f"Transcript overlap {name} must not be negative.")
        self.overlap_factory = self._overlap_factory
        self.candidate_intervals = self._candidate_intervals
        self.match_subject = self._match_subject

    def _candidate_intervals(self, transcript: Transcript) -> tuple[Interval, ...]:
        """Return the transcript span expanded for flanking consequences.

        Which side each distance applies to follows the transcript: a reverse-strand transcript's
        5' end is its highest coordinate, so the upstream distance extends past ``end`` there.
        """
        before, after = self.upstream_distance, self.downstream_distance
        if transcript.strand == Strand.REVERSE:
            before, after = after, before
        return (
            Interval(
                start=max(1, transcript.start - before),
                end=transcript.end + after,
            ),
        )

    def _match_subject(self, transcript: Transcript) -> Interval:
        """Return the transcript influence interval used for final association."""
        return self._candidate_intervals(transcript)[0]

    def _overlap_factory(
        self,
        transcript: Transcript,
        variation_allele: VariationAllele,
    ) -> TranscriptVariationAllele:
        """Construct a transcript overlap using this finder's readthrough policy."""
        return transcript_overlap(transcript, variation_allele, readthrough=self.readthrough)
