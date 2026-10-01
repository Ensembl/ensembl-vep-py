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
"""Generic feature-overlap finding."""

from __future__ import annotations

from collections.abc import Iterator
from typing import Generic

from attrs import define

from ensembl.vep.caches.interval_index import Interval
from ensembl.vep.models import Feature
from ensembl.vep.models.overlap import FeatureOverlap
from ensembl.vep.models.variation import VariationAllele
from ensembl.vep.overlaps.base import (
    CandidateLookup,
    CandidateIntervals,
    FeatureMatcher,
    FeatureT,
    MatchSubject,
    OverlapFactory,
    OverlapT,
    VariationAlleleBuffer,
)

__all__ = ["FeatureOverlapFinder"]


def _feature_overlap(feature: Feature, variation_allele: VariationAllele) -> FeatureOverlap:
    """Construct a generic feature overlap."""
    return FeatureOverlap(feature=feature, overlapping_feature=variation_allele)


def _feature_intervals(feature: Feature) -> tuple[Interval, ...]:
    """Return the feature's own interval as one buffer query."""
    return (Interval(start=feature.start, end=feature.end),)


def _match_feature(feature: Feature) -> Feature:
    """Return the feature itself as the final matcher subject."""
    return feature


@define(kw_only=True, auto_attribs=True)
class FeatureOverlapFinder(Generic[FeatureT, OverlapT]):
    """Find overlaps between source features and variation alleles.

    Candidate lookup owns source, cache, and index selection. This keeps the finder independent of
    the feature source and cache implementations used by a workflow track.

    Attributes:
        candidates (CandidateLookup[FeatureT]): Fetch candidates for each input buffer.
        candidate_intervals (CandidateIntervals[FeatureT]): Build buffer-index queries for a feature.
        match_subject (MatchSubject[FeatureT]): Build the final matcher subject for a feature.
        matcher (FeatureMatcher): Decide whether each candidate pair is associated.
        overlap_factory (OverlapFactory[FeatureT, OverlapT]): Construct matched overlap models.
    """

    candidates: CandidateLookup[FeatureT]
    matcher: FeatureMatcher
    candidate_intervals: CandidateIntervals[FeatureT] = _feature_intervals
    match_subject: MatchSubject[FeatureT] = _match_feature
    overlap_factory: OverlapFactory[FeatureT, OverlapT] = _feature_overlap  # type: ignore[assignment]

    def find(self, input_buffer: VariationAlleleBuffer) -> Iterator[OverlapT]:
        """Yield overlaps for every matching candidate pair."""
        for feature in self.candidates(input_buffer):
            for variation_allele in input_buffer.get_objects_by_intervals(self.candidate_intervals(feature)):
                if self.matcher.matches(self.match_subject(feature), variation_allele):
                    yield self.overlap_factory(feature, variation_allele)
