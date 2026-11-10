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
"""Tests for overlap finders."""

from __future__ import annotations

from collections.abc import Iterable
from typing import cast

import pytest

from ensembl.vep.caches.interval_tree import Interval
from ensembl.vep.common import Strand
from ensembl.vep.matcher.base import IntervalLike
from ensembl.vep.models import Feature
from ensembl.vep.models.custom_feature import CustomFeature
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import (
    CustomVariationAllele,
    FeatureOverlap,
    TranscriptStructuralVariationAllele,
    TranscriptVariationAllele,
)
from ensembl.vep.models.regulatory import Regulatory
from ensembl.vep.models.variation import StructuralVariationAllele, Variation, VariationAllele
from ensembl.vep.overlaps import (
    ColocatedVariantOverlapFinder,
    CustomOverlapFinder,
    FeatureOverlapFinder,
    RegulatoryOverlapFinder,
    TranscriptOverlapFinder,
)
from ensembl.vep.overlaps.base import CandidateLookup, VariationAlleleBuffer


class _Buffer:
    """Minimal variation-allele buffer for finder tests."""

    def __init__(self, objects: Iterable[VariationAllele]) -> None:
        self.objects = tuple(objects)
        self.queries: list[tuple[Interval, ...]] = []

    def get_objects_by_intervals(self, intervals: Iterable[Interval]) -> list[VariationAllele]:
        """Return alleles inside one or more inclusive query intervals."""
        query_intervals = tuple(intervals)
        self.queries.append(query_intervals)
        return [
            variation_allele
            for variation_allele in self.objects
            if any(interval.start <= variation_allele.start <= interval.end for interval in query_intervals)
        ]


class _ContainingMatcher:
    """Match a variation allele contained by a feature."""

    def matches(self, feature: IntervalLike, variation_allele: IntervalLike) -> bool:
        """Return whether the allele lies within the feature."""
        return feature.start <= variation_allele.start <= feature.end


def _features(_input_buffer: VariationAlleleBuffer) -> list[Feature]:
    """Return two candidate features."""
    return [
        Feature(start=1, end=10, strand=Strand.FORWARD),
        Feature(start=20, end=30, strand=Strand.FORWARD),
    ]


class TestFeatureOverlapFinder:
    """Candidate pairing and overlap construction."""

    def test_yields_only_matched_feature_pairs(self) -> None:
        """Apply the matcher before constructing generic overlaps."""
        variation_allele = VariationAllele(start=5, end=5, strand=Strand.FORWARD)
        finder: FeatureOverlapFinder[Feature, FeatureOverlap] = FeatureOverlapFinder(
            candidates=cast(CandidateLookup[Feature], _features),
            matcher=_ContainingMatcher(),
            overlap_factory=lambda feature, allele: FeatureOverlap(
                feature=feature,
                overlapping_feature=allele,
            ),
        )

        overlaps = list(finder.find(_Buffer([variation_allele])))

        assert len(overlaps) == 1
        assert overlaps[0].feature.start == 1
        assert overlaps[0].overlapping_feature is variation_allele

    def test_looks_up_candidates_once_per_buffer(self) -> None:
        """Use indexed lookup once for each candidate feature."""
        calls = 0

        def candidates(input_buffer: VariationAlleleBuffer) -> list[Feature]:
            nonlocal calls
            del input_buffer
            calls += 1
            return _features(_Buffer(()))

        finder: FeatureOverlapFinder[Feature, FeatureOverlap] = FeatureOverlapFinder(
            candidates=cast(CandidateLookup[Feature], candidates),
            matcher=_ContainingMatcher(),
            overlap_factory=lambda feature, allele: FeatureOverlap(
                feature=feature,
                overlapping_feature=allele,
            ),
        )

        input_buffer = _Buffer(
            [
                VariationAllele(start=5, end=5, strand=Strand.FORWARD),
                VariationAllele(start=25, end=25, strand=Strand.FORWARD),
            ]
        )
        overlaps = list(finder.find(input_buffer))

        assert calls == 1
        assert len(overlaps) == 2
        assert input_buffer.queries == [
            (Interval(start=1, end=10),),
            (Interval(start=20, end=30),),
        ]


class TestTranscriptOverlapFinder:
    """Transcript-specific overlap construction."""

    @pytest.fixture
    def transcript(self) -> Transcript:
        """Return a simple transcript candidate."""
        return Transcript(stable_id="ENST00000000001", start=1, end=10, strand=Strand.FORWARD)

    def test_constructs_a_transcript_variation_allele(self, transcript: Transcript) -> None:
        """Use the transcript model for a small variation allele."""
        finder = TranscriptOverlapFinder(
            candidates=lambda input_buffer: [transcript],
            matcher=_ContainingMatcher(),
        )

        overlap = next(finder.find(_Buffer([VariationAllele(start=5, end=5, strand=Strand.FORWARD)])))

        assert isinstance(overlap, TranscriptVariationAllele)
        assert overlap.feature is transcript

    def test_constructs_a_structural_transcript_overlap(self, transcript: Transcript) -> None:
        """Use the structural overlap model for a structural variation allele."""
        finder = TranscriptOverlapFinder(
            candidates=lambda input_buffer: [transcript],
            matcher=_ContainingMatcher(),
        )

        overlap = next(
            finder.find(_Buffer([StructuralVariationAllele(start=5, end=5, strand=Strand.FORWARD)]))
        )

        assert isinstance(overlap, TranscriptStructuralVariationAllele)

    def test_queries_the_configured_flank(self, transcript: Transcript) -> None:
        """Match an upstream allele through the transcript influence interval."""
        transcript = Transcript(
            stable_id=transcript.stable_id,
            start=10,
            end=20,
            strand=Strand.FORWARD,
        )
        finder = TranscriptOverlapFinder(
            candidates=lambda input_buffer: [transcript],
            matcher=_ContainingMatcher(),
            upstream_distance=5,
            downstream_distance=5,
        )
        input_buffer = _Buffer([VariationAllele(start=5, end=5, strand=Strand.FORWARD)])

        overlaps = list(finder.find(input_buffer))

        assert len(overlaps) == 1
        assert overlaps[0].feature is transcript
        assert input_buffer.queries == [(Interval(start=5, end=25),)]

    @pytest.mark.parametrize("strand", [Strand.FORWARD, Strand.REVERSE])
    @pytest.mark.parametrize(
        ("distance", "expected"),
        [
            pytest.param(4999, True, id="inside-flank"),
            pytest.param(5000, True, id="flank-boundary"),
            pytest.param(5001, False, id="outside-flank"),
        ],
    )
    def test_matches_only_within_the_configured_flank(
        self,
        strand: Strand,
        distance: int,
        expected: bool,
    ) -> None:
        """Use the full configured flank on either transcript strand."""
        transcript = Transcript(
            stable_id="ENST00000000001",
            start=10000,
            end=20000,
            strand=strand,
        )
        finder = TranscriptOverlapFinder(
            candidates=lambda input_buffer: [transcript],
            matcher=_ContainingMatcher(),
        )
        variation_allele = VariationAllele(
            start=transcript.end + distance,
            end=transcript.end + distance,
            strand=Strand.FORWARD,
        )

        overlaps = list(finder.find(_Buffer([variation_allele])))

        assert bool(overlaps) is expected


class TestAsymmetricDistances:
    """``--distance bp,downstream``: each side measured from the end it belongs to."""

    @staticmethod
    def _interval(transcript: Transcript, upstream: int, downstream: int) -> Interval:
        """Return the influence interval one pair of distances gives a transcript."""
        finder = TranscriptOverlapFinder(
            candidates=lambda _buffer: [transcript],
            matcher=_ContainingMatcher(),
            upstream_distance=upstream,
            downstream_distance=downstream,
        )
        # pylint: disable=protected-access
        return finder._candidate_intervals(transcript)[0]

    def test_a_forward_transcript_extends_upstream_before_its_start(self) -> None:
        """5' is the lower coordinate on the forward strand."""
        transcript = Transcript(stable_id="ENST1", start=1000, end=2000, strand=Strand.FORWARD)

        interval = self._interval(transcript, 500, 100)

        assert (interval.start, interval.end) == (500, 2100)

    def test_a_reverse_transcript_extends_upstream_past_its_end(self) -> None:
        """5' is the higher coordinate on the reverse strand, so the distances swap sides."""
        transcript = Transcript(stable_id="ENST2", start=1000, end=2000, strand=Strand.REVERSE)

        interval = self._interval(transcript, 500, 100)

        assert (interval.start, interval.end) == (900, 2500)

    def test_one_distance_still_applies_both_ways(self) -> None:
        """The single-value form is the symmetric window it has always been."""
        transcript = Transcript(stable_id="ENST3", start=1000, end=2000, strand=Strand.FORWARD)

        interval = self._interval(transcript, 300, 300)

        assert (interval.start, interval.end) == (700, 2300)

    def test_the_window_is_clamped_to_the_first_base(self) -> None:
        """A distance wider than the transcript's start cannot run off the region."""
        transcript = Transcript(stable_id="ENST4", start=10, end=20, strand=Strand.FORWARD)

        assert self._interval(transcript, 5000, 5000).start == 1

    def test_a_downstream_only_allele_is_dropped_by_a_zero_distance(self) -> None:
        """Setting one side to zero pairs on that side with nothing."""
        transcript = Transcript(stable_id="ENST5", start=1000, end=2000, strand=Strand.FORWARD)
        finder = TranscriptOverlapFinder(
            candidates=lambda _buffer: [transcript],
            matcher=_ContainingMatcher(),
            upstream_distance=500,
            downstream_distance=0,
        )

        downstream_allele = VariationAllele(start=2050, end=2050, strand=Strand.FORWARD)
        upstream_allele = VariationAllele(start=600, end=600, strand=Strand.FORWARD)

        assert not list(finder.find(_Buffer([downstream_allele])))
        assert list(finder.find(_Buffer([upstream_allele])))


class TestOtherOverlapFinders:
    """Feature kinds using generic overlap models."""

    @pytest.mark.parametrize(
        ("finder", "feature"),
        [
            pytest.param(
                RegulatoryOverlapFinder,
                Regulatory(start=1, end=10, strand=Strand.FORWARD),
                id="regulatory",
            ),
            pytest.param(
                CustomOverlapFinder,
                CustomFeature(start=1, end=10, strand=Strand.FORWARD),
                id="custom",
            ),
            pytest.param(
                ColocatedVariantOverlapFinder,
                Variation(start=1, end=10, strand=Strand.FORWARD),
                id="colocated",
            ),
        ],
    )
    def test_constructs_a_generic_feature_overlap(
        self,
        finder: type[FeatureOverlapFinder[Feature, FeatureOverlap]],
        feature: Feature,
    ) -> None:
        """Construct the shared overlap model for non-transcript features."""
        overlap_finder = finder(
            candidates=lambda input_buffer: [feature],
            matcher=_ContainingMatcher(),
        )

        overlap = next(overlap_finder.find(_Buffer([VariationAllele(start=5, end=5, strand=Strand.FORWARD)])))

        assert isinstance(overlap, FeatureOverlap)
        assert overlap.feature is feature

    def test_custom_finder_constructs_a_typed_custom_overlap(self) -> None:
        """Retain custom record geometry in its dedicated overlap model."""
        feature = CustomFeature(start=1, end=10, strand=Strand.FORWARD, record={"score": 2})
        finder = CustomOverlapFinder(
            candidates=lambda input_buffer: [feature],
            matcher=_ContainingMatcher(),
        )

        overlap = next(finder.find(_Buffer([VariationAllele(start=5, end=5, strand=Strand.FORWARD)])))

        assert isinstance(overlap, CustomVariationAllele)
        assert overlap.feature.record == {"score": 2}
