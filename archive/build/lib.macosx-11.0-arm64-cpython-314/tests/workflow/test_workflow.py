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
"""Tests for workflow composition."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from threading import Barrier
from typing import Any

import pytest

from ensembl.vep.annotators.base import Annotator
from ensembl.vep.common import Strand
from ensembl.vep.caches.interval_index import Interval
from ensembl.vep.matcher import OverlapMatcher
from ensembl.vep.models import Feature
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import FeatureOverlap
from ensembl.vep.models.variation import VariationAllele
from ensembl.vep.overlaps import FeatureOverlapFinder
from ensembl.vep.workflow import (
    AnnotationTrack,
    OverlapAnnotationTrack,
    TrackResult,
    TrackSpec,
    WorkflowPlan,
    WorkflowSpec,
    available_track_factories,
    available_workflows,
    build_track,
    build_workflow,
    register_track,
)
from ensembl.vep.workflow.tracks import registry as track_registry_module


class _Buffer:
    """Minimal variation-allele buffer for workflow tests."""

    def __init__(self, objects: Iterable[VariationAllele]) -> None:
        self.objects = tuple(objects)

    def get_objects_by_intervals(self, intervals: Iterable[Interval]) -> list[VariationAllele]:
        """Return every test allele for each query."""
        del intervals
        return list(self.objects)


class _AlwaysMatcher:
    """Match every candidate pair."""

    def matches(self, feature: Feature, variation_allele: VariationAllele) -> bool:
        """Return true for all pairs."""
        del feature, variation_allele
        return True


class _LabelsAnnotator(Annotator[str]):
    """Return an identifier for each overlap."""

    name = "test.labels"

    def annotate_feature_overlap(self, overlap: FeatureOverlap) -> list[str]:
        """Return the feature start as one label."""
        return [str(overlap.feature.start)]

    def annotate_feature_overlaps(self, overlaps: Sequence[FeatureOverlap]) -> list[list[str]]:
        """Annotate a sequence of overlaps."""
        return [self.annotate_feature_overlap(overlap) for overlap in overlaps]


class _StaticTrack:
    """Return preconstructed results without source or annotator dependencies."""

    def __init__(self, identifier: str, results: list[str]) -> None:
        self.identifier = identifier
        self._results = results

    def annotate(self, input_buffer: _Buffer) -> list[TrackResult[str]]:
        """Return the configured result values."""
        del input_buffer
        overlap = FeatureOverlap(
            feature=Feature(start=1, end=1, strand=Strand.FORWARD),
            overlapping_feature=Feature(start=1, end=1, strand=Strand.FORWARD),
        )
        return [
            TrackResult(track_identifier=self.identifier, overlap=overlap, annotations=(result,))
            for result in self._results
        ]


class _BlockingTrack(_StaticTrack):
    """Wait until another track starts to prove concurrent execution."""

    def __init__(self, identifier: str, barrier: Barrier) -> None:
        super().__init__(identifier, [identifier])
        self._barrier = barrier

    def annotate(self, input_buffer: _Buffer) -> list[TrackResult[str]]:
        """Synchronise both tracks before producing the configured result."""
        self._barrier.wait(timeout=1)
        return super().annotate(input_buffer)


class _FailingTrack(_StaticTrack):
    """Raise an error from a track worker."""

    def annotate(self, input_buffer: _Buffer) -> list[TrackResult[str]]:
        """Raise the configured worker error."""
        del input_buffer
        raise RuntimeError("track failed")


def _overlap_finder() -> FeatureOverlapFinder[Feature, FeatureOverlap]:
    """Return a finder with one matching feature."""
    return FeatureOverlapFinder(
        candidates=lambda input_buffer: [Feature(start=1, end=10, strand=Strand.FORWARD)],
        matcher=_AlwaysMatcher(),
        overlap_factory=lambda feature, allele: FeatureOverlap(feature=feature, overlapping_feature=allele),
    )


class TestOverlapAnnotationTrack:
    """One finder and annotator composed into a track."""

    def test_preserves_overlap_order_and_annotation_values(self) -> None:
        """Pair batch annotator results with the overlap that produced them."""
        track = OverlapAnnotationTrack(
            identifier="test",
            overlap_finder=_overlap_finder(),
            annotator=_LabelsAnnotator(),
        )

        results = track.annotate(_Buffer([VariationAllele(start=5, end=5, strand=Strand.FORWARD)]))

        assert [result.annotations for result in results] == [("1",)]
        assert results[0].track_identifier == "test"

    def test_rejects_a_misaligned_batch_result(self) -> None:
        """Reject annotators that lose the correspondence to input overlaps."""

        class _BrokenAnnotator(_LabelsAnnotator):
            name = "test.broken"

            def annotate_feature_overlaps(self, overlaps: Sequence[FeatureOverlap]) -> list[list[str]]:
                del overlaps
                return []

        track = OverlapAnnotationTrack(
            identifier="test",
            overlap_finder=_overlap_finder(),
            annotator=_BrokenAnnotator(),
        )

        with pytest.raises(ValueError, match="returned 0 result groups for 1 overlaps"):
            track.annotate(_Buffer([VariationAllele(start=5, end=5, strand=Strand.FORWARD)]))

    def test_rejects_a_nonpositive_batch_size(self) -> None:
        """Require tracks to retain a bounded positive overlap batch."""
        with pytest.raises(ValueError, match="batch_size must be positive"):
            OverlapAnnotationTrack(
                identifier="test",
                overlap_finder=_overlap_finder(),
                annotator=_LabelsAnnotator(),
                batch_size=0,
            )


class TestWorkflowPlan:
    """Ordered runtime tracks and registered workflow builders."""

    def test_runs_tracks_in_declared_order(self) -> None:
        """Concatenate track results in the workflow order."""
        first = OverlapAnnotationTrack(
            identifier="first", overlap_finder=_overlap_finder(), annotator=_LabelsAnnotator()
        )
        second = OverlapAnnotationTrack(
            identifier="second", overlap_finder=_overlap_finder(), annotator=_LabelsAnnotator()
        )
        workflow = WorkflowPlan(name="default", tracks=(first, second))

        results = workflow.annotate(_Buffer([VariationAllele(start=5, end=5, strand=Strand.FORWARD)]))

        assert [result.track_identifier for result in results] == ["first", "second"]

    def test_rejects_duplicate_track_identifiers(self) -> None:
        """Require track result origins to remain unambiguous."""
        track = OverlapAnnotationTrack(
            identifier="duplicate", overlap_finder=_overlap_finder(), annotator=_LabelsAnnotator()
        )

        with pytest.raises(ValueError, match="unique identifiers"):
            WorkflowPlan(name="default", tracks=(track, track))

    def test_threaded_results_match_serial_order(self) -> None:
        """Retain track and result order when independent tracks use threads."""
        workflow = WorkflowPlan(
            name="default",
            tracks=(
                _StaticTrack("first", ["one", "two"]),
                _StaticTrack("second", ["three"]),
            ),
        )
        input_buffer = _Buffer([])

        serial = workflow.annotate(input_buffer)
        threaded = workflow.annotate(input_buffer, threads=2)

        assert [(result.track_identifier, result.annotations) for result in threaded] == [
            (result.track_identifier, result.annotations) for result in serial
        ]
        assert [result.track_identifier for result in threaded] == ["first", "first", "second"]

    def test_runs_independent_tracks_concurrently(self) -> None:
        """Allow independently owned tracks to overlap their blocking work."""
        barrier = Barrier(2)
        workflow = WorkflowPlan(
            name="default",
            tracks=(
                _BlockingTrack("first", barrier),
                _BlockingTrack("second", barrier),
            ),
        )

        results = workflow.annotate(_Buffer([]), threads=2)

        assert [result.track_identifier for result in results] == ["first", "second"]

    def test_propagates_track_worker_errors(self) -> None:
        """Raise a track error in the caller after worker completion."""
        workflow = WorkflowPlan(
            name="default",
            tracks=(
                _FailingTrack("failed", []),
                _StaticTrack("other", ["value"]),
            ),
        )

        with pytest.raises(RuntimeError, match="track failed"):
            workflow.annotate(_Buffer([]), threads=2)

    def test_caps_threads_to_the_number_of_tracks(self) -> None:
        """Use the available tracks when more worker threads are requested."""
        workflow = WorkflowPlan(
            name="default",
            tracks=(
                _StaticTrack("first", ["one"]),
                _StaticTrack("second", ["two"]),
            ),
        )

        results = workflow.annotate(_Buffer([]), threads=10)

        assert [result.annotations for result in results] == [("one",), ("two",)]

    def test_rejects_a_nonpositive_thread_count(self) -> None:
        """Require at least one workflow thread."""
        workflow = WorkflowPlan(name="default", tracks=())

        with pytest.raises(ValueError, match="thread count must be positive"):
            workflow.annotate(_Buffer([]), threads=0)

    def test_default_workflow_builds_registered_tracks(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Build registered tracks through the default workflow factory."""

        @register_track("test.feature")
        def build_feature_track(specification: TrackSpec, context: dict[str, object]) -> AnnotationTrack[Any]:
            del context
            return OverlapAnnotationTrack(
                identifier=specification.identifier,
                overlap_finder=_overlap_finder(),
                annotator=_LabelsAnnotator(),
            )

        monkeypatch.setitem(
            getattr(track_registry_module, "_TRACK_FACTORIES"),
            "test.feature",
            build_feature_track,
        )
        specification = WorkflowSpec(
            name="default",
            tracks=(TrackSpec(identifier="feature", kind="test.feature"),),
        )

        workflow = build_workflow(specification, {})

        assert workflow.name == "default"
        assert workflow.tracks[0].identifier == "feature"

    def test_haplo_workflow_builds_registered_tracks(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Use the same track-factory contract in haplotype mode."""

        def build_feature_track(specification: TrackSpec, context: dict[str, object]) -> AnnotationTrack[Any]:
            del context
            return OverlapAnnotationTrack(
                identifier=specification.identifier,
                overlap_finder=_overlap_finder(),
                annotator=_LabelsAnnotator(),
            )

        monkeypatch.setitem(
            getattr(track_registry_module, "_TRACK_FACTORIES"),
            "test.haplo",
            build_feature_track,
        )
        specification = WorkflowSpec(
            name="haplo",
            tracks=(TrackSpec(identifier="feature", kind="test.haplo"),),
        )

        workflow = build_workflow(specification, {})

        assert workflow.name == "haplo"
        assert workflow.tracks[0].identifier == "feature"


class TestWorkflowRegistry:
    """Track and workflow factory lookup."""

    def test_exposes_builtin_workflows(self) -> None:
        """Register default and haplotype builders during package import."""
        assert available_workflows() == ("default", "haplo")

    def test_unknown_factories_describe_available_names(self) -> None:
        """Report registered names when a factory cannot be resolved."""
        specification = TrackSpec(identifier="unknown", kind="missing")

        with pytest.raises(KeyError, match="Unknown annotation track 'missing'"):
            build_track(specification, {})
        with pytest.raises(KeyError, match="Unknown workflow 'missing'"):
            build_workflow(WorkflowSpec(name="missing", tracks=()), {})

    def test_lists_registered_track_factories_in_order(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Sort track factory names for stable configuration diagnostics."""
        factories = getattr(track_registry_module, "_TRACK_FACTORIES")
        monkeypatch.setitem(factories, "test.zeta", lambda specification, context: None)
        monkeypatch.setitem(factories, "test.alpha", lambda specification, context: None)

        names = available_track_factories()

        assert names.index("test.alpha") < names.index("test.zeta")

    def test_builds_the_builtin_transcript_consequence_track(self) -> None:
        """Inject worker-owned lookup, matcher, and annotator dependencies."""
        transcript = Transcript(
            stable_id="ENST00000000001",
            start=1,
            end=10,
            strand=Strand.FORWARD,
        )
        specification = TrackSpec(
            identifier="transcript",
            kind="transcript.consequence",
            options={"batch_size": 1},
        )
        context = {
            "transcript_candidates": lambda input_buffer: [transcript],
            "transcript_matcher": _AlwaysMatcher(),
            "consequence_annotator": _LabelsAnnotator(),
        }

        track = build_track(specification, context)
        results = track.annotate(_Buffer([VariationAllele(start=5, end=5, strand=Strand.FORWARD)]))

        assert track.batch_size == 1
        assert results[0].annotations == ("1",)

    @pytest.mark.parametrize(
        ("flank", "distance", "expected"),
        [
            pytest.param(6000, 6000, True, id="above-default-boundary"),
            pytest.param(6000, 6001, False, id="above-default-outside"),
            pytest.param(4000, 4000, True, id="below-default-boundary"),
            pytest.param(4000, 4001, False, id="below-default-outside"),
        ],
    )
    def test_transcript_track_matches_the_configured_flank_boundary(
        self,
        flank: int,
        distance: int,
        expected: bool,
    ) -> None:
        """Use the track option instead of the finder's fixed default flank."""
        transcript = Transcript(
            stable_id="ENST00000000001",
            start=10000,
            end=20000,
            strand=Strand.FORWARD,
        )
        specification = TrackSpec(
            identifier="transcript",
            kind="transcript.consequence",
            options={"upstream_distance": flank, "downstream_distance": flank},
        )
        context = {
            "transcript_candidates": lambda input_buffer: [transcript],
            "transcript_matcher": OverlapMatcher(),
            "consequence_annotator": _LabelsAnnotator(),
        }
        variation_allele = VariationAllele(
            start=transcript.end + distance,
            end=transcript.end + distance,
            strand=Strand.FORWARD,
        )

        track = build_track(specification, context)
        results = track.annotate(_Buffer([variation_allele]))

        assert bool(results) is expected

    @pytest.mark.parametrize("name", ["upstream_distance", "downstream_distance"])
    def test_transcript_track_rejects_a_noninteger_distance(self, name: str) -> None:
        """Reject a pairing distance that cannot define an exact coordinate boundary."""
        specification = TrackSpec(
            identifier="transcript",
            kind="transcript.consequence",
            options={name: 5000.0},
        )

        with pytest.raises(TypeError, match=f"{name} must be an integer"):
            build_track(specification, {})

    @pytest.mark.parametrize(
        ("option", "value", "message"),
        [
            pytest.param(
                "readthrough",
                object(),
                "readthrough must be a Readthrough value or name",
                id="readthrough",
            ),
            pytest.param("batch_size", 1.5, "batch_size must be an integer", id="batch-size"),
        ],
    )
    def test_transcript_track_rejects_invalid_option_types(
        self,
        option: str,
        value: object,
        message: str,
    ) -> None:
        """Reject options that cannot construct the transcript track."""
        specification = TrackSpec(
            identifier="transcript",
            kind="transcript.consequence",
            options={option: value},
        )

        with pytest.raises(TypeError, match=message):
            build_track(specification, {})

    def test_transcript_track_requires_every_worker_dependency(self) -> None:
        """Name a missing worker dependency when track construction cannot continue."""
        specification = TrackSpec(
            identifier="transcript",
            kind="transcript.consequence",
        )

        with pytest.raises(KeyError, match="transcript_candidates"):
            build_track(specification, {})

    @pytest.mark.parametrize("name", ["upstream_distance", "downstream_distance"])
    def test_transcript_track_rejects_a_negative_distance(self, name: str) -> None:
        """Reject a pairing distance that would contract the transcript influence interval."""
        specification = TrackSpec(
            identifier="transcript",
            kind="transcript.consequence",
            options={name: -1},
        )

        with pytest.raises(ValueError, match=f"{name} must not be negative"):
            build_track(specification, {})
