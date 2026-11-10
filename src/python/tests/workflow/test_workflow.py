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
from types import SimpleNamespace
from typing import Any, cast

import pytest

from ensembl.vep.annotators.base import Annotator
from ensembl.vep.common import Strand
from ensembl.vep.caches.interval_tree import Interval
from ensembl.vep.config import Config
from ensembl.vep.models import Feature
from ensembl.vep.models.overlap import FeatureOverlap
from ensembl.vep.models.variation import VariationAllele
from ensembl.vep.overlaps import FeatureOverlapFinder
from ensembl.vep.workflow import (
    AnnotationTrack,
    OverlapAnnotationTrack,
    TrackResources,
    TrackResult,
    Workflow,
    available_track_builders,
    available_workflows,
    build_tracks,
    build_workflow,
    register_track,
)
from ensembl.vep.workflow.build import default as default_build_module
from ensembl.vep.workflow.build import haplo as haplo_build_module
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

    def begin_payload(self, input_buffer: _Buffer) -> None:
        """Prime nothing: these tracks hold no source."""
        del input_buffer

    def attach(self, variations: Sequence[Any], results: Sequence[TrackResult[str]]) -> None:
        """Attach nothing: these tracks are never asked for output."""
        del variations, results

    def close(self) -> None:
        """Release nothing: these tracks hold no source."""

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
        matcher=cast(Any, _AlwaysMatcher()),
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


class TestWorkflow:
    """Ordered runtime tracks and registered workflow builders."""

    def test_runs_tracks_in_declared_order(self) -> None:
        """Concatenate track results in the workflow order."""
        first = OverlapAnnotationTrack(
            identifier="first", overlap_finder=_overlap_finder(), annotator=_LabelsAnnotator()
        )
        second = OverlapAnnotationTrack(
            identifier="second", overlap_finder=_overlap_finder(), annotator=_LabelsAnnotator()
        )
        workflow = Workflow(name="default", tracks=cast(Any, (first, second)))

        results = workflow.annotate(_Buffer([VariationAllele(start=5, end=5, strand=Strand.FORWARD)]))

        assert [result.track_identifier for result in results] == ["first", "second"]

    def test_rejects_duplicate_track_identifiers(self) -> None:
        """Require track result origins to remain unambiguous."""
        track = OverlapAnnotationTrack(
            identifier="duplicate", overlap_finder=_overlap_finder(), annotator=_LabelsAnnotator()
        )

        with pytest.raises(ValueError, match="unique identifiers"):
            Workflow(name="default", tracks=cast(Any, (track, track)))

    def test_threaded_results_match_serial_order(self) -> None:
        """Retain track and result order when independent tracks use threads."""
        workflow = Workflow(
            name="default",
            tracks=cast(
                Any,
                (
                    _StaticTrack("first", ["one", "two"]),
                    _StaticTrack("second", ["three"]),
                ),
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
        workflow = Workflow(
            name="default",
            tracks=cast(
                Any,
                (
                    _BlockingTrack("first", barrier),
                    _BlockingTrack("second", barrier),
                ),
            ),
        )

        results = workflow.annotate(_Buffer([]), threads=2)

        assert [result.track_identifier for result in results] == ["first", "second"]

    def test_propagates_track_worker_errors(self) -> None:
        """Raise a track error in the caller after worker completion."""
        workflow = Workflow(
            name="default",
            tracks=cast(
                Any,
                (
                    _FailingTrack("failed", []),
                    _StaticTrack("other", ["value"]),
                ),
            ),
        )

        with pytest.raises(RuntimeError, match="track failed"):
            workflow.annotate(_Buffer([]), threads=2)

    def test_caps_threads_to_the_number_of_tracks(self) -> None:
        """Use the available tracks when more worker threads are requested."""
        workflow = Workflow(
            name="default",
            tracks=cast(
                Any,
                (
                    _StaticTrack("first", ["one"]),
                    _StaticTrack("second", ["two"]),
                ),
            ),
        )

        results = workflow.annotate(_Buffer([]), threads=10)

        assert [result.annotations for result in results] == [("one",), ("two",)]

    def test_rejects_a_nonpositive_thread_count(self) -> None:
        """Require at least one workflow thread."""
        workflow = Workflow(name="default", tracks=())

        with pytest.raises(ValueError, match="thread count must be positive"):
            workflow.annotate(_Buffer([]), threads=0)


class TestWorkflowBuilding:
    """Building a workflow from the registered track builders."""

    @staticmethod
    def _resources() -> TrackResources:
        """Return resources no stub track reads."""
        return TrackResources(
            source_config=Config().source,
            sequence_source=None,
            synonyms=cast(Any, SimpleNamespace()),
        )

    @staticmethod
    def _register_feature_track(monkeypatch: pytest.MonkeyPatch, name: str) -> None:
        """Register one stub track builder under ``name``."""

        def build_feature_tracks(
            config: Config,
            resources: TrackResources,
        ) -> tuple[AnnotationTrack[Any], ...]:
            del config, resources
            return (
                OverlapAnnotationTrack(
                    identifier="feature",
                    overlap_finder=_overlap_finder(),
                    annotator=_LabelsAnnotator(),
                ),
            )

        monkeypatch.setitem(
            getattr(track_registry_module, "_TRACK_BUILDERS"),
            name,
            build_feature_tracks,
        )

    def test_default_workflow_builds_registered_tracks(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Build registered tracks through the default workflow builder."""
        self._register_feature_track(monkeypatch, "test.feature")
        monkeypatch.setattr(default_build_module, "DEFAULT_TRACKS", ("test.feature",))

        workflow = build_workflow(Config(), self._resources())

        assert workflow.name == "default"
        assert workflow.tracks[0].identifier == "feature"

    def test_haplo_workflow_builds_registered_tracks(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Use the same track-builder contract in haplotype mode."""
        self._register_feature_track(monkeypatch, "test.haplo")
        monkeypatch.setattr(haplo_build_module, "HAPLO_TRACKS", ("test.haplo",))

        workflow = build_workflow(Config(), self._resources(), name="haplo")

        assert workflow.name == "haplo"
        assert workflow.tracks[0].identifier == "feature"

    def test_a_builder_contributing_nothing_leaves_the_workflow_empty(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A kind the configuration does not ask for adds no track."""
        monkeypatch.setitem(
            getattr(track_registry_module, "_TRACK_BUILDERS"),
            "test.absent",
            lambda _config, _resources: (),
        )
        monkeypatch.setattr(default_build_module, "DEFAULT_TRACKS", ("test.absent",))

        workflow = build_workflow(Config(), self._resources())

        assert workflow.tracks == ()


class TestRegistries:
    """Track and workflow builder lookup."""

    def test_exposes_builtin_workflows(self) -> None:
        """Register default and haplotype builders during package import."""
        assert available_workflows() == ("default", "haplo")

    def test_exposes_builtin_tracks(self) -> None:
        """Register every built-in kind of track during package import."""
        assert available_track_builders() == (
            "colocated",
            "consequence.regulatory",
            "consequence.transcript",
            "custom",
        )

    def test_unknown_builders_describe_available_names(self) -> None:
        """Report registered names when a builder cannot be resolved."""
        resources = TestWorkflowBuilding._resources()  # pylint: disable=protected-access

        with pytest.raises(KeyError, match="Unknown annotation track 'missing'"):
            build_tracks("missing", Config(), resources)
        with pytest.raises(KeyError, match="Unknown workflow 'missing'"):
            build_workflow(Config(), resources, name="missing")

    def test_lists_registered_track_builders_in_order(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Sort track builder names for stable configuration diagnostics."""
        builders = getattr(track_registry_module, "_TRACK_BUILDERS")
        monkeypatch.setitem(builders, "test.zeta", lambda _config, _resources: ())
        monkeypatch.setitem(builders, "test.alpha", lambda _config, _resources: ())

        names = available_track_builders()

        assert names.index("test.alpha") < names.index("test.zeta")

    def test_registering_returns_the_builder_unchanged(self) -> None:
        """The decorator hands back what it registered, so a module can also call it."""

        def build_nothing(_config: Config, _resources: TrackResources) -> tuple[AnnotationTrack[Any], ...]:
            return ()

        try:
            assert register_track("test.identity")(build_nothing) is build_nothing
        finally:
            getattr(track_registry_module, "_TRACK_BUILDERS").pop("test.identity", None)
