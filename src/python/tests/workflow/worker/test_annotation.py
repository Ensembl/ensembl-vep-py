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
"""Tests for the annotation worker."""

from __future__ import annotations

import pickle
from collections.abc import Sequence
from functools import partial
from types import SimpleNamespace
from typing import Any, cast

import pytest

from ensembl.vep.common import CoordSystemLevel, Strand
from ensembl.vep.config import Config
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.models.variation import Variation, VariationAllele
from ensembl.vep.output.factory import AnnotatedVariant
from ensembl.vep.stats import RunSummaryDelta
from ensembl.vep.workflow.base import Workflow
from ensembl.vep.workflow.tracks.base import TrackResult
from ensembl.vep.workflow.worker import AnnotationBufferResult, AnnotationBufferWorker, Worker


class _RecordingTrack:
    """Record every group of results handed to it."""

    def __init__(self, identifier: str) -> None:
        """Initialise the record of attached groups."""
        self.identifier = identifier
        self.attached: list[tuple[int, ...]] = []
        self.payloads: list[object] = []
        self.closed = False

    def begin_payload(self, input_buffer: object) -> None:
        """Record the payload this track was primed with."""
        self.payloads.append(input_buffer)

    def annotate(self, input_buffer: object) -> list[TrackResult[Any]]:
        """Return nothing: these tests are about dispatch, not annotation."""
        del input_buffer
        return []

    def attach(self, variations: Sequence[Variation], results: Sequence[TrackResult[Any]]) -> None:
        """Record how many results this track was given."""
        del variations
        self.attached.append(tuple(len(result.annotations) for result in results))

    def close(self) -> None:
        """Record that the workflow released this track."""
        self.closed = True


def _worker(tracks: tuple[_RecordingTrack, ...]) -> AnnotationBufferWorker:
    """Return a worker over a workflow of recording tracks, with nothing else built."""
    # pylint: disable=protected-access
    worker = object.__new__(AnnotationBufferWorker)
    worker.workflow = Workflow(name="default", tracks=cast(Any, tracks))
    worker._track_identifiers = frozenset(track.identifier for track in tracks)
    return worker


def _result(track_identifier: str, annotations: tuple[str, ...]) -> TrackResult[Any]:
    """Return one track result carrying ``annotations``."""
    return TrackResult(
        track_identifier=track_identifier,
        overlap=cast(Any, SimpleNamespace()),
        annotations=annotations,
    )


class TestProcessSafety:
    """What crosses a worker-process boundary."""

    def test_the_worker_factory_is_process_safe(self) -> None:
        """The factory is what crosses to a worker process, so it is what has to pickle."""
        worker_factory = partial(AnnotationBufferWorker, Config())

        assert pickle.loads(pickle.dumps(worker_factory)).func is AnnotationBufferWorker

    def test_the_worker_runs_payloads_through_the_executor_contract(self) -> None:
        """An executor knows a worker only by ``run``, which is what makes it interchangeable."""
        assert isinstance(object.__new__(AnnotationBufferWorker), Worker)

    def test_worker_results_are_process_safe(self) -> None:
        """Serialise detached annotations without carrying their live model graph."""
        variant = AnnotatedVariant(
            id="rs1",
            seq_region_name="1",
            start=10,
            end=10,
            strand=1,
            allele_string="A/G",
            annotations=({"Consequence": ["missense_variant"]},),
        )
        result = AnnotationBufferResult(
            variants=(variant,),
            summary_delta=RunSummaryDelta(variants=1),
        )

        restored = pickle.loads(pickle.dumps(result))

        assert restored == result


class TestResultDispatch:
    """Which track is handed which results."""

    def test_each_track_is_given_its_own_results(self) -> None:
        """Results are routed by the identifier of the track that produced them."""
        transcript = _RecordingTrack("transcript")
        regulatory = _RecordingTrack("regulatory")
        worker = _worker((transcript, regulatory))

        # pylint: disable=protected-access
        worker._attach_results(
            (),
            [_result("transcript", ("missense",)), _result("regulatory", ("promoter", "tfbs"))],
        )

        assert transcript.attached == [(1,)]
        assert regulatory.attached == [(2,)]

    def test_a_track_with_no_results_is_still_asked_to_attach(self) -> None:
        """A window that found nothing for a track releases what the last one attached."""
        transcript = _RecordingTrack("transcript")
        regulatory = _RecordingTrack("regulatory")
        worker = _worker((transcript, regulatory))

        worker._attach_results((), [_result("transcript", ("missense",))])  # pylint: disable=protected-access

        assert regulatory.attached == [()]

    def test_results_from_an_unknown_track_are_refused(self) -> None:
        """A result no track in the workflow owns is reported rather than dropped."""
        worker = _worker((_RecordingTrack("transcript"),))

        with pytest.raises(ValueError, match="no track regulatory"):
            worker._attach_results((), [_result("regulatory", ())])  # pylint: disable=protected-access


class TestWorkflowLifetime:
    """What a workflow does to its tracks."""

    def test_beginning_a_payload_primes_every_track(self) -> None:
        """Each track is primed once with the whole payload."""
        transcript = _RecordingTrack("transcript")
        regulatory = _RecordingTrack("regulatory")
        workflow = Workflow(name="default", tracks=cast(Any, (transcript, regulatory)))
        payload = SimpleNamespace()

        workflow.begin_payload(cast(Any, payload))

        assert transcript.payloads == [payload]
        assert regulatory.payloads == [payload]

    def test_leaving_the_workflow_closes_every_track(self) -> None:
        """A workflow used as a context manager releases each track's source."""
        transcript = _RecordingTrack("transcript")
        regulatory = _RecordingTrack("regulatory")

        with Workflow(name="default", tracks=cast(Any, (transcript, regulatory))):
            pass

        assert transcript.closed
        assert regulatory.closed

    def test_tracks_must_have_unique_identifiers(self) -> None:
        """Two tracks under one identifier cannot form a workflow."""
        with pytest.raises(ValueError, match="unique identifiers"):
            Workflow(name="default", tracks=cast(Any, (_RecordingTrack("x"), _RecordingTrack("x"))))

    def test_a_named_track_is_found_by_identifier(self) -> None:
        """A workflow resolves its own tracks by identifier."""
        transcript = _RecordingTrack("transcript")
        workflow = Workflow(name="default", tracks=cast(Any, (transcript,)))

        assert workflow.track("transcript") is transcript
        with pytest.raises(KeyError, match="no track 'regulatory'"):
            workflow.track("regulatory")


def _allele(seq_region_name: str, start: int) -> VariationAllele:
    """Return one allele placed on a sequence region."""
    return VariationAllele(
        start=start,
        end=start,
        strand=Strand.FORWARD,
        slice=Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name=seq_region_name,
            start=1,
            end=250_000_000,
            strand=Strand.FORWARD,
        ),
    )


class TestPayloadBuffer:
    """The buffer every track is primed from."""

    def test_the_whole_batch_is_buffered(self) -> None:
        """Hold every allele of the batch, since the prefetch primes from all of them."""
        alleles = [_allele("1", start) for start in (10, 20, 30)]

        buffer = AnnotationBufferWorker._payload_buffer(alleles)  # pylint: disable=protected-access

        assert list(buffer.objects) == alleles
        assert buffer.is_empty

    def test_an_empty_batch_leaves_an_empty_buffer(self) -> None:
        """A batch whose variants are all on unknown regions buffers nothing."""

        buffer = AnnotationBufferWorker._payload_buffer([])  # pylint: disable=protected-access

        assert buffer.objects == ()

    def test_unsorted_input_is_refused(self) -> None:
        """The refill policy is what rejects an unsorted input, so driving it has to reach it."""
        # Each start is far enough behind the one before it to count, and six of them is one
        # more than the policy tolerates.
        alleles = [_allele("1", start) for start in range(7_000, 0, -1_000)]

        with pytest.raises(RuntimeError, match="appears to be unsorted"):
            AnnotationBufferWorker._payload_buffer(alleles)  # pylint: disable=protected-access
