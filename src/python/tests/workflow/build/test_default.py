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
"""Tests for building a workflow from a configuration."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import pytest

from ensembl.vep.common import Readthrough
from ensembl.vep.config import Config
from ensembl.vep.config.sections.common.annotation import AnnotationConfig
from ensembl.vep.config.sections.common.source import SourceConfig
from ensembl.vep.overlaps import TranscriptOverlapFinder
from ensembl.vep.workflow import build_workflow, workflow_for_command
from ensembl.vep.workflow.build.default import DEFAULT_TRACKS
from ensembl.vep.workflow.tracks.base import TrackResources
from ensembl.vep.workflow.tracks.consequence import transcript as transcript_module
from ensembl.vep.workflow.tracks.consequence.transcript import TranscriptConsequenceTrack


@pytest.fixture(name="stub_transcript_source")
def fixture_stub_transcript_source(monkeypatch: pytest.MonkeyPatch) -> None:
    """Build the transcript track over a stub source instead of a real file."""
    monkeypatch.setattr(
        transcript_module,
        "build_transcript_source",
        lambda _config, **_shared: SimpleNamespace(transcript_pool=None, transcript_pool_poisoned=False),
    )


def _resources(config: Config) -> TrackResources:
    """Return track resources over no reference sequence."""
    return TrackResources(
        source_config=config.source,
        sequence_source=None,
        synonyms=cast(Any, SimpleNamespace()),
    )


class TestDefaultWorkflow:
    """What a default run's workflow holds."""

    @pytest.mark.parametrize(
        ("flank", "batch_size", "readthrough"),
        [
            pytest.param(4000, 17, "none", id="below-default-flank"),
            pytest.param(6000, 31, "utr", id="above-default-flank"),
        ],
    )
    def test_configuration_reaches_the_runtime_track(
        self,
        flank: int,
        batch_size: int,
        readthrough: str,
        stub_transcript_source: None,
    ) -> None:
        """Carry source and annotation options into the track the builder makes."""
        del stub_transcript_source
        config = Config(
            command="run.default",
            annotation=AnnotationConfig(batch_size=batch_size, readthrough=readthrough),
            source=SourceConfig(gff=Path("transcripts.gff3"), distance=flank),
        )

        workflow = build_workflow(config, _resources(config))

        (built,) = workflow.tracks
        track = cast(TranscriptConsequenceTrack, built)
        finder = cast(TranscriptOverlapFinder, track.overlap_finder)
        assert workflow.name == "default"
        assert track.identifier == "transcript"
        assert track.batch_size == batch_size
        assert finder.upstream_distance == flank
        assert finder.downstream_distance == flank
        assert finder.readthrough is Readthrough[readthrough.upper()]

    def test_empty_command_builds_the_python_api_default(self, stub_transcript_source: None) -> None:
        """Treat a directly constructed Config as the default run mode."""
        del stub_transcript_source
        config = Config(source=SourceConfig(gff=Path("transcripts.gff3")))

        workflow = build_workflow(config, _resources(config))

        assert workflow.name == "default"
        assert [track.identifier for track in workflow.tracks] == ["transcript"]

    def test_only_configured_kinds_contribute_tracks(self, stub_transcript_source: None) -> None:
        """the kinds nothing asks for contribute no tracks at all."""
        del stub_transcript_source
        config = Config(source=SourceConfig(gff=Path("transcripts.gff3")))

        workflow = build_workflow(config, _resources(config))

        assert len(DEFAULT_TRACKS) == 4
        assert len(workflow.tracks) == 1

    def test_custom_annotation_is_refused_rather_than_ignored(self, stub_transcript_source: None) -> None:
        """a run asking for an unimplemented kind is told, not quietly annotated without it."""
        del stub_transcript_source
        config = Config(source=SourceConfig(gff=Path("transcripts.gff3"), custom=("features.bed",)))

        with pytest.raises(NotImplementedError, match="custom annotation sources"):
            build_workflow(config, _resources(config))


class TestCommandSelection:
    """Which workflow a command annotates with."""

    @pytest.mark.parametrize("command", ["", "run.default"])
    def test_default_commands_select_the_default_workflow(self, command: str) -> None:
        """both spellings of a default run resolve to the same workflow."""
        assert workflow_for_command(command) == "default"

    @pytest.mark.parametrize("command", ["run.haplo", "catalogue.show"])
    def test_unsupported_commands_fail_clearly(self, command: str) -> None:
        """Reject modes without a complete source-to-output workflow."""
        with pytest.raises(ValueError, match="Only 'run.default' is implemented"):
            workflow_for_command(command)

    def test_an_unknown_workflow_name_lists_what_exists(self) -> None:
        """asking for a workflow nobody registered names the ones that are."""
        config = Config()

        with pytest.raises(KeyError, match="Available workflows: default, haplo"):
            build_workflow(config, _resources(config), name="nonexistent")
