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
"""Tests for workflow compilation."""

from __future__ import annotations

from types import SimpleNamespace
from typing import cast

import pytest

from ensembl.vep.common import Readthrough
from ensembl.vep.config import Config
from ensembl.vep.config.sections.common.annotation import AnnotationConfig
from ensembl.vep.config.sections.common.source import SourceConfig
from ensembl.vep.matcher import OverlapMatcher
from ensembl.vep.overlaps import TranscriptOverlapFinder
from ensembl.vep.workflow import compile_workflow, compile_workflow_spec


def _worker_context() -> dict[str, object]:
    """Return the worker-owned dependencies required by the transcript track."""
    return {
        "transcript_candidates": lambda _input_buffer: (),
        "transcript_matcher": OverlapMatcher(),
        "consequence_annotator": SimpleNamespace(name="test.consequence"),
    }


class TestWorkflowCompilation:
    """Configuration-to-plan compilation."""

    @pytest.mark.parametrize(
        ("flank", "batch_size", "readthrough"),
        [
            pytest.param(4000, 17, "none", id="below-default-flank"),
            pytest.param(6000, 31, "utr", id="above-default-flank"),
        ],
    )
    def test_default_config_controls_the_specification_and_runtime_track(
        self,
        flank: int,
        batch_size: int,
        readthrough: str,
    ) -> None:
        """Carry source and annotation options into the registered track factory."""
        config = Config(
            command="run.default",
            annotation=AnnotationConfig(batch_size=batch_size, readthrough=readthrough),
            source=SourceConfig(distance=flank),
        )

        specification = compile_workflow_spec(config)
        workflow = compile_workflow(config, _worker_context())

        track_specification = specification.tracks[0]
        assert specification.name == "default"
        assert track_specification.identifier == "transcript"
        assert track_specification.kind == "transcript.consequence"
        assert track_specification.options == {
            "batch_size": batch_size,
            "upstream_distance": flank,
            "downstream_distance": flank,
            "readthrough": readthrough,
        }

        track = workflow.tracks[0]
        finder = cast(TranscriptOverlapFinder, track.overlap_finder)
        assert workflow.name == "default"
        assert track.batch_size == batch_size
        assert finder.upstream_distance == flank
        assert finder.readthrough is Readthrough[readthrough.upper()]

    def test_empty_command_compiles_the_python_api_default(self) -> None:
        """Treat a directly constructed Config as the default run mode."""
        specification = compile_workflow_spec(Config())

        assert specification.name == "default"
        assert specification.tracks[0].kind == "transcript.consequence"

    @pytest.mark.parametrize("command", ["run.haplo", "catalogue.show"])
    def test_unsupported_commands_fail_clearly(self, command: str) -> None:
        """Reject modes without a complete source-to-output workflow."""
        with pytest.raises(ValueError, match="Only 'run.default' is implemented"):
            compile_workflow_spec(Config(command=command))
