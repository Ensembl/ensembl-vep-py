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
"""Compile configuration into annotation workflows."""

from __future__ import annotations

from collections.abc import Mapping

from ensembl.vep.config import Config
from ensembl.vep.workflow.plans import WorkflowPlan, WorkflowSpec, build_workflow
from ensembl.vep.workflow.tracks import TrackSpec

__all__ = ["compile_workflow", "compile_workflow_spec"]

_DEFAULT_COMMANDS = frozenset(("", "run.default"))


def compile_workflow_spec(config: Config) -> WorkflowSpec:
    """Return the workflow specification supported by ``config``.

    An empty command represents direct Python construction of the default runner. Haplotype and
    other future run modes remain unavailable until their complete source-to-output path exists.

    Raises:
        ValueError: If the selected command has no implemented annotation workflow.
    """
    if config.command not in _DEFAULT_COMMANDS:
        raise ValueError(
            f"Cannot compile an annotation workflow for command '{config.command}'. "
            "Only 'run.default' is implemented."
        )

    annotation = config.annotation
    return WorkflowSpec(
        name="default",
        tracks=(
            TrackSpec(
                identifier="transcript",
                kind="transcript.consequence",
                options={
                    "batch_size": annotation.batch_size,
                    "upstream_distance": config.source.upstream_distance,
                    "downstream_distance": config.source.downstream_distance,
                    "readthrough": annotation.readthrough,
                },
            ),
        ),
    )


def compile_workflow(config: Config, worker_context: Mapping[str, object]) -> WorkflowPlan:
    """Build the configured workflow through the registered factories."""
    return build_workflow(compile_workflow_spec(config), worker_context)
