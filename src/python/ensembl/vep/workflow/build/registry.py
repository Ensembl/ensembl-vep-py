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
"""Registry for workflow builders."""

from __future__ import annotations

from collections.abc import Callable

from ensembl.vep.config import Config
from ensembl.vep.workflow.base import Workflow
from ensembl.vep.workflow.tracks.base import TrackResources

__all__ = [
    "available_workflows",
    "build_workflow",
    "register_workflow",
    "workflow_builder",
    "workflow_for_command",
]

WorkflowBuilder = Callable[[Config, TrackResources], Workflow]
"""Build one runtime workflow from a configuration and the resources its tracks share."""

_WORKFLOW_BUILDERS: dict[str, WorkflowBuilder] = {}

_DEFAULT_COMMANDS = frozenset(("", "run.default"))
"""The commands the default workflow serves. An empty command is direct Python construction."""


def register_workflow(name: str) -> Callable[[WorkflowBuilder], WorkflowBuilder]:
    """Register a workflow builder under ``name``."""

    def decorator(builder: WorkflowBuilder) -> WorkflowBuilder:
        _WORKFLOW_BUILDERS[name] = builder
        return builder

    return decorator


def workflow_builder(name: str) -> WorkflowBuilder:
    """Return the workflow builder registered under ``name``.

    Raises:
        KeyError: If no builder is registered under that name.
    """
    try:
        return _WORKFLOW_BUILDERS[name]
    except KeyError as error:
        available = ", ".join(available_workflows()) or "<none>"
        raise KeyError(f"Unknown workflow '{name}'. Available workflows: {available}") from error


def workflow_for_command(command: str) -> str:
    """Return the name of the workflow a command annotates with.

    Cheap enough to call before any worker starts, which is where an unimplemented command should
    be reported rather than inside every worker process.

    Raises:
        ValueError: If the command has no implemented annotation workflow.
    """
    if command not in _DEFAULT_COMMANDS:
        raise ValueError(
            f"Cannot build an annotation workflow for command '{command}'. "
            "Only 'run.default' is implemented."
        )
    return "default"


def build_workflow(config: Config, resources: TrackResources, *, name: str | None = None) -> Workflow:
    """Build the workflow a configuration describes.

    Args:
        config (Config): The run's configuration.
        resources (TrackResources): What every track in the workflow shares.
        name (str | None): A registered workflow name, or None for the configured command's own.

    Returns:
        Workflow: The runtime workflow.
    """
    chosen = name if name is not None else workflow_for_command(config.command)
    return workflow_builder(chosen)(config, resources)


def available_workflows() -> tuple[str, ...]:
    """Return registered workflow names in sorted order."""
    return tuple(sorted(_WORKFLOW_BUILDERS))
