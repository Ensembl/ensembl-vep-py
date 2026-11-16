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
"""Registry for workflow-plan factories."""

from __future__ import annotations

from collections.abc import Callable, Mapping

from ensembl.vep.workflow.plans.base import WorkflowPlan, WorkflowSpec

__all__ = [
    "available_workflows",
    "build_workflow",
    "register_workflow",
    "workflow_factory",
]

WorkflowFactory = Callable[[WorkflowSpec, Mapping[str, object]], WorkflowPlan]
"""Build one runtime workflow from its specification and worker context."""

_WORKFLOW_FACTORIES: dict[str, WorkflowFactory] = {}


def register_workflow(name: str) -> Callable[[WorkflowFactory], WorkflowFactory]:
    """Register a workflow factory under ``name``."""

    def decorator(factory: WorkflowFactory) -> WorkflowFactory:
        _WORKFLOW_FACTORIES[name] = factory
        return factory

    return decorator


def workflow_factory(name: str) -> WorkflowFactory:
    """Return the workflow factory registered under ``name``."""
    try:
        return _WORKFLOW_FACTORIES[name]
    except KeyError as error:
        available = ", ".join(available_workflows()) or "<none>"
        raise KeyError(f"Unknown workflow '{name}'. Available workflows: {available}") from error


def build_workflow(specification: WorkflowSpec, context: Mapping[str, object]) -> WorkflowPlan:
    """Build a workflow from its registered factory."""
    return workflow_factory(specification.name)(specification, context)


def available_workflows() -> tuple[str, ...]:
    """Return registered workflow names in sorted order."""
    return tuple(sorted(_WORKFLOW_FACTORIES))
