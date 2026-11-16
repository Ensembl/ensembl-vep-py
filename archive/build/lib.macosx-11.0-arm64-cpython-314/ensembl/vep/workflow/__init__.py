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
"""Compose configured annotation tracks into executable workflows."""

from __future__ import annotations

from ensembl.vep.workflow.execution import ExecutionSpec, ordered_map
from ensembl.vep.workflow.plans import (
    WorkflowPlan,
    WorkflowSpec,
    available_workflows,
    build_default_workflow,
    build_haplo_workflow,
    build_workflow,
    register_workflow,
)
from ensembl.vep.workflow.tracks import (
    AnnotationTrack,
    OverlapAnnotationTrack,
    TrackResult,
    TrackSpec,
    available_track_factories,
    build_track,
    build_transcript_consequence_track,
    register_track,
)
from ensembl.vep.workflow.compiler import compile_workflow, compile_workflow_spec

__all__ = [
    "AnnotationTrack",
    "ExecutionSpec",
    "OverlapAnnotationTrack",
    "TrackResult",
    "TrackSpec",
    "WorkflowPlan",
    "WorkflowSpec",
    "available_track_factories",
    "available_workflows",
    "build_default_workflow",
    "build_haplo_workflow",
    "build_track",
    "build_transcript_consequence_track",
    "build_workflow",
    "compile_workflow",
    "compile_workflow_spec",
    "ordered_map",
    "register_track",
    "register_workflow",
]
