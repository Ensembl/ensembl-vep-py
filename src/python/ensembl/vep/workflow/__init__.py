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
"""Executable workflows composed from annotation tracks."""

from __future__ import annotations

from ensembl.vep.workflow.base import Workflow
from ensembl.vep.workflow.build import (
    available_workflows,
    build_default_workflow,
    build_haplo_workflow,
    build_workflow,
    register_workflow,
    workflow_for_command,
)
from ensembl.vep.workflow.executor import ExecutionSpec, Executor, PayloadIterator
from ensembl.vep.workflow.tracks import (
    AnnotationTrack,
    OverlapAnnotationTrack,
    TrackResources,
    TrackResult,
    available_track_builders,
    build_tracks,
    build_transcript_overlap_finder,
    build_transcript_tracks,
    register_track,
    transcript_candidate_window,
)
from ensembl.vep.workflow.worker import AnnotationBufferResult, AnnotationBufferWorker, Worker

__all__ = [
    "AnnotationBufferResult",
    "AnnotationBufferWorker",
    "AnnotationTrack",
    "ExecutionSpec",
    "Executor",
    "OverlapAnnotationTrack",
    "PayloadIterator",
    "TrackResources",
    "TrackResult",
    "Worker",
    "Workflow",
    "available_track_builders",
    "available_workflows",
    "build_default_workflow",
    "build_haplo_workflow",
    "build_tracks",
    "build_transcript_overlap_finder",
    "build_transcript_tracks",
    "build_workflow",
    "register_track",
    "register_workflow",
    "transcript_candidate_window",
]
