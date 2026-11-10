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
"""The default annotation workflow builder."""

from __future__ import annotations

from typing import Any

from ensembl.vep.config import Config
from ensembl.vep.workflow.base import Workflow
from ensembl.vep.workflow.build.registry import register_workflow
from ensembl.vep.workflow.tracks.base import AnnotationTrack, TrackResources
from ensembl.vep.workflow.tracks.registry import build_tracks

__all__ = ["DEFAULT_TRACKS", "build_default_workflow"]

DEFAULT_TRACKS = (
    "consequence.transcript",
    "consequence.regulatory",
    "colocated",
    "custom",
)
"""Every kind of track a default run can annotate with, in output order.

Each kind decides from the configuration how many tracks it contributes, so a kind nothing asks
for contributes none and the order of those that remain is this one. Transcript consequences come
first because that is the order VEP writes them in.
"""


@register_workflow("default")
def build_default_workflow(config: Config, resources: TrackResources) -> Workflow:
    """Build the standard per-variant workflow from the configured tracks.

    Args:
        config (Config): The run's configuration.
        resources (TrackResources): What every track in the workflow shares.

    Returns:
        Workflow: The workflow, holding one live track per configured kind.
    """
    tracks: list[AnnotationTrack[Any]] = []
    for name in DEFAULT_TRACKS:
        tracks.extend(build_tracks(name, config, resources))
    return Workflow(name="default", tracks=tuple(tracks))
