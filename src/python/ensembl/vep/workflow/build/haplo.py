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
"""The haplotype annotation workflow builder."""

from __future__ import annotations

from typing import Any

from ensembl.vep.config import Config
from ensembl.vep.workflow.base import Workflow
from ensembl.vep.workflow.build.registry import register_workflow
from ensembl.vep.workflow.tracks.base import AnnotationTrack, TrackResources
from ensembl.vep.workflow.tracks.registry import build_tracks

__all__ = ["HAPLO_TRACKS", "build_haplo_workflow"]

HAPLO_TRACKS = ("consequence.transcript",)
"""The kinds of track a haplotype run annotates with.

Haplotype behaviour belongs to the tracks themselves rather than to runner or execution code, so
this workflow differs from the default one only in which kinds it asks for.
"""


@register_workflow("haplo")
def build_haplo_workflow(config: Config, resources: TrackResources) -> Workflow:
    """Build the haplotype-aware workflow."""
    tracks: list[AnnotationTrack[Any]] = []
    for name in HAPLO_TRACKS:
        tracks.extend(build_tracks(name, config, resources))
    return Workflow(name="haplo", tracks=tuple(tracks))
