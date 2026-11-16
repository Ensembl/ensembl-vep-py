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
"""Build the haplotype annotation workflow."""

from __future__ import annotations

from collections.abc import Mapping

from ensembl.vep.workflow.plans.base import WorkflowPlan, WorkflowSpec
from ensembl.vep.workflow.plans.registry import register_workflow
from ensembl.vep.workflow.tracks.registry import build_track

__all__ = ["build_haplo_workflow"]


@register_workflow("haplo")
def build_haplo_workflow(specification: WorkflowSpec, context: Mapping[str, object]) -> WorkflowPlan:
    """Build the haplotype-aware workflow.

    Haplotype behaviour belongs to the track factories selected in ``specification``, rather than
    in runner or execution code.
    """
    tracks = tuple(build_track(track_specification, context) for track_specification in specification.tracks)
    return WorkflowPlan(name=specification.name, tracks=tracks)
