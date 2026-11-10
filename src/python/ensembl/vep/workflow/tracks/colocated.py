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
"""Co-located variant tracks.

A placeholder: ``ColocatedVariantOverlapFinder`` exists, but there is no known-variant source to
read and no option asking for one, so no run configures this track yet.
"""

from __future__ import annotations

from typing import Any

from ensembl.vep.config import Config
from ensembl.vep.workflow.tracks.base import AnnotationTrack, TrackResources
from ensembl.vep.workflow.tracks.registry import register_track

__all__ = ["build_colocated_tracks"]


@register_track("colocated")
def build_colocated_tracks(
    config: Config,
    resources: TrackResources,
) -> tuple[AnnotationTrack[Any], ...]:
    """Build no tracks: nothing configures co-located variants yet."""
    del config, resources
    return ()
