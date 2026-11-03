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
"""Custom annotation tracks.

A placeholder: ``CustomOverlapFinder`` and the ``CustomFeature`` model exist, but no reader parses
a ``--custom`` spec, so a run that asks for one is told rather than quietly annotating without it.
``--custom`` is repeatable, so this is the builder that returns one track per occurrence once the
readers exist.
"""

from __future__ import annotations

from typing import Any

from ensembl.vep.config import Config
from ensembl.vep.workflow.tracks.base import AnnotationTrack, TrackResources
from ensembl.vep.workflow.tracks.registry import register_track

__all__ = ["build_custom_tracks"]


@register_track("custom")
def build_custom_tracks(
    config: Config,
    resources: TrackResources,
) -> tuple[AnnotationTrack[Any], ...]:
    """Build no tracks, or refuse a run that asked for custom annotation.

    Raises:
        NotImplementedError: If ``--custom`` was given.
    """
    del resources
    if config.source.custom:
        raise NotImplementedError("custom annotation sources are not implemented yet")
    return ()
