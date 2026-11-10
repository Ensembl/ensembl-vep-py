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
"""Registry for annotation-track builders."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

from ensembl.vep.config import Config
from ensembl.vep.workflow.tracks.base import AnnotationTrack, TrackResources

__all__ = [
    "available_track_builders",
    "build_tracks",
    "register_track",
    "track_builder",
]

TrackBuilder = Callable[[Config, TrackResources], tuple[AnnotationTrack[Any], ...]]
"""Build every track of one kind that a configuration asks for.

A builder returns as many tracks as the configuration names: none when the kind was not asked
for, one for a kind configured at most once, and one per occurrence for a repeatable option.
"""

_TRACK_BUILDERS: dict[str, TrackBuilder] = {}


def register_track(name: str) -> Callable[[TrackBuilder], TrackBuilder]:
    """Register a track builder under ``name``."""

    def decorator(builder: TrackBuilder) -> TrackBuilder:
        _TRACK_BUILDERS[name] = builder
        return builder

    return decorator


def track_builder(name: str) -> TrackBuilder:
    """Return the track builder registered under ``name``.

    Raises:
        KeyError: If no builder is registered under that name.
    """
    try:
        return _TRACK_BUILDERS[name]
    except KeyError as error:
        available = ", ".join(available_track_builders()) or "<none>"
        raise KeyError(f"Unknown annotation track '{name}'. Available tracks: {available}") from error


def build_tracks(
    name: str,
    config: Config,
    resources: TrackResources,
) -> tuple[AnnotationTrack[Any], ...]:
    """Build every track of one registered kind."""
    return track_builder(name)(config, resources)


def available_track_builders() -> tuple[str, ...]:
    """Return registered track builder names in sorted order."""
    return tuple(sorted(_TRACK_BUILDERS))
