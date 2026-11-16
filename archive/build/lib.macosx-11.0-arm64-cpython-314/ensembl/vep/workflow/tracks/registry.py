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
"""Registry for annotation-track factories."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any

from ensembl.vep.workflow.tracks.base import AnnotationTrack, TrackSpec

__all__ = [
    "available_track_factories",
    "build_track",
    "register_track",
    "track_factory",
]

TrackFactory = Callable[[TrackSpec, Mapping[str, object]], AnnotationTrack[Any]]
"""Build one runtime track from its specification and worker context."""

_TRACK_FACTORIES: dict[str, TrackFactory] = {}


def register_track(name: str) -> Callable[[TrackFactory], TrackFactory]:
    """Register a track factory under ``name``."""

    def decorator(factory: TrackFactory) -> TrackFactory:
        _TRACK_FACTORIES[name] = factory
        return factory

    return decorator


def track_factory(name: str) -> TrackFactory:
    """Return the track factory registered under ``name``."""
    try:
        return _TRACK_FACTORIES[name]
    except KeyError as error:
        available = ", ".join(available_track_factories()) or "<none>"
        raise KeyError(f"Unknown annotation track '{name}'. Available tracks: {available}") from error


def build_track(specification: TrackSpec, context: Mapping[str, object]) -> AnnotationTrack[Any]:
    """Build a track from its registered factory."""
    return track_factory(specification.kind)(specification, context)


def available_track_factories() -> tuple[str, ...]:
    """Return registered track factory names in sorted order."""
    return tuple(sorted(_TRACK_FACTORIES))
