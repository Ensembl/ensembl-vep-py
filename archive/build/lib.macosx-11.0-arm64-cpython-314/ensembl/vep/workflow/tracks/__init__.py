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
"""Annotation tracks and their registry."""

from __future__ import annotations

from ensembl.vep.workflow.tracks.base import (
    AnnotationTrack,
    OverlapAnnotationTrack,
    TrackResult,
    TrackSpec,
)
from ensembl.vep.workflow.tracks.registry import (
    available_track_factories,
    build_track,
    register_track,
    track_factory,
)
from ensembl.vep.workflow.tracks.transcript_consequence import build_transcript_consequence_track

__all__ = [
    "AnnotationTrack",
    "OverlapAnnotationTrack",
    "TrackResult",
    "TrackSpec",
    "available_track_factories",
    "build_track",
    "build_transcript_consequence_track",
    "register_track",
    "track_factory",
]
