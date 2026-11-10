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
    TrackResources,
    TrackResult,
)
from ensembl.vep.workflow.tracks.colocated import build_colocated_tracks
from ensembl.vep.workflow.tracks.consequence import (
    TranscriptCandidateLookup,
    TranscriptConsequenceTrack,
    build_regulatory_tracks,
    build_transcript_overlap_finder,
    build_transcript_tracks,
    transcript_candidate_window,
)
from ensembl.vep.workflow.tracks.custom import build_custom_tracks
from ensembl.vep.workflow.tracks.registry import (
    available_track_builders,
    build_tracks,
    register_track,
    track_builder,
)

__all__ = [
    "AnnotationTrack",
    "OverlapAnnotationTrack",
    "TrackResources",
    "TrackResult",
    "TranscriptCandidateLookup",
    "TranscriptConsequenceTrack",
    "available_track_builders",
    "build_colocated_tracks",
    "build_custom_tracks",
    "build_regulatory_tracks",
    "build_tracks",
    "build_transcript_overlap_finder",
    "build_transcript_tracks",
    "register_track",
    "track_builder",
    "transcript_candidate_window",
]
