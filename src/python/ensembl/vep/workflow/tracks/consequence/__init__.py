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
"""Consequence annotation tracks for overlapping features."""

from __future__ import annotations

from ensembl.vep.workflow.tracks.consequence.regulatory import build_regulatory_tracks
from ensembl.vep.workflow.tracks.consequence.transcript import (
    TranscriptCandidateLookup,
    TranscriptConsequenceTrack,
    build_transcript_overlap_finder,
    build_transcript_tracks,
    transcript_candidate_window,
)

__all__ = [
    "TranscriptCandidateLookup",
    "TranscriptConsequenceTrack",
    "build_regulatory_tracks",
    "build_transcript_overlap_finder",
    "build_transcript_tracks",
    "transcript_candidate_window",
]
