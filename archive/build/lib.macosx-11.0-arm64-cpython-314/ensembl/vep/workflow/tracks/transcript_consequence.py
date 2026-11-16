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
"""Transcript consequence annotation track."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, cast

from ensembl.vep.annotators.base import Annotator
from ensembl.vep.common import DEFAULT_DISTANCE, Readthrough
from ensembl.vep.models.overlap import OverlapConsequence
from ensembl.vep.overlaps import CandidateLookup, FeatureMatcher, TranscriptOverlapFinder
from ensembl.vep.workflow.tracks.base import OverlapAnnotationTrack, TrackSpec
from ensembl.vep.workflow.tracks.registry import register_track

__all__ = ["build_transcript_consequence_track"]


def _context_value(context: Mapping[str, object], name: str) -> object:
    """Return a required worker-context value."""
    try:
        return context[name]
    except KeyError as error:
        raise KeyError(f"Track requires worker context value '{name}'.") from error


def _distance(specification: TrackSpec, name: str) -> int:
    """Return one pairing distance from a track specification.

    Raises:
        TypeError: If the option is not an integer.
        ValueError: If it is negative.
    """

    distance = specification.options.get(name, DEFAULT_DISTANCE)
    if isinstance(distance, bool) or not isinstance(distance, int):
        raise TypeError(f"Transcript consequence track {name} must be an integer.")
    if distance < 0:
        raise ValueError(f"Transcript consequence track {name} must not be negative.")
    return distance


@register_track("transcript.consequence")
def build_transcript_consequence_track(
    specification: TrackSpec,
    context: Mapping[str, object],
) -> OverlapAnnotationTrack[OverlapConsequence]:
    """Build the built-in transcript consequence track.

    ``context`` must provide ``transcript_candidates``, ``transcript_matcher``, and
    ``consequence_annotator``. The worker supplies the candidate callable so it can use the
    configured source, feature cache, and region index.
    """
    readthrough = specification.options.get("readthrough", Readthrough.NONE)
    if isinstance(readthrough, str):
        readthrough = Readthrough[readthrough.upper()]
    if not isinstance(readthrough, Readthrough):
        raise TypeError("Transcript consequence track readthrough must be a Readthrough value or name.")

    batch_size = specification.options.get("batch_size", 256)
    if not isinstance(batch_size, int):
        raise TypeError("Transcript consequence track batch_size must be an integer.")

    upstream_distance = _distance(specification, "upstream_distance")
    downstream_distance = _distance(specification, "downstream_distance")

    candidates = cast(CandidateLookup[Any], _context_value(context, "transcript_candidates"))
    matcher = cast(FeatureMatcher, _context_value(context, "transcript_matcher"))
    annotator = cast(Annotator[OverlapConsequence], _context_value(context, "consequence_annotator"))
    return OverlapAnnotationTrack(
        identifier=specification.identifier,
        overlap_finder=TranscriptOverlapFinder(
            candidates=candidates,
            matcher=matcher,
            readthrough=readthrough,
            upstream_distance=upstream_distance,
            downstream_distance=downstream_distance,
        ),
        annotator=annotator,
        batch_size=batch_size,
    )
