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
"""Workflow plan models and execution."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import Any

from attrs import define

from ensembl.vep.overlaps.base import VariationAlleleBuffer
from ensembl.vep.workflow.tracks.base import AnnotationTrack, TrackResult, TrackSpec

__all__ = ["WorkflowPlan", "WorkflowSpec"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class WorkflowSpec:
    """Describe one worker-instantiable annotation workflow.

    Attributes:
        name (str): Registered workflow implementation name.
        tracks (tuple[TrackSpec, ...]): Ordered configured tracks.
    """

    name: str
    tracks: tuple[TrackSpec, ...]


@define(kw_only=True, auto_attribs=True)
class WorkflowPlan:
    """Run ordered annotation tracks for an input buffer.

    Attributes:
        name (str): Workflow implementation name.
        tracks (tuple[AnnotationTrack[Any], ...]): Ordered runtime tracks.
    """

    name: str
    tracks: tuple[AnnotationTrack[Any], ...]

    def __attrs_post_init__(self) -> None:
        identifiers = [track.identifier for track in self.tracks]
        if len(set(identifiers)) != len(identifiers):
            raise ValueError("Workflow tracks must have unique identifiers.")

    def annotate(
        self,
        input_buffer: VariationAlleleBuffer,
        *,
        threads: int = 1,
    ) -> list[TrackResult[Any]]:
        """Return ordered results from every track.

        Args:
            input_buffer (VariationAlleleBuffer): Read-only input for every track.
            threads (int): Maximum independently runnable tracks.

        Returns:
            list[TrackResult[Any]]: Results ordered by declared track order.

        Raises:
            ValueError: If ``threads`` is less than one.
        """
        if threads < 1:
            raise ValueError("Workflow thread count must be positive.")
        if threads == 1 or len(self.tracks) < 2:
            return self._annotate_serial(input_buffer)

        workers = min(threads, len(self.tracks))
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(track.annotate, input_buffer) for track in self.tracks]
            results: list[TrackResult[Any]] = []
            for future in futures:
                results.extend(future.result())
            return results

    def _annotate_serial(self, input_buffer: VariationAlleleBuffer) -> list[TrackResult[Any]]:
        """Return ordered results without creating worker threads."""
        results: list[TrackResult[Any]] = []
        for track in self.tracks:
            results.extend(track.annotate(input_buffer))
        return results
