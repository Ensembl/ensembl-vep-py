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
"""Input-buffer models."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Generic, Protocol, TypeVar

from attrs import define, field

from ensembl.vep.caches.interval_index import Interval, IntervalIndex, IntervalTreeIndex

__all__ = ["InputBuffer", "IntervalObject"]


class IntervalObject(Protocol):
    """Describe an object that can be indexed by an interval."""

    start: int
    end: int


ObjectT = TypeVar("ObjectT", bound=IntervalObject)


@define(kw_only=True, auto_attribs=True)
class InputBuffer(Generic[ObjectT]):
    """Group input objects and index them for overlap discovery."""

    objects: tuple[ObjectT, ...] = field(converter=tuple)
    interval_index: IntervalIndex[ObjectT] = field(factory=IntervalTreeIndex, repr=False)
    """Candidate index used for input-object lookup."""

    def __attrs_post_init__(self) -> None:
        """Index input objects after construction."""

        indexed_objects = ((obj, Interval(start=obj.start, end=obj.end)) for obj in self.objects)
        self.interval_index.extend(indexed_objects)

    def get_objects_by_intervals(self, intervals: Iterable[Interval]) -> list[ObjectT]:
        """Return input objects that meet any candidate interval."""

        return self.interval_index.get_objects_by_intervals(intervals)

    def clear(self) -> None:
        """Clear indexed input objects."""

        self.interval_index.clear()
