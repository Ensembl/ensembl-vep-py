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
"""Input buffers: a batch of input objects, indexed for overlap discovery."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from typing import Generic, TypeVar

from attrs import define, field

from ensembl.vep.caches.interval_tree import Interval, IntervalTree
from ensembl.vep.matcher.base import IntervalLike

__all__ = ["InputBuffer"]

ObjectT = TypeVar("ObjectT", bound=IntervalLike)
"""What object type an input buffer holds; anything with inclusive coordinates."""


@define(kw_only=True, auto_attribs=True)
class InputBuffer(Generic[ObjectT]):
    """Indexed batch of input objects for overlap discovery.

    The index is built on the first query and discarded whenever the batch changes, so a buffer
    that is only ever filled costs nothing to index.

    Attributes:
        objects (tuple[ObjectT, ...]): The objects in the current batch, in input order.
    """

    objects: tuple[ObjectT, ...] = field(default=())
    _interval_tree: IntervalTree[ObjectT] | None = field(default=None, init=False, repr=False)
    _exhausted: bool = field(default=False, init=False)

    @property
    def is_empty(self) -> bool:
        """Return whether the input this batch was drawn from is spent."""

        return self._exhausted

    def next(self, items: Iterator[ObjectT]) -> None:
        """Reset the batch and refill it from ``items``.

        Args:
            items (Iterator[ObjectT]): The objects left to buffer. The iterator stays with the
                caller, so a buffer holds nothing that cannot cross a process boundary.

        Raises:
            NotImplementedError: Always; a refill policy defines this.
        """

        raise NotImplementedError

    def get_objects_by_intervals(self, intervals: Iterable[Interval]) -> list[ObjectT]:
        """Return the objects meeting any candidate interval, in input order.

        Args:
            intervals (Iterable[Interval]): The candidate intervals to query.

        Returns:
            list[ObjectT]: Every object meeting at least one interval, reported once.
        """

        if not self.objects:
            return []

        tree = self._tree()
        positions: set[int] = set()
        for interval in intervals:
            positions.update(tree.get_overlapping_positions(interval.lookup_start, interval.lookup_end))
        return [self.objects[position] for position in sorted(positions)]

    def get_overlapping_objects(self, start: int, end: int) -> list[ObjectT]:
        """Return the objects overlapping one inclusive interval, in input order.

        Args:
            start (int): The 1-based inclusive start coordinate.
            end (int): The inclusive end coordinate.

        Returns:
            list[ObjectT]: The overlapping objects, empty when the batch is.
        """

        if not self.objects:
            return []
        return self._tree().get_overlapping_features(start, end)

    def clear(self) -> None:
        """Drop the index built for the batch, leaving the batch itself alone."""

        self._interval_tree = None

    def _fill(self, objects: Sequence[ObjectT], *, exhausted: bool = False) -> None:
        """Replace the batch, discarding the index built for the previous one.

        Args:
            objects (Sequence[ObjectT]): The new batch, in input order.
            exhausted (bool): Whether the input the batch was drawn from is now spent.
        """

        self.objects = tuple(objects)
        self._interval_tree = None
        self._exhausted = exhausted

    def _tree(self) -> IntervalTree[ObjectT]:
        """Return the index for the current batch, building it on first query."""

        if self._interval_tree is None:
            self._interval_tree = IntervalTree(features=list(self.objects))
        return self._interval_tree
