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
"""Fixed-size input buffer: one sequence region, at most ``buffer_size`` objects."""

from __future__ import annotations

from collections import deque
from collections.abc import Iterator
from typing import TypeVar

from attrs import define, field

from ensembl.vep.buffer.input_buffer.base import InputBuffer
from ensembl.vep.buffer.registry import register_buffer
from ensembl.vep.models import Feature

__all__ = [
    "DEFAULT_BUFFER_SIZE",
    "DEFAULT_MAX_NOT_ORDERED_VARIANTS",
    "DEFAULT_MAX_NOT_ORDERED_VARIANTS_DISTANCE",
    "FixedSizeInputBuffer",
]

DEFAULT_BUFFER_SIZE = 5000
"""How many objects one batch holds at most."""

DEFAULT_MAX_NOT_ORDERED_VARIANTS = 5
"""How many out-of-order objects are tolerated before the input is called unsorted."""

DEFAULT_MAX_NOT_ORDERED_VARIANTS_DISTANCE = 50
"""How far back an object may start before it counts as out of order."""

ObjectT = TypeVar("ObjectT", bound=Feature)
"""What object type the buffer holds."""


@register_buffer("variation.fixed")
@define(kw_only=True, auto_attribs=True)
class FixedSizeInputBuffer(InputBuffer[ObjectT]):
    """Fixed-size input batch restricted to one sequence region.

    An object read past the end of a batch, or on the next sequence region, is held back for the
    following refill rather than dropped.

    Attributes:
        buffer_size (int): How many objects one batch holds at most.
        max_not_ordered_variants (int): How many out-of-order objects to tolerate.
        max_not_ordered_variants_distance (int): How far back an object may start before it counts
            as out of order, which leaves room for deletions and minimised alleles.
    """

    buffer_size: int = field(default=DEFAULT_BUFFER_SIZE)
    max_not_ordered_variants: int = field(default=DEFAULT_MAX_NOT_ORDERED_VARIANTS)
    max_not_ordered_variants_distance: int = field(default=DEFAULT_MAX_NOT_ORDERED_VARIANTS_DISTANCE)
    _pre_buffer: deque[ObjectT] = field(factory=deque, init=False, repr=False)
    _not_ordered_count: int = field(default=0, init=False)

    def next(self, items: Iterator[ObjectT]) -> None:
        """Reset the batch and refill it from the held-back objects, then from ``items``.

        Args:
            items (Iterator[ObjectT]): The objects left to buffer, in input order.

        Raises:
            RuntimeError: If the input looks unsorted.
            ValueError: If an object carries no slice to group it by.
        """

        batch: list[ObjectT] = []
        seq_region_name: str | None = None

        while self._pre_buffer and len(batch) < self.buffer_size:
            held = self._pre_buffer.popleft()
            held_seq_region_name = self._seq_region_name(held)
            if seq_region_name is None:
                seq_region_name = held_seq_region_name
            elif held_seq_region_name != seq_region_name:
                self._pre_buffer.appendleft(held)
                self._fill(batch)
                return
            batch.append(held)

        if len(batch) == self.buffer_size:
            self._fill(batch)
            return

        previous_start: int | None = None
        for item in items:
            item_seq_region_name = self._seq_region_name(item)
            if seq_region_name is None:
                seq_region_name = item_seq_region_name
            elif item_seq_region_name != seq_region_name:
                self._pre_buffer.append(item)
                self._fill(batch)
                return

            if len(batch) == self.buffer_size:
                self._pre_buffer.append(item)
                self._fill(batch)
                return

            batch.append(item)
            if previous_start is not None:
                self._check_for_unsorted_objects(item, previous_start)
            previous_start = item.start

        self._fill(batch, exhausted=True)

    def _check_for_unsorted_objects(self, item: ObjectT, previous_start: int) -> None:
        """Count an object that starts before the one before it, and give up on too many.

        Args:
            item (ObjectT): The object just buffered.
            previous_start (int): Where the object before it started.

        Raises:
            RuntimeError: If more than ``max_not_ordered_variants`` are far enough out of order.
        """

        if item.start - previous_start >= 0:
            return
        if abs(item.start - previous_start) <= self.max_not_ordered_variants_distance:
            return

        self._not_ordered_count += 1
        if self._not_ordered_count > self.max_not_ordered_variants:
            raise RuntimeError(
                "Exiting the program. The input file appears to be unsorted. "
                "Please sort by chromosome and by location and re-submit."
            )

    @staticmethod
    def _seq_region_name(item: ObjectT) -> str:
        """Return the sequence region an object belongs to.

        Args:
            item (ObjectT): The object to place.

        Returns:
            str: The name of its sequence region.

        Raises:
            ValueError: If the object carries no slice.
        """

        if item.slice is None:
            raise ValueError("a buffered object needs a slice to be grouped by sequence region")
        return item.slice.seq_region_name
