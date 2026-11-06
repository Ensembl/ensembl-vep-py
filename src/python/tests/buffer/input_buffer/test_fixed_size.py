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
"""Tests for the fixed-size input buffer."""

# pylint: disable=redefined-outer-name,protected-access

from __future__ import annotations

from collections import deque

import pytest

from ensembl.vep.buffer.input_buffer import FixedSizeInputBuffer
from ensembl.vep.common import Strand
from ensembl.vep.models.variation import Variation


@pytest.fixture
def buffer() -> FixedSizeInputBuffer[Variation]:
    """Return a buffer holding five objects at most, intolerant of unsorted input."""
    return FixedSizeInputBuffer(buffer_size=5, max_not_ordered_variants=2)


class TestInitialization:
    """How a fixed-size buffer starts out."""

    def test_initialization_without_keywords(self) -> None:
        """Refuse positional construction."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = FixedSizeInputBuffer(5)  # type: ignore[arg-type, misc]

    def test_initial_attributes(self, buffer: FixedSizeInputBuffer[Variation]) -> None:
        """Start with an empty batch, nothing held back, and input to come."""
        assert buffer.objects == ()
        assert buffer._pre_buffer == deque()
        assert not buffer.is_empty


class TestNext:
    """Refilling the batch."""

    def test_empty_input(self, buffer: FixedSizeInputBuffer[Variation]) -> None:
        """Report spent input on the first refill when there is none."""
        buffer.next(iter(()))

        assert buffer.objects == ()
        assert buffer.is_empty

    def test_fills_to_the_buffer_size_then_to_the_region_boundary(
        self, buffer: FixedSizeInputBuffer[Variation], variation_list: list[Variation]
    ) -> None:
        """Take at most ``buffer_size`` objects, and never two sequence regions at once."""
        variations = iter(variation_list)

        buffer.next(variations)
        assert [variation.start for variation in buffer.objects] == [1, 2, 3, 4, 5]
        assert [held.start for held in buffer._pre_buffer] == [6]

        buffer.next(variations)
        assert [variation.start for variation in buffer.objects] == [6, 7, 8]
        assert [held.start for held in buffer._pre_buffer] == [10]
        assert not buffer.is_empty

        buffer.next(variations)
        assert [variation.start for variation in buffer.objects] == [10, 11, 12, 13]
        assert buffer._pre_buffer == deque()
        assert buffer.is_empty

    def test_each_batch_holds_one_sequence_region(
        self, buffer: FixedSizeInputBuffer[Variation], variation_list: list[Variation]
    ) -> None:
        """Group each batch by the sequence region its first object is on."""
        variations = iter(variation_list)
        buffer.next(variations)
        buffer.next(variations)

        regions = {variation.slice.seq_region_name for variation in buffer.objects if variation.slice}

        assert regions == {"1"}

    def test_stops_when_the_held_back_objects_fill_the_batch(
        self, buffer: FixedSizeInputBuffer[Variation], variation_list: list[Variation]
    ) -> None:
        """Fill a batch from held-back objects alone, reading no further input."""
        buffer._pre_buffer = deque(variation_list)

        buffer.next(iter(()))

        assert len(buffer.objects) == 5
        assert len(buffer._pre_buffer) == 7

    def test_stops_at_a_region_change_among_the_held_back_objects(
        self, buffer: FixedSizeInputBuffer[Variation], variation_list: list[Variation]
    ) -> None:
        """Stop a batch where the held-back objects change sequence region."""
        buffer._pre_buffer = deque(variation_list[0:3] + variation_list[8:12])

        buffer.next(iter(()))

        assert len(buffer.objects) == 3
        assert len(buffer._pre_buffer) == 4

    def test_rejects_unsorted_input(
        self, buffer: FixedSizeInputBuffer[Variation], variation_list_out_of_order: list[Variation]
    ) -> None:
        """Give up on input that jumps backwards too often."""
        with pytest.raises(RuntimeError, match="The input file appears to be unsorted"):
            buffer.next(iter(variation_list_out_of_order))

    def test_tolerates_unsorted_input_below_the_count_threshold(
        self, buffer: FixedSizeInputBuffer[Variation], variation_list_out_of_order: list[Variation]
    ) -> None:
        """Accept backward jumps while they stay within the tolerated count."""
        buffer.max_not_ordered_variants = 5

        buffer.next(iter(variation_list_out_of_order))

        assert len(buffer.objects) == 5

    def test_tolerates_unsorted_input_below_the_distance_threshold(
        self, buffer: FixedSizeInputBuffer[Variation], variation_list_out_of_order: list[Variation]
    ) -> None:
        """Ignore a backward jump shorter than the tolerated distance."""
        buffer.max_not_ordered_variants_distance = 100

        buffer.next(iter(variation_list_out_of_order))

        assert len(buffer.objects) == 5

    def test_refuses_an_object_without_a_slice(self, buffer: FixedSizeInputBuffer[Variation]) -> None:
        """Refuse an object that cannot be placed on a sequence region."""
        unplaced = Variation(start=1, end=1, strand=Strand.FORWARD, allele_string="A/T")

        with pytest.raises(ValueError, match="needs a slice"):
            buffer.next(iter([unplaced]))


class TestOverlapQueries:
    """Querying a refilled batch."""

    def test_empty_batch_overlaps_nothing(self, buffer: FixedSizeInputBuffer[Variation]) -> None:
        """Answer an overlap query against an empty batch with nothing."""
        assert buffer.get_overlapping_objects(start=1, end=10) == []

    def test_returns_the_overlapping_objects(
        self, buffer: FixedSizeInputBuffer[Variation], variation_list: list[Variation]
    ) -> None:
        """Return the buffered objects an interval meets."""
        buffer.next(iter(variation_list))

        assert buffer.get_overlapping_objects(start=2, end=3) == variation_list[1:3]
        assert buffer.get_overlapping_objects(start=6, end=6) == []

    def test_reads_a_reversed_interval_as_an_insertion_point(
        self, buffer: FixedSizeInputBuffer[Variation], variation_list: list[Variation]
    ) -> None:
        """Swap a reversed interval rather than refusing it."""
        buffer.next(iter(variation_list))

        assert buffer.get_overlapping_objects(start=3, end=2) == variation_list[1:3]

    def test_a_refill_discards_the_previous_index(
        self, buffer: FixedSizeInputBuffer[Variation], variation_list: list[Variation]
    ) -> None:
        """Index the current batch, not the one before it."""
        variations = iter(variation_list)
        buffer.next(variations)
        buffer.get_overlapping_objects(start=1, end=5)

        buffer.next(variations)

        assert buffer.get_overlapping_objects(start=1, end=5) == []
        assert buffer.get_overlapping_objects(start=6, end=8) == variation_list[5:8]
