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
"""Tests for input buffering."""

from __future__ import annotations

import pickle

import pytest
from attrs import define

from ensembl.vep.buffer.input_buffer import InputBuffer
from ensembl.vep.caches import Interval


@define(kw_only=True, auto_attribs=True)
class _Input:
    """Minimal input object with interval coordinates."""

    start: int
    end: int


class TestInputBuffer:
    """Input-object indexing."""

    def test_returns_objects_for_query_intervals(self) -> None:
        """Retrieve input objects that meet a feature interval."""
        first = _Input(start=1, end=10)
        second = _Input(start=20, end=30)
        buffer = InputBuffer(objects=(first, second))

        objects = buffer.get_objects_by_intervals([Interval(start=5, end=20)])

        assert objects == [first, second]

    def test_reports_each_object_once_across_query_intervals(self) -> None:
        """Report an object met by several intervals once, in input order."""
        first = _Input(start=1, end=10)
        second = _Input(start=5, end=30)
        buffer = InputBuffer(objects=(first, second))

        objects = buffer.get_objects_by_intervals([Interval(start=6, end=7), Interval(start=8, end=9)])

        assert objects == [first, second]

    def test_returns_boundary_objects_for_an_insertion_query(self) -> None:
        """Retrieve both input flanks for an insertion candidate query."""
        left = _Input(start=1, end=9)
        right = _Input(start=10, end=20)
        buffer = InputBuffer(objects=(left, right))

        objects = buffer.get_objects_by_intervals([Interval(start=10, end=9)])

        assert objects == [left, right]

    def test_clamps_an_insertion_before_the_first_base(self) -> None:
        """Retrieve the first input object for an insertion before the first base."""
        first = _Input(start=1, end=1)
        buffer = InputBuffer(objects=(first,))

        assert buffer.get_objects_by_intervals([Interval(start=1, end=0)]) == [first]

    def test_returns_objects_overlapping_one_interval(self) -> None:
        """Retrieve input objects that overlap a plain coordinate pair."""
        item = _Input(start=1, end=10)
        buffer = InputBuffer(objects=(item,))

        assert buffer.get_overlapping_objects(start=5, end=5) == [item]
        assert buffer.get_overlapping_objects(start=11, end=12) == []

    def test_an_empty_batch_meets_nothing(self) -> None:
        """Answer queries against an empty batch without indexing anything."""
        buffer: InputBuffer[_Input] = InputBuffer()

        assert buffer.get_objects_by_intervals([Interval(start=1, end=10)]) == []
        assert buffer.get_overlapping_objects(start=1, end=10) == []

    def test_clear_drops_the_index_and_keeps_the_batch(self) -> None:
        """Clear the index without mutating buffered input objects."""
        item = _Input(start=1, end=10)
        buffer = InputBuffer(objects=(item,))
        buffer.get_objects_by_intervals([Interval(start=1, end=10)])

        buffer.clear()

        assert buffer.objects == (item,)
        assert buffer.get_objects_by_intervals([Interval(start=1, end=10)]) == [item]

    def test_is_empty_until_the_input_is_spent(self) -> None:
        """Report a buffer as empty only once its input is exhausted."""
        buffer: InputBuffer[_Input] = InputBuffer()

        assert not buffer.is_empty

    def test_refilling_needs_a_policy(self) -> None:
        """Leave the refill policy to a subclass."""
        buffer: InputBuffer[_Input] = InputBuffer()

        with pytest.raises(NotImplementedError):
            buffer.next(iter(()))

    def test_round_trips_through_pickle(self) -> None:
        """Preserve indexed candidates across worker-spec serialisation."""
        item = _Input(start=1, end=10)
        buffer = InputBuffer(objects=(item,))
        buffer.get_objects_by_intervals([Interval(start=1, end=10)])

        restored = pickle.loads(pickle.dumps(buffer))

        assert restored.get_objects_by_intervals([Interval(start=1, end=10)]) == [item]
