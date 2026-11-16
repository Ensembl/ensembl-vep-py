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
from collections.abc import Iterable

from attrs import define

from ensembl.vep.buffer import InputBuffer
from ensembl.vep.caches import Interval


@define(kw_only=True, auto_attribs=True)
class _Input:
    """Minimal input object with interval coordinates."""

    start: int
    end: int


@define(kw_only=True, auto_attribs=True)
class _DummyIndex:
    """Minimal injected index for buffer construction tests."""

    indexed: list[tuple[_Input, Interval]]
    cleared: bool = False

    def add(self, value: _Input, interval: Interval) -> None:
        """Add one indexed object."""

        self.indexed.append((value, interval))

    def extend(self, values: Iterable[tuple[_Input, Interval]]) -> None:
        """Record construction values."""

        self.indexed.extend(values)

    def get_objects_by_intervals(self, intervals: Iterable[Interval]) -> list[_Input]:
        """Return injected objects independently of query coordinates."""

        del intervals
        return [item for item, _ in self.indexed]

    def clear(self) -> None:
        """Record a clear request."""

        self.cleared = True


class TestInputBuffer:
    """Input-object indexing."""

    def test_returns_objects_for_query_intervals(self) -> None:
        """Retrieve input objects that meet a feature interval."""

        first = _Input(start=1, end=10)
        second = _Input(start=20, end=30)
        buffer = InputBuffer(objects=[first, second])

        objects = buffer.get_objects_by_intervals([Interval(start=5, end=20)])

        assert objects == [first, second]

    def test_returns_boundary_objects_for_an_insertion_query(self) -> None:
        """Retrieve both input flanks for an insertion candidate query."""

        left = _Input(start=1, end=9)
        right = _Input(start=10, end=20)
        buffer = InputBuffer(objects=[left, right])

        objects = buffer.get_objects_by_intervals([Interval(start=10, end=9)])

        assert objects == [left, right]

    def test_clear_removes_the_input_index(self) -> None:
        """Clear candidates without mutating buffered input objects."""

        item = _Input(start=1, end=10)
        buffer = InputBuffer(objects=[item])
        buffer.clear()

        assert buffer.objects == (item,)
        assert buffer.get_objects_by_intervals([Interval(start=1, end=10)]) == []

    def test_round_trips_through_pickle(self) -> None:
        """Preserve indexed candidates across worker-spec serialisation."""

        item = _Input(start=1, end=10)
        buffer = InputBuffer(objects=[item])

        restored = pickle.loads(pickle.dumps(buffer))

        assert restored.get_objects_by_intervals([Interval(start=1, end=10)]) == [item]

    def test_uses_an_injected_interval_index(self) -> None:
        """Use the supplied index implementation instead of the default tree."""

        item = _Input(start=1, end=10)
        index = _DummyIndex(indexed=[])
        buffer = InputBuffer(objects=[item], interval_index=index)

        assert buffer.get_objects_by_intervals([Interval(start=100, end=100)]) == [item]
        buffer.clear()
        assert index.cleared
