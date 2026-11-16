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
"""Tests for interval lookup indexes."""

from __future__ import annotations

import pytest

from ensembl.vep.caches import Interval, IntervalTreeIndex


class TestInterval:
    """Interval coordinate validation."""

    def test_identifies_insertion_point(self) -> None:
        """Identify an inverted one-base span as an insertion."""

        interval = Interval(start=10, end=9)

        assert interval.is_insertion
        assert interval.lookup_start == 9
        assert interval.lookup_end == 10

    @pytest.mark.parametrize(
        ("start", "end"),
        [
            pytest.param(0, 0, id="zero_start"),
            pytest.param(10, 8, id="more_than_one_inverted_base"),
        ],
    )
    def test_rejects_invalid_coordinates(self, start: int, end: int) -> None:
        """Reject coordinates outside the normal or insertion representation."""

        with pytest.raises(ValueError):
            Interval(start=start, end=end)


class TestIntervalTreeIndex:
    """Candidate retrieval."""

    def test_returns_intersecting_objects_once(self) -> None:
        """Return each object once for overlapping query intervals."""

        index = IntervalTreeIndex[str]()
        index.extend(
            [
                ("left", Interval(start=1, end=10)),
                ("middle", Interval(start=10, end=20)),
                ("right", Interval(start=30, end=40)),
            ]
        )

        objects = index.get_objects_by_intervals([Interval(start=10, end=10), Interval(start=15, end=15)])

        assert objects == ["left", "middle"]

    def test_expands_insertion_for_candidate_retrieval(self) -> None:
        """Find features on both flanks of an insertion point."""

        index = IntervalTreeIndex[str]()
        index.extend(
            [
                ("left", Interval(start=1, end=9)),
                ("right", Interval(start=10, end=20)),
            ]
        )

        objects = index.get_objects_by_intervals([Interval(start=10, end=9)])

        assert objects == ["left", "right"]

    def test_clamps_insertion_before_first_base_for_candidate_retrieval(self) -> None:
        """Find the first feature for an insertion before the first base."""

        index = IntervalTreeIndex[str]()
        index.add("first", Interval(start=1, end=1))

        objects = index.get_objects_by_intervals([Interval(start=1, end=0)])

        assert objects == ["first"]

    def test_clear_removes_indexed_objects(self) -> None:
        """Remove all values after clearing the index."""

        index = IntervalTreeIndex[str]()
        index.add("feature", Interval(start=1, end=10))
        index.clear()

        assert index.get_objects_by_intervals([Interval(start=1, end=10)]) == []
