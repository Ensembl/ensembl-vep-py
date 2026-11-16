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
"""Interval lookup indexes."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Generic, Protocol, TypeVar, runtime_checkable

from attrs import define, field, frozen

__all__ = ["Interval", "IntervalIndex", "IntervalTreeIndex"]

ObjectT = TypeVar("ObjectT")


@frozen(kw_only=True, auto_attribs=True)
class Interval:
    """Represent an inclusive genomic interval or interbase insertion point.

    ``start=end+1`` denotes an insertion before ``start``. The interval index
    expands that point to its two flanking bases only for candidate retrieval.
    Matchers retain the zero-width insertion semantics.
    """

    start: int
    end: int

    def __attrs_post_init__(self) -> None:
        """Validate the interval representation."""

        if self.start < 1:
            raise ValueError("interval start must be at least 1")
        if self.end < self.start - 1:
            raise ValueError("interval end must not be less than start - 1")

    @property
    def is_insertion(self) -> bool:
        """Return whether the interval is a zero-width insertion point."""

        return self.start == self.end + 1

    @property
    def lookup_start(self) -> int:
        """Return the inclusive start used for candidate retrieval."""

        return max(1, self.end) if self.is_insertion else self.start

    @property
    def lookup_end(self) -> int:
        """Return the inclusive end used for candidate retrieval."""

        return self.start if self.is_insertion else self.end


@runtime_checkable
class IntervalIndex(Protocol[ObjectT]):
    """Define an interchangeable interval-candidate index."""

    def add(self, value: ObjectT, interval: Interval) -> None:
        """Add an object under an interval."""

        raise NotImplementedError

    def extend(self, values: Iterable[tuple[ObjectT, Interval]]) -> None:
        """Add multiple objects under their intervals."""

        raise NotImplementedError

    def get_objects_by_intervals(self, intervals: Iterable[Interval]) -> list[ObjectT]:
        """Return objects meeting any candidate interval."""

        raise NotImplementedError

    def clear(self) -> None:
        """Clear indexed objects."""

        raise NotImplementedError


@frozen(kw_only=True, auto_attribs=True)
class _Entry(Generic[ObjectT]):
    """Store one indexed object and its candidate-retrieval bounds."""

    value: ObjectT
    interval: Interval
    order: int


@define(kw_only=True, auto_attribs=True)
class _Node(Generic[ObjectT]):
    """Store one centred interval-tree node."""

    center: int
    entries: tuple[_Entry[ObjectT], ...]
    left: _Node[ObjectT] | None = None
    right: _Node[ObjectT] | None = None


@define(kw_only=True, auto_attribs=True)
class IntervalTreeIndex(Generic[ObjectT]):
    """Index objects by inclusive intervals.

    The tree is rebuilt lazily after writes. This keeps index construction cheap
    for an input buffer while serving repeated feature queries efficiently.
    """

    _entries: list[_Entry[ObjectT]] = field(factory=list, init=False, repr=False)
    _root: _Node[ObjectT] | None = field(default=None, init=False, repr=False)
    _dirty: bool = field(default=False, init=False, repr=False)

    def add(self, value: ObjectT, interval: Interval) -> None:
        """Add an object under an interval."""

        self._entries.append(_Entry(value=value, interval=interval, order=len(self._entries)))
        self._dirty = True

    def extend(self, values: Iterable[tuple[ObjectT, Interval]]) -> None:
        """Add multiple objects under their intervals."""

        for value, interval in values:
            self._entries.append(_Entry(value=value, interval=interval, order=len(self._entries)))
        self._dirty = True

    def get_objects_by_intervals(self, intervals: Iterable[Interval]) -> list[ObjectT]:
        """Return objects that intersect any candidate-retrieval interval."""

        self._ensure_tree()
        if self._root is None:
            return []

        found: dict[int, _Entry[ObjectT]] = {}
        for interval in intervals:
            for entry in self._query(self._root, interval):
                entry_id = id(entry.value)
                found.setdefault(entry_id, entry)
        return [entry.value for entry in sorted(found.values(), key=lambda entry: entry.order)]

    def clear(self) -> None:
        """Clear indexed objects."""

        self._entries.clear()
        self._root = None
        self._dirty = False

    def _ensure_tree(self) -> None:
        """Build the tree after indexed objects change."""

        if not self._dirty:
            return
        self._root = self._build(self._entries)
        self._dirty = False

    @classmethod
    def _build(cls, entries: list[_Entry[ObjectT]]) -> _Node[ObjectT] | None:
        """Build a centred tree from entries."""

        if not entries:
            return None

        centers = sorted((entry.interval.lookup_start + entry.interval.lookup_end) // 2 for entry in entries)
        center = centers[len(centers) // 2]
        left_entries: list[_Entry[ObjectT]] = []
        right_entries: list[_Entry[ObjectT]] = []
        center_entries: list[_Entry[ObjectT]] = []

        for entry in entries:
            if entry.interval.lookup_end < center:
                left_entries.append(entry)
            elif entry.interval.lookup_start > center:
                right_entries.append(entry)
            else:
                center_entries.append(entry)

        return _Node(
            center=center,
            entries=tuple(center_entries),
            left=cls._build(left_entries),
            right=cls._build(right_entries),
        )

    @classmethod
    def _query(cls, node: _Node[ObjectT] | None, interval: Interval) -> list[_Entry[ObjectT]]:
        """Return matching entries beneath one node."""

        if node is None:
            return []

        found: list[_Entry[ObjectT]] = []
        query_start = interval.lookup_start
        query_end = interval.lookup_end
        for entry in node.entries:
            if entry.interval.lookup_start <= query_end and entry.interval.lookup_end >= query_start:
                found.append(entry)

        if query_start <= node.center:
            found.extend(cls._query(node.left, interval))
        if query_end >= node.center:
            found.extend(cls._query(node.right, interval))
        return found
