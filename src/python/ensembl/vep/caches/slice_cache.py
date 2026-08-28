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
"""Slice cache."""

from __future__ import annotations

from threading import RLock
from typing import TYPE_CHECKING, Protocol, runtime_checkable
from weakref import WeakValueDictionary

from attrs import define, field, frozen

from ensembl.vep.common.enums import CoordSystemLevel, Strand
from ensembl.vep.models.sequence.coord_system import CoordSystem
from ensembl.vep.models.sequence.slice import Slice

if TYPE_CHECKING:
    from ensembl.vep.sources.sequence_source import SequenceSource

_CHROMOSOME_COORD_SYSTEM = CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME, default=True)


@frozen(kw_only=True, auto_attribs=True)
class SliceKey:
    """Identity of a cached slice."""

    sequence_source_id: int | None
    seq_region_name: str
    start: int
    end: int
    strand: Strand


@runtime_checkable
class SliceCache(Protocol):
    """Interface for caches that construct and share Slice objects."""

    def get_slice(
        self,
        seq_region_name: str,
        sequence_source: SequenceSource | None = None,
        *,
        start: int = 1,
        end: int | None = None,
        strand: Strand = Strand.FORWARD,
    ) -> Slice:
        """Return the matching slice, constructing it on a cache miss."""

        raise NotImplementedError

    def clear(self) -> None:
        """Clear cached state."""

        raise NotImplementedError


@define(kw_only=True, auto_attribs=True)
class WeakSliceCache:
    """Share slices without extending their lifetime.

    Cached slices remain available while a feature or another consumer holds a strong reference.
    The cache entry is removed automatically after the last strong reference is released.
    """

    _cache: WeakValueDictionary[SliceKey, Slice] = field(factory=WeakValueDictionary, init=False, repr=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    def get_slice(
        self,
        seq_region_name: str,
        sequence_source: SequenceSource | None = None,
        *,
        start: int = 1,
        end: int | None = None,
        strand: Strand = Strand.FORWARD,
    ) -> Slice:
        """Return the matching slice, constructing it on a cache miss.

        If ``end`` is omitted, the complete sequence-region length is loaded from
        ``sequence_source``. A slice cannot be constructed when neither an explicit end nor a
        resolvable sequence-region length is available.

        Args:
            seq_region_name: Name of the sequence region.
            sequence_source: Sequence provider used to resolve a missing end and attached to a
                newly built slice.
            start: 1-based inclusive slice start.
            end: 1-based inclusive slice end. Defaults to the complete sequence-region length.
            strand: Slice orientation relative to the sequence region.

        Returns:
            The shared Slice for the requested identity.

        Raises:
            ValueError: If the coordinates are invalid or the sequence-region length cannot be resolved.
        """

        resolved_end = self._resolve_end(seq_region_name, sequence_source, end)
        if start < 1:
            raise ValueError("slice start must be 1-based")
        if resolved_end < start:
            raise ValueError("slice start must be less than or equal to slice end")

        key = SliceKey(
            sequence_source_id=id(sequence_source) if sequence_source is not None else None,
            seq_region_name=seq_region_name,
            start=start,
            end=resolved_end,
            strand=strand,
        )

        with self._lock:
            cached = self._cache.get(key)
            if cached is not None:
                return cached

        loaded = self._load_slice(
            seq_region_name,
            sequence_source,
            start=start,
            end=resolved_end,
            strand=strand,
        )

        with self._lock:
            cached = self._cache.get(key)
            if cached is not None:
                return cached
            self._cache[key] = loaded

        return loaded

    @staticmethod
    def _resolve_end(
        seq_region_name: str,
        sequence_source: SequenceSource | None,
        end: int | None,
    ) -> int:
        """Resolve an explicit end or load the complete sequence-region length."""

        if end is not None:
            return end
        if sequence_source is None:
            raise ValueError(
                f"Cannot build slice for {seq_region_name}: no end or sequence source was provided"
            )

        seq_region_length = sequence_source.get_seq_region_length(seq_region_name)
        if seq_region_length is None:
            raise ValueError(f"Cannot build slice for {seq_region_name}: sequence-region length is unknown")
        return seq_region_length

    @staticmethod
    def _load_slice(
        seq_region_name: str,
        sequence_source: SequenceSource | None,
        *,
        start: int,
        end: int,
        strand: Strand,
    ) -> Slice:
        """Construct a slice for a cache miss."""

        return Slice(
            coord_system=_CHROMOSOME_COORD_SYSTEM,
            seq_region_name=seq_region_name,
            start=start,
            end=end,
            strand=strand,
            sequence_source=sequence_source,
        )

    def clear(self) -> None:
        """Clear all live cache entries."""

        with self._lock:
            self._cache.clear()
