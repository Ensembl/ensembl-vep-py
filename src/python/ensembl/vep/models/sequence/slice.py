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
"""Slice class."""

from __future__ import annotations

from copy import deepcopy
from functools import cached_property
from threading import RLock
from typing import TYPE_CHECKING

from attrs import define, field

from ensembl.vep.common import Strand, UnsafeSetterMixin, invalidates
from ensembl.vep.models.sequence.coord_system import CoordSystem

if TYPE_CHECKING:
    from ensembl.vep.sources.sequence_source import SequenceSource


@define(kw_only=True, auto_attribs=True)
class Slice(UnsafeSetterMixin):
    """Represents a genomic region on a sequence.

    Attributes:
        coord_system (CoordSystem): Coordinate system.
        seq_region_name (str): Sequence region name.
        start (int): Start position relative to sequence region.
        end (int): End position relative to sequence region.
        strand (Strand): Slice orientation relative to sequence region.
        seq (str | None): Sequence of the slice.
        seq_region_length (int | None): Sequence region length.
        sequence_source (SequenceSource | None): Lazy sequence provider.
    """

    coord_system: CoordSystem
    seq_region_name: str = field(on_setattr=invalidates("_seq", "_seq_region_length"))
    start: int = field(on_setattr=invalidates("_seq"))
    end: int = field(on_setattr=invalidates("_seq"))
    strand: Strand = field(on_setattr=invalidates("_seq"))
    sequence_source: SequenceSource | None = field(
        default=None, repr=False, on_setattr=invalidates("_seq", "_seq_region_length")
    )
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    @property
    def seq(self) -> str | None:
        """Public accessor for the _seq

        Returns:
            str | None: The resolved sequence, or None if unavailable.
        """
        return self._seq

    @cached_property
    def _seq(self) -> str | None:
        """Resolve the sequence lazily from the attached sequence source.

        Returns:
            str | None: The resolved sequence, or None if unavailable.
        """

        with self._lock:
            if self.sequence_source is None:
                return None

            return self.sequence_source.get_sequence_by_start_end_strand(
                self.seq_region_name,
                self.start,
                self.end,
                strand=self.strand,
            )

    @property
    def seq_region_length(self) -> int | None:
        """Public accessor for the _seq_region_length.

        Returns:
            int | None: The sequence region length, or None if unavailable.
        """
        return self._seq_region_length

    @cached_property
    def _seq_region_length(self) -> int | None:
        """Resolve the sequence region length from the attached sequence source.

        Returns:
            int | None: The sequence region length, or None if unavailable.
        """

        with self._lock:
            if self.sequence_source is None:
                return None

            return self.sequence_source.get_seq_region_length(self.seq_region_name)

    def __deepcopy__(self, memo: dict[int, object]) -> Slice:
        """Create a deepcopy without copying the internal lock.

        Args:
            memo (dict[int, object]): Deepcopy memoization map.

        Returns:
            Slice: The copied slice instance.
        """

        copied = type(self)(
            coord_system=deepcopy(self.coord_system, memo),
            seq_region_name=self.seq_region_name,
            start=self.start,
            end=self.end,
            strand=self.strand,
            sequence_source=self.sequence_source,
        )

        memo[id(self)] = copied
        return copied

    def __eq__(self, other: object) -> bool:
        """Compare two Slice objects and define if they are identical."""

        if type(self) is not type(other):
            return False

        return (
            self.coord_system == other.coord_system
            and self.seq_region_name == other.seq_region_name
            and self.start == other.start
            and self.end == other.end
            and self.strand == other.strand
        )

    def is_circular(self) -> bool:
        """Indicate whether the slice is circular.

        Returns:
            bool: False for regular slices.
        """

        return False
