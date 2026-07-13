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
from typing import Any, TYPE_CHECKING

from attrs import define, field

from ensembl.vep.common import Strand, invalidate, invalidates
from ensembl.vep.models.sequence.coord_system import CoordSystem

if TYPE_CHECKING:
    from ensembl.vep.sources.sequence_source import SequenceSource


@define(kw_only=True, auto_attribs=True)
class Slice:
    """Represents a genomic region on a sequence.

    Attributes:
        coord_system (CoordSystem): Coordinate system.
        seq_region_name (str): Sequence region name.
        start (int): Start position relative to sequence region.
        end (int): End position relative to sequence region.
        strand (Strand): Slice orientation relative to sequence region.
        seq_region_length (int | None): Sequence region length.
        sequence_source (SequenceSource | None): Lazy sequence provider.
        seq (Any | None): Manually assigned sequence override.
    """

    coord_system: CoordSystem
    seq_region_name: str = field(on_setattr=invalidates("_seq", "_resolved_seq_region_length"))
    start: int = field(on_setattr=invalidates("_seq"))
    end: int = field(on_setattr=invalidates("_seq"))
    strand: Strand = field(on_setattr=invalidates("_seq"))
    _seq_region_length: int | None = field(default=None, repr=False, alias="seq_region_length")
    _seq_override: Any | None = field(default=None, repr=False, alias="seq")
    sequence_source: SequenceSource | None = field(
        default=None, repr=False, on_setattr=invalidates("_seq", "_resolved_seq_region_length")
    )
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    @property
    def seq(self) -> Any | None:
        """Return the resolved sequence or a manually assigned override.

        Returns:
            Any | None: The resolved sequence, or None if unavailable.
        """

        if self._seq_override is not None:
            return self._seq_override

        return self._seq

    @seq.setter
    def seq(self, value: Any | None) -> None:
        """Set or clear the sequence override.

        Args:
            value (Any | None): Sequence override value.
        """

        self._seq_override = value
        invalidate(self, "_seq")

    @property
    def seq_region_length(self) -> int | None:
        """Return the resolved sequence region length.

        Returns:
            int | None: The sequence region length, or None if unavailable.
        """

        if self._seq_region_length is not None:
            return self._seq_region_length

        return self._resolved_seq_region_length

    @seq_region_length.setter
    def seq_region_length(self, value: int | None) -> None:
        """Set the explicit sequence region length.

        Args:
            value (int | None): Sequence region length.
        """

        self._seq_region_length = value
        invalidate(self, "_resolved_seq_region_length")

    @cached_property
    def _resolved_seq_region_length(self) -> int | None:
        """Resolve the sequence region length from the attached source.

        Returns:
            int | None: The sequence region length, or None if unavailable.
        """

        if self.sequence_source is None:
            return None

        return self.sequence_source.get_seq_region_length(self.seq_region_name)

    @cached_property
    def _seq(self) -> Any | None:
        """Resolve the sequence lazily from the attached source.

        Returns:
            Any | None: The resolved sequence, or None if unavailable.
        """

        if self.sequence_source is None:
            return None

        return self.sequence_source.get_sequence_by_start_end_strand(
            self.seq_region_name,
            self.start,
            self.end,
            strand=self.strand,
        )

    def __deepcopy__(self, memo: dict[int, Any]) -> Slice:
        """Create a deepcopy without copying the internal lock.

        Args:
            memo (dict[int, Any]): Deepcopy memoization map.

        Returns:
            Slice: The copied slice instance.
        """

        copied = type(self)(
            coord_system=deepcopy(self.coord_system, memo),
            seq_region_name=self.seq_region_name,
            start=self.start,
            end=self.end,
            strand=self.strand,
            seq_region_length=deepcopy(self._seq_region_length, memo),
            sequence_source=self.sequence_source,
            seq=deepcopy(self._seq_override, memo),
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
