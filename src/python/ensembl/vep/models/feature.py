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
"""Feature class."""

from __future__ import annotations

import logging
from typing import Any
from functools import cached_property
from threading import Lock
from attrs import define, field

from ensembl.vep.models.sequence import Slice
from ensembl.vep.models.sequence.coord_system import CoordSystem
from ensembl.vep.common import Strand, invalidates

logger = logging.getLogger(__name__)


@define(kw_only=True, auto_attribs=True)
class Feature:
    """
    Represents an Ensembl-specific sequence feature describing a location on a sequence in an arbitrary coordinate system.

    Args:
        start (int): Start position relative to slice.
        end (int): End position relative to slice.
        strand (Strand): Feature orientation relative to slice.
        slice (Slice | None): Slice object describing the feature context.
        seq (Any | None): Associated biological sequence.

    Attributes:
        start (int): Start position relative to slice.
        end (int): End position relative to slice.
        strand (Strand): Feature orientation relative to slice.
        slice (Slice | None): Slice object describing the feature context.
        seq (Any | None): Associated biological sequence.
    """

    start: int = field(on_setattr=invalidates("__length", "__seq_region_start"))
    end: int = field(on_setattr=invalidates("__length", "__seq_region_end"))
    strand: Strand = field(on_setattr=invalidates("__seq_region_strand"))
    slice: Slice | None = field(
        default=None, on_setattr=invalidates("__seq_region_start", "__seq_region_end", "__seq_region_strand")
    )
    seq: Any | None = field(default=None)
    _lock: Lock = field(factory=Lock, init=False, repr=False)

    @property
    def length(self) -> int:
        """
        Public accessor for the __length.

        Returns:
            int: The length of the feature.
        """
        return self.__length

    @cached_property
    def __length(self) -> int:
        """
        Returns the length of the feature.

        Returns:
            int: The length of the feature.
        """
        with self._lock:
            return self.end - self.start + 1

    @property
    def seq_region_start(self) -> int:
        """
        Public accessor for the __seq_region_start.

        Returns:
            int: The start position in the sequence region.
        """
        return self.__seq_region_start

    @cached_property
    def __seq_region_start(self) -> int:
        """
        Returns the start position in the sequence region.

        Returns:
            int: The start position in the sequence region.
        """
        with self._lock:
            if self.slice is None:
                logging.warning("Cannot obtain seq_region_start for a feature without attached slice.")
                return None
            return self.slice.start + self.start

    @property
    def seq_region_end(self) -> int:
        """
        Public accessor for the __seq_region_end.

        Returns:
            int: The end position in the sequence region.
        """
        return self.__seq_region_end

    @cached_property
    def __seq_region_end(self) -> int:
        """
        Returns the end position in the sequence region.

        Returns:
            int: The end position in the sequence region.
        """
        with self._lock:
            if self.slice is None:
                logging.warning("Cannot obtain seq_region_end for a feature without attached slice.")
                return None
            return self.slice.start + self.end

    @property
    def seq_region_strand(self) -> Strand:
        """
        Public accessor for the __seq_region_strand.

        Returns:
            Strand: The strand in the sequence region.
        """
        return self.__seq_region_strand

    @cached_property
    def __seq_region_strand(self) -> Strand:
        """
        Returns the strand in the sequence region.

        Returns:
            Strand: The strand in the sequence region.
        """
        with self._lock:
            if self.slice is None:
                logging.warning("Cannot obtain seq_region_strand for a feature without attached slice.")
                return None
        return self.slice.strand * self.strand

    def __eq__(self, other: object) -> bool:
        """
        Compares two features using various criteria.

        Args:
            other (object): Another object to compare against.

        Returns:
            bool: True if equal, False otherwise.
        """
        if type(self) is not type(other):
            return False

        return (
            self.start == other.start
            and self.end == other.end
            and self.strand == other.strand
            and self.slice == other.slice
        )

    def transform(self, coord_system: CoordSystem) -> Feature | None:
        """
        Transforms the feature to a new coordinate system.

        Args:
            coord_system (CoordSystem): The coordinate system to transform to.

        Returns:
            Feature | None: The transformed feature or None if transformation is not possible.
        """
        return None

    def transfer(self, slice: Slice = None) -> Feature | None:
        """
        Transfers the feature to a new slice.

        Args:
            slice (Slice, optional): The slice to transfer to. If None, uses current slice.

        Returns:
            Feature | None: The transferred feature or None if transfer is not possible.
        """
        return None

    def overlaps(self, other: Feature) -> bool | None:
        """
        Checks if this feature overlaps with another feature.

        Args:
            other (Feature): Another feature to check overlap with.

        Returns:
            bool | None: True if overlaps, False if not, None if slices are missing.
        """
        if self.slice is None or other.slice is None:
            logging.warning("Cannot compare two features is any of them is without attached slice.")
            return None

        if self.slice.seq_region_name != other.slice.seq_region_name:
            logging.warning(
                f"Cannot compare two features on different sequence regions. {self.__class__.__name__} is on "
                f"{self.slice.seq_region_name} and {other.__class__.__name__} is on {other.slice.seq_region_name}."
            )
            return False

        return self.seq_region_start <= other.seq_region_end and self.seq_region_end >= other.seq_region_start

    def feature_slice(self) -> Slice | None:
        """
        Returns a slice that covers the area of this feature.

        The feature start will be at 1 on the returned slice, and it will have the length of this feature.

        Returns:
            Slice | None: A Slice object covering the feature area, or None if the feature has no attached slice.

        Raises:
            Warning: If the feature has no attached slice, a warning is logged and None is returned.
        """
        slice = self.slice

        if slice is None:
            logger.warning("Cannot obtain feature_slice for a feature without attached slice.")
            return None

        return Slice(
            seq_region_name=slice.seq_region_name,
            seq_region_length=slice.seq_region_length,
            coord_system=slice.coord_system,
            start=self.seq_region_start,
            end=self.seq_region_end,
            strand=self.seq_region_strand,
        )
