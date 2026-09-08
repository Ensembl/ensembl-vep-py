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
from collections.abc import Mapping
from functools import cached_property
from threading import RLock
from attrs import define, field

from ensembl.vep.models.sequence import Slice, CoordSystem
from ensembl.vep.common import Strand, UnsafeSetterMixin, invalidates, reverse_complement

logger = logging.getLogger(__name__)


@define(kw_only=True, auto_attribs=True)
class Feature(UnsafeSetterMixin):
    """
    Represents a genomic feature in a specific locus.

    Attributes:
        start: Start position relative to slice (1-based, inclusive).
        end: End position relative to slice (1-based, inclusive).
        strand (Strand): Feature orientation relative to slice.
        slice (Slice | None): Slice object describing the feature context.
        seq (str | None): Associated biological sequence.
    """

    start: int = field(on_setattr=invalidates("_length", "_seq_region_start", "_seq"))
    end: int = field(on_setattr=invalidates("_length", "_seq_region_end", "_seq"))
    strand: Strand = field(on_setattr=invalidates("_seq_region_strand"))
    slice: Slice | None = field(
        default=None,
        on_setattr=invalidates("_seq_region_start", "_seq_region_end", "_seq_region_strand", "_seq"),
    )
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    @property
    def length(self) -> int:
        """
        Public accessor for the _length.

        Returns:
            int: The length of the feature.
        """
        return self._length

    @cached_property
    def _length(self) -> int:
        """
        Returns the length of the feature.

        Returns:
            int: The length of the feature.
        """
        with self._lock:
            return self.end - self.start + 1

    @property
    def seq_region_start(self) -> int | None:
        """
        Public accessor for the _seq_region_start.

        Returns:
            The start position in the sequence region,
            or None if the feature has no attached slice.
        """
        return self._seq_region_start

    @cached_property
    def _seq_region_start(self) -> int | None:
        """
        Returns the start position in the sequence region.

        Returns:
            The start position in the sequence region,
            or None if the feature has no attached slice.
        """
        with self._lock:
            if self.slice is None:
                logger.warning("Cannot obtain seq_region_start for a feature without attached slice.")
                return None
            return self.slice.start + self.start - 1

    @property
    def seq_region_end(self) -> int | None:
        """
        Public accessor for the _seq_region_end.

        Returns:
            int | None: The end position in the sequence region, or None if the feature has no attached slice.
        """
        return self._seq_region_end

    @cached_property
    def _seq_region_end(self) -> int | None:
        """
        Returns the end position in the sequence region.

        Returns:
            int | None: The end position in the sequence region, or None if the feature has no attached slice.
        """
        with self._lock:
            if self.slice is None:
                logger.warning("Cannot obtain seq_region_end for a feature without attached slice.")
                return None
            return self.slice.start + self.end - 1

    @property
    def seq_region_strand(self) -> Strand | None:
        """
        Public accessor for the _seq_region_strand.

        Returns:
            Strand | None: The strand in the sequence region, or None if the feature has no attached slice.
        """
        return self._seq_region_strand

    @cached_property
    def _seq_region_strand(self) -> Strand | None:
        """
        Returns the strand in the sequence region.

        Returns:
            Strand | None: The strand in the sequence region, or None if the feature has no attached slice.
        """
        with self._lock:
            if self.slice is None:
                logger.warning("Cannot obtain seq_region_strand for a feature without attached slice.")
                return None
            return Strand(self.slice.strand.value * self.strand.value)

    @property
    def seq(self) -> str | None:
        """
        Public accessor for the _seq.

        Returns:
            str | None: The associated biological sequence, or None if unavailable.
        """
        return self._seq

    @cached_property
    def _seq(self) -> str | None:
        """
        Returns the associated biological sequence.

        Here, ``start``/``end`` are always ascending order in the attached slice and ``strand`` here
        is the feature's orientation relative to the slice. Therefore when feature's strand disagrees
        with the slice's strand, the slice sequence must be reverse-complemented before being returned.

        Returns:
            str | None: The associated biological sequence, in this feature's own 5'->3'
                orientation, or None if unavailable.
        """
        with self._lock:
            if self.slice is None:
                logger.warning("Cannot obtain seq for a feature without attached slice.")
                return None

            slice_seq = self.slice.seq
            if slice_seq is None:
                return None

            feature_seq = slice_seq[self.start - 1 : self.end]
            if self.strand == Strand.REVERSE:
                return reverse_complement(feature_seq)
            return feature_seq

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
        raise NotImplementedError("Feature transformation is not implemented yet.")

    def transfer(self, _genomic_slice: Slice) -> Feature | None:
        """
        Transfers the feature to a new slice.

        Args:
            slice (Slice, optional): The slice to transfer to. If None, uses current slice.

        Returns:
            Feature | None: The transferred feature or None if transfer is not possible.
        """
        raise NotImplementedError("Feature transfer is not implemented yet.")

    def overlaps(self, other: Feature) -> bool | None:
        """
        Checks if this feature overlaps with another feature.

        Args:
            other (Feature): Another feature to check overlap with.

        Returns:
            bool | None: True if overlaps, False if not, None if slices are missing.
        """
        if self.slice is None or other.slice is None:
            logger.warning("Cannot compare two features is any of them is without attached slice.")
            return None

        if self.slice.seq_region_name != other.slice.seq_region_name:
            logger.warning(
                "Cannot compare two features on different sequence regions."
                f" {self.__class__.__name__} is on {self.slice.seq_region_name}"
                f" and {other.__class__.__name__} is on {other.slice.seq_region_name}."
            )
            return False

        if (
            self.seq_region_start is not None
            and other.seq_region_start is not None
            and self.seq_region_end is not None
            and other.seq_region_end is not None
        ):
            return (
                self.seq_region_start <= other.seq_region_end
                and self.seq_region_end >= other.seq_region_start
            )

        if self.slice != other.slice:
            return self.start <= other.end and self.end >= other.start

        return None

    def feature_slice(self) -> Slice | None:
        """
        Returns a slice that covers the area of this feature.

        The feature start will be at 1 on the returned slice, and it will have the length of this feature.
        If the feature has no attached slice, a warning will be logged and None will be returned.

        Returns:
            A Slice object covering the feature area,
            or None if the feature has no attached slice.
        """
        genomic_slice = self.slice

        if genomic_slice is None:
            logger.warning("Cannot obtain feature_slice for a feature without attached slice.")
            return None

        if self.seq_region_start is None or self.seq_region_end is None or self.seq_region_strand is None:
            logger.warning(
                (
                    "Cannot obtain feature_slice for a feature with undefined"
                    " seq_region_start, seq_region_end, or seq_region_strand."
                )
            )
            return None

        return Slice(
            seq_region_name=genomic_slice.seq_region_name,
            coord_system=genomic_slice.coord_system,
            start=self.seq_region_start,
            end=self.seq_region_end,
            strand=self.seq_region_strand,
            sequence_source=genomic_slice.sequence_source,
        )
