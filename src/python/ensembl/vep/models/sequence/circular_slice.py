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
"""CircularSlice class."""

from __future__ import annotations

from functools import cached_property
import logging

from attrs import define

from ensembl.vep.common import Strand, reverse_complement
from ensembl.vep.models.sequence.slice import Slice

logger = logging.getLogger(__name__)


@define(kw_only=True, auto_attribs=True)
class CircularSlice(Slice):
    """Represents a circular genomic region."""

    @cached_property
    def _seq(self) -> str | None:
        """Resolve the sequence for a circular slice.

        Returns:
            str | None: The resolved sequence, or None if unavailable.
        """

        with self._lock:

            if self.sequence_source is None:
                return None

            seq_region_length = self.seq_region_length
            if seq_region_length is None:
                logger.warning(
                    "Cannot resolve circular slice on %s without seq_region_length.",
                    self.seq_region_name,
                )
                return None

            spans = self._sequence_spans(seq_region_length)
            if not spans:
                return None

            sequence = ""
            for start, end in spans:
                piece = self.sequence_source.get_sequence_by_start_end_strand(
                    self.seq_region_name,
                    start,
                    end,
                    strand=Strand.FORWARD,
                )
                if piece is None:
                    piece = "N" * (end - start + 1)
                sequence += piece

            if self.strand == Strand.REVERSE:
                return reverse_complement(sequence)

            return sequence

    def _sequence_spans(self, seq_region_length: int) -> list[tuple[int, int]]:
        """Return linear spans that cover the circular slice.

        Args:
            seq_region_length (int): The length of the circular sequence region.

        Returns:
            list[tuple[int, int]]: Linear spans used to assemble the slice sequence.
        """

        if self.start >= 1 and self.end >= self.start and self.end <= seq_region_length:
            return [(self.start, self.end)]

        normalized_start = ((self.start - 1) % seq_region_length) + 1
        normalized_end = ((self.end - 1) % seq_region_length) + 1

        if normalized_start <= normalized_end and self.start <= seq_region_length and self.end <= seq_region_length:
            return [(normalized_start, normalized_end)]

        if normalized_start <= seq_region_length:
            return [(normalized_start, seq_region_length), (1, normalized_end)]

        return [(normalized_start, seq_region_length), (1, normalized_end)]

    def is_circular(self) -> bool:
        """Indicates whether the slice is circular.

        Returns:
            bool: True, indicating the slice is circular.
        """
        return True

    # TODO: We needed to override __eq__ of attrs regeneration which compares internal field like _lock; 
    # is there a generic way to do this?
    def __eq__(self, other: object) -> bool:
        """Compare two CircularSlice objects and define if they are identical."""

        return Slice.__eq__(self, other)
