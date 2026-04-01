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

from typing import Any
from attrs import define, field

from ensembl.vep.common import Strand


@define(kw_only=True, auto_attribs=True)
class Slice:
    """
    Represents a genomic region on a sequence (chromosome, contig, patch, etc.).

    Args:
        coord_system (str): Coordinate system name.
        seq_region_name (str): Sequence region name.
        seq_region_length (str): Sequence region length.
        start (int): Start position relative to sequence region.
        end (int): End position relative to sequence region.
        strand (Strand): Slice orientation relative to sequence region.
        seq (Any | None): Associated biological sequence.

    Attributes:
        coord_system (str): Coordinate system name.
        seq_region_name (str): Sequence region name.
        seq_region_length (str): Sequence region length.
        start (int): Start position relative to sequence region.
        end (int): End position relative to sequence region.
        strand (Strand): Slice orientation relative to sequence region.
        seq (Any | None): Associated biological sequence.
    """

    coord_system: str
    seq_region_name: str
    start: int
    end: int
    strand: Strand
    seq_region_length: int | None = field(default=None)
    seq: Any | None = field(default=None)

    def __eq__(self, other: object) -> bool:
        """Compare two Slice using various criteria.

        Args:
            other (object): Another object to compare against.

        Returns:
            bool: True if equal, False otherwise.
        """
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
        """Indicates whether the slice is circular.

        By default, slices are not circular, but this method can be overridden in subclasses to provide specific behavior.

        Returns:
            bool: False by default, indicating the slice is not circular.
        """
        return False
