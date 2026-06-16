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
"""CoordSystem class."""

from __future__ import annotations

from attrs import define, field

from ensembl.vep.common import CoordSystemLevel


@define(kw_only=True, auto_attribs=True)

__all__ = ["CoordSystem"]


class CoordSystem:
    """
    Represents a coordinate system in the Ensembl database.

    Attributes:
        name (str): Name of the coordinate system.
        level (CoordSystemLevel): Level of the coordinate system (e.g. chromosome, contig, toplevel).
        default (bool): Whether this is the default coordinate system.
        alias_to (str | None): Name of another coordinate system this is an alias to, if applicable.
        rank (int | None): Rank of the coordinate system.
    """

    name: str
    level: CoordSystemLevel
    default: bool = field(default=False)
    alias_to: str | None = field(default=None)
    rank: int | None = field(default=None)

    @property
    def is_sequence_level(self) -> bool:
        """
        Returns True if the coordinate system is sequence-level (e.g. contig), False otherwise.

        Returns:
            bool: True if sequence-level, False otherwise.
        """
        return self.level == CoordSystemLevel.CONTIG

    @property
    def is_chromosome_level(self) -> bool:
        """
        Returns True if the coordinate system is chromosome-level, False otherwise.

        Returns:
            bool: True if chromosome-level, False otherwise.
        """
        return self.level == CoordSystemLevel.CHROMOSOME

    @property
    def is_toplevel(self) -> bool:
        """
        Returns True if the coordinate system is toplevel, False otherwise.

        Returns:
            bool: True if toplevel, False otherwise.
        """
        return self.level == CoordSystemLevel.TOPLEVEL

    def __eq__(self, other: object) -> bool:
        """Compare two CoordSystem objects for equality.

        Args:
            other (object): Another object to compare against.

        Returns:
            bool: True if equal, False otherwise.
        """
        if type(self) is not type(other):
            return False

        return self.name == other.name and self.level == other.level and self.rank == other.rank
