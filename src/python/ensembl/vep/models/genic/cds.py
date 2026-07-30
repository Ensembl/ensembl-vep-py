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
"""CDS class."""

from __future__ import annotations

from typing import TYPE_CHECKING
from functools import cached_property
from threading import RLock
from attrs import define, field

from ensembl.vep.models import Feature
from ensembl.vep.common import invalidates

if TYPE_CHECKING:
    from ensembl.vep.models.genic import Transcript, Translation


@define(kw_only=True, auto_attribs=True)

__all__ = ["CDS"]


class CDS(Feature):
    """
    Represents a CDS.

    Attributes:
        seq (str | None): The nucleotide sequence of the CDS.
        transcript (Transcript | None): The transcript to which this CDS belongs.
        translation (Translation | None): The translation associated with this CDS.
    """

    seq: str | None = field(default=None, on_setattr=invalidates("_nucleotide_length", "_protein_length"))
    transcript: Transcript | None = field(default=None, init=False)
    translation: Translation | None = field(default=None, init=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    @property
    def nucleotide_length(self) -> int | None:
        """Public accessor for _nucleotide_length.

        Returns:
            int | None: Length in nucleotides, or None if not applicable.
        """
        return self._nucleotide_length

    @cached_property
    def _nucleotide_length(self) -> int | None:
        """Returns the length of the CDS in nucleotides.

        Returns:
            int | None: Length in nucleotides, or None if not applicable.
        """
        with self._lock:
            return len(self.seq) if self.seq is not None else None

    @property
    def protein_length(self) -> int | None:
        """Public accessor for _protein_length.

        Returns:
            int | None: Length in amino acids, or None if not applicable.
        """
        return self._protein_length

    @cached_property
    def _protein_length(self) -> int | None:
        """Returns the length of the CDS in amino acids.

        Returns:
            int | None: Length in amino acids, or None if not applicable.
        """
        with self._lock:
            nucleotide_length = self._nucleotide_length
            return (nucleotide_length // 3) if nucleotide_length is not None else None

    def __eq__(self, other: object) -> bool:
        """Compare two CDS using various criteria.

        Args:
            other (object): Another object to compare against.

        Returns:
            bool: True if equal, False otherwise.
        """
        if not isinstance(other, CDS):
            return NotImplemented

        return (
            self.start == other.start
            and self.end == other.end
            and self.strand == other.strand
            and self.seq == other.seq
            and self.transcript == other.transcript
            and self.translation == other.translation
        )
