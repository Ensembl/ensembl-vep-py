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
"""Exon class."""

from __future__ import annotations

import logging
from functools import cached_property
from threading import RLock
from attrs import define, field

from ensembl.vep.models import Feature
from ensembl.vep.common import ExonPhase, invalidates

logger = logging.getLogger(__name__)


@define(kw_only=True, auto_attribs=True)

__all__ = ["Exon"]


class Exon(Feature):
    """
    Represents an exon.

    Attributes:
        stable_id (str):
            Ensembl stable exon identifier (e.g. "ENSE00003666667").
        version (int | None):
            Version of the stable identifier.
        phase (ExonPhase | None):
            The phase of the exon, which is the number of bases of the first codon
            that are on the previous exon/intron. The value can be 0, 1, or 2.
        end_phase (ExonPhase | None):
            Number of bases of the last codon if it is incomplete. The value can be 0, 1, or 2.
    """

    stable_id: str
    version: int | None = field(default=None, init=False)
    phase: ExonPhase | None = field(default=None, init=False, on_setattr=invalidates("_frame"))
    end_phase: ExonPhase | None = field(default=None, init=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    @property
    def frame(self) -> int | None:
        """Public accessor for _frame.

        Returns:
            int | None: The frame of the exon, or None if not applicable.
        """
        return self._frame

    @cached_property
    def _frame(self) -> int | None:
        """
        Returns the frame of the exon.

        Frame tells where the exon starts in relation to the reading frame as given by the phase in GFF.

        Returns:
            int | None: The frame of the exon.
        """
        if self.phase is None or self.phase == ExonPhase.UNKNOWN:
            logger.warning(
                f"Cannot calculate frame for exon {self.stable_id} because phase is None or UNKNOWN"
            )
            return None

        # TBD: why we use module of start position? instead would it be simple convert - 0->0, 1->2, 2->1?
        if self.phase == ExonPhase.ZERO:
            return self.start % 3
        if self.phase == ExonPhase.ONE:
            return (self.start + 2) % 3
        if self.phase == ExonPhase.TWO:
            return (self.start + 1) % 3

        logger.warning(
            f"Cannot calculate frame for exon {self.stable_id} because phase is invalid: {self.phase}"
        )
        return None
