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

from __future__ import annotations

import logging
from attrs import define, field
from typing import TYPE_CHECKING
from functools import cached_property
from threading import Lock

from ensembl.vep.models import Feature
from ensembl.vep.common import ExonPhase, invalidates

logger = logging.getLogger(__name__)


@define(kw_only=True, auto_attribs=True)
class Exon(Feature):
    """
    Represents an exon.

    Attributes:
        stable_id:  Ensembl stable exon identifier (e.g. "ENSE00003666667").
        version:    Version of the stable identifier.
    """
    stable_id: str
    version: int | None = field(default=None, init=False)
    phase: ExonPhase | None = field(default=None, init=False, on_setattr=invalidates("__frame"))
    end_phase: ExonPhase | None = field(default=None, init=False)
    _lock: Lock = field(factory=Lock, init=False, repr=False)

    @property
    def frame(self) -> int | None:
        """Public accessor for __frame.

        Returns:
            int | None: The frame of the exon, or None if not applicable.
        """
        return self.__frame

    @cached_property
    def __frame(self) -> int | None:
        """
        Returns the frame of the exon, which is the number of bases from the start of the exon to the next codon boundary. This is 0-based and can be 0, 1, or 2. If the phase is None, then the frame is also None.

        Returns:
            int | None: The frame of the exon.
        """
        if self.phase is None or self.phase == ExonPhase.UNKNOWN:
            logger.warning(f"Cannot calculate frame for exon {self.stable_id} because phase is None or UNKNOWN")
            return None
        
        if self.phase == ExonPhase.ZERO:
            return self.start%3
        elif self.phase == ExonPhase.ONE:
            return (self.start + 2) % 3
        elif self.phase == ExonPhase.TWO:
            return (self.start + 1) % 3
        else:
            logger.warning(f"Cannot calculate frame for exon {self.stable_id} because phase is invalid: {self.phase}")
            return None