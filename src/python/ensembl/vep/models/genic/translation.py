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
"""Translation class."""

from __future__ import annotations

from typing import TYPE_CHECKING
from functools import cached_property
from threading import RLock
from attrs import define, field

from ensembl.vep.common import invalidates

if TYPE_CHECKING:
    from ensembl.vep.models.genic import Transcript


@define(kw_only=True, auto_attribs=True)
class Translation:
    """
    Represents a translation.

    Attributes:
        start (int):
            The start position of the translation on the transcript,
            relative to the start of the transcript.
        end (int):
            The end position of the translation on the transcript,
            relative to the start of the transcript.
        stable_id (str):
            Ensembl stable transcript identifier (e.g. "ENST00000366667").
        version (int | None):
            Version of the stable identifier.
        seq (str | None):
            The amino acid sequence of the translation.
        transcript (Transcript | None):
            The transcript to which this translation belongs.
        exon_coding_start (int | None):
            The position of the first coding base in the start exon, relative to the start of
            the exon (1-based).
        exon_coding_end (int | None):
            The position of the last coding base in the end exon, relative to the start of the
            exon (1-based).
    """

    start: int = field(on_setattr=invalidates("_length"))
    end: int = field(on_setattr=invalidates("_length"))
    stable_id: str
    version: int | None = field(default=None, init=False)
    seq: str | None = field(default=None, init=False)
    transcript: Transcript | None = field(default=None, init=False)
    exon_coding_start: int | None = field(default=None)
    exon_coding_end: int | None = field(default=None)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    @property
    def length(self) -> int:
        """Public accessor for _length.

        Returns:
            int: The length of the translation.
        """
        return self._length

    @cached_property
    def _length(self) -> int:
        """
        Returns the length of the feature.

        Returns:
            int: The length of the feature.
        """
        return self.end - self.start + 1
