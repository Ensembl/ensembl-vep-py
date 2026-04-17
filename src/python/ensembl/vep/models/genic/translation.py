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
from threading import Lock
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
    """

    start: int = field(on_setattr=invalidates("__length"))
    end: int = field(on_setattr=invalidates("__length"))
    stable_id: str
    version: int | None = field(default=None, init=False)
    seq: str | None = field(default=None, init=False)
    transcript: Transcript | None = field(default=None, init=False)
    _lock: Lock = field(factory=Lock, init=False, repr=False)

    @property
    def length(self) -> int:
        """Public accessor for __length.

        Returns:
            int: The length of the translation.
        """
        return self.__length

    @cached_property
    def __length(self) -> int:
        """
        Returns the length of the feature.

        Returns:
            int: The length of the feature.
        """
        return self.end - self.start + 1

    @property
    def exon_coding_start(self) -> int:
        """
        The position of the first coding base in the start exon, relative to the start of the exon.

        Returns:
            The position of the first coding base in the start exon,
            relative to the start of the exon.
        """
        return 1

    @property
    def exon_coding_end(self) -> int | None:
        """
        Returns the position of the last coding base in the end exon,
        relative to the start of the exon.

        Returns:
            The position of the last coding base in the end exon,
            relative to the start of the exon.
        """
        return self.__exon_coding_end

    @cached_property
    def __exon_coding_end(self) -> int | None:
        """
        Returns the position of the last coding base in the end exon,
        relative to the start of the exon.

        Returns:
            The position of the last coding base in the end exon,
            relative to the start of the exon.
        """
        if (
            self.transcript is None
            or self.transcript.spliced_exons is None
            or len(self.transcript.spliced_exons) == 0
            or self.transcript.spliced_exons[-1].seq is None
        ):
            return None

        return len(self.transcript.spliced_exons[-1].seq)
