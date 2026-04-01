from __future__ import annotations

from attrs import define, field
from typing import TYPE_CHECKING
from functools import cached_property
from threading import Lock

from ensembl.vep.common import invalidates

if TYPE_CHECKING:
    from ensembl.vep.models.genic import Transcript


@define(kw_only=True, auto_attribs=True)
class Translation():
    """
    Represents a transcript.

    Attributes:
        stable_id:  Ensembl stable transcript identifier (e.g. "ENST00000366667").
        version:    Version of the stable identifier.
        exon_coding_start: The position of the first coding base in the start exon, relative to the start of the exon. This is 1-based and defaults to 1.
        exon_coding_end: The position of the last coding base in the end exon, relative to the start of the exon. This is 1-based and defaults to length of the end exon.
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
        Returns the position of the first coding base in the start exon, relative to the start of the exon. This is 1-based and defaults to 1.

        Returns:
            int: The position of the first coding base in the start exon, relative to the start of the exon.
        """
        return 1

    @property
    def exon_coding_end(self) -> int | None:
        """
        Returns the position of the last coding base in the end exon, relative to the start of the exon. This is 1-based and defaults to length of the end exon.

        Returns:
            int | None: The position of the last coding base in the end exon, relative to the start of the exon.
        """
        return self.__exon_coding_end

    @cached_property
    def __exon_coding_end(self) -> int | None:
        """
        Returns the position of the last coding base in the end exon, relative to the start of the exon. This is 1-based and defaults to the length of the end exon.

        Returns:
            int | None: The position of the last coding base in the end exon, relative to the start of the exon.
        """
        if (
                self.transcript is None
                or self.transcript.spliced_exons is None
                or len(self.transcript.spliced_exons) == 0
        ):
            return None
        
        return len(self.transcript.spliced_exons[-1].seq)