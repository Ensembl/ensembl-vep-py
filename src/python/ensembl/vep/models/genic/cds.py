from __future__ import annotations

from attrs import define, field
from typing import TYPE_CHECKING
from functools import cached_property
from threading import Lock

from ensembl.vep.models import Feature
from ensembl.vep.common import invalidates
if TYPE_CHECKING:
    from ensembl.vep.models.genic import Transcript, Translation

@define(kw_only=True, auto_attribs=True)
class CDS(Feature):
    """
    Represents a CDS.

    Attributes:
        transcript: The transcript to which this CDS belongs.
        translation: The translation associated with this CDS.
    """
    seq: str | None = field(default=None, on_setattr=invalidates("__neucleotide_length", "__protein_length"))
    transcript: Transcript | None = field(default=None, init=False)
    translation: Translation | None = field(default=None, init=False)
    _lock: Lock = field(factory=Lock, init=False, repr=False)

    @property
    def neucleotide_length(self) -> int | None:
        """Public accessor for __neucleotide_length.

        Returns:
            int | None: Length of the CDS in nucleotides, or None if not applicable.
        """
        return self.__neucleotide_length
    
    @cached_property
    def __neucleotide_length(self) -> int | None:
        """Returns the length of the CDS in nucleotides.

        Returns:
            int | None: Length of the CDS in nucleotides, or None if not applicable.
        """
        with self._lock:
            return len(self.seq) if self.seq is not None else None
        
    @property
    def protein_length(self) -> int | None:
        """Public accessor for __protein_length.

        Returns:
            int | None: Length of the CDS in amino acids, or None if not applicable.
        """
        return self.__protein_length
    
    @cached_property
    def __protein_length(self) -> int | None:
        """Returns the length of the CDS in amino acids.

        Returns:
            int | None: Length of the CDS in amino acids, or None if not applicable.
        """
        with self._lock:
            neucleotide_length = self.__neucleotide_length
            return (neucleotide_length // 3) if neucleotide_length is not None else None
        
    def __eq__(self, other: CDS) -> bool:
        """Compare two CDS using various criteria.
        
        Args:
            other: Another CDS object to compare against.

        Returns:
            bool: True if equal, False otherwise.
        """
        if not isinstance(other, CDS):
            return NotImplemented
        
        return (
            self.start == other.start and
            self.end == other.end and
            self.strand == other.strand and
            self.seq == other.seq and
            self.transcript == other.transcript and
            self.translation == other.translation
        )