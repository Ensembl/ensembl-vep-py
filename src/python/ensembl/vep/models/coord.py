from __future__ import annotations

from attrs import define, field
from functools import cached_property
from threading import Lock

from ensembl.vep.common import CoordSystemLevel, CachedAttributeGuardMixin


@define(kw_only=True, auto_attribs=True)
class CoordSystem(CachedAttributeGuardMixin):
    """
    Represents a coordinate system in the Ensembl database.

    Args:
        name (str): Name of the coordinate system.
        level (CoordSystemLevel): Level of the coordinate system (e.g. chromosome, contig, toplevel).
        default (bool): Whether this is the default coordinate system.
        alias_to (str | None): Name of another coordinate system this is an alias to, if applicable.
        rank (int | None): Rank of the coordinate system.

    Attributes:
        name (str): Name of the coordinate system.
        level (CoordSystemLevel): Level of the coordinate system (e.g. chromosome, contig, toplevel).
        default (bool): Whether this is the default coordinate system.
        alias_to (str | None): Name of another coordinate system this is an alias to, if applicable.
        rank (int | None): Rank of the coordinate system.
    """
    name: str
    level: CoordSystemLevel
    default: bool = False
    alias_to: str | None = None
    rank: int | None = None
    _lock: Lock = field(factory=Lock, init=False, repr=False)
    _DERIVED = frozenset({
        "is_sequence_level",
        "is_chromosome_level", 
        "is_toplevel"
    })
    _INVALIDATES = {
        "level":[
            "is_sequence_level",
            "is_chromosome_level", 
            "is_toplevel"
        ]
    }

    @cached_property
    def is_sequence_level(self) -> bool:
        """
        Returns True if the coordinate system is sequence-level (e.g. contig), False otherwise.

        Returns:
            bool: True if sequence-level, False otherwise.
        """
        with self._lock:
            return self.level == CoordSystemLevel.CONTIG
    
    @cached_property
    def is_chromosome_level(self) -> bool:
        """
        Returns True if the coordinate system is chromosome-level, False otherwise.

        Returns:
            bool: True if chromosome-level, False otherwise.
        """
        with self._lock:
            return self.level == CoordSystemLevel.CHROMOSOME
    
    @cached_property
    def is_toplevel(self) -> bool:
        """
        Returns True if the coordinate system is toplevel, False otherwise.

        Returns:
            bool: True if toplevel, False otherwise.
        """
        with self._lock:
            return self.level == CoordSystemLevel.TOPLEVEL
        
    def __eq__(self, other: "CoordSystem") -> bool:
        """Compare two CoordSystem objects for equality.
        
        Args:
            other (CoordSystem): Another CoordSystem to compare against.
        
        Returns:
            bool: True if equal, False otherwise.
        """
        if type(self) is not type(other):
            return False

        return (
            self.name == other.name
            and self.level == other.level
            and self.rank == other.rank
        )