from __future__ import annotations

from attrs import define, field
from typing import Any

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

    def __eq__(self, other: "Slice") -> bool:
        """Compare two Slice using various criteria.
        
        Args:
            other (Slice): Another slice to compare against.
        
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