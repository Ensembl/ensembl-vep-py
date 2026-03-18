from attrs import define

from ensembl.vep.models.slice import Slice


@define(kw_only=True, auto_attribs=True)
class CircularSlice(Slice):
    """
    Represents a circular genomic region on a sequence (chromosome, contig, patch, etc.).
    """

    def is_circular(self) -> bool:
        """Indicates whether the slice is circular.

        Returns:
            bool: True, indicating the slice is circular.
        """
        return True