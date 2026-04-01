from __future__ import annotations

from attr import define, field
from typing import TYPE_CHECKING
from functools import cached_property
from threading import Lock

from ensembl.vep.common import invalidates
from ensembl.vep.models.overlap import FeatureOverlap

if TYPE_CHECKING:
    from ensembl.vep.models.genic import Transcript
    from ensembl.vep.models.variation import VariationAllele
    from ensembl.vep.common import Locus


@define(kw_only=True, auto_attribs=True)
class TranscriptVariationAllele(FeatureOverlap):
    """
    Represents the overlap between a transcript and a variation allele.
    
    Attributes:
        feature: The transcript that is getting overlapped.
        overlapping_feature: The variation allele that is overlapping.
        is_reference: Whether the overlapping variation allele is the reference allele.
    """
    feature: Transcript
    overlapping_feature: VariationAllele
    is_reference: bool = field(default=False, init=False)
    cdna_locus: Locus = field(init=False, on_setattr=invalidates("hgvs_transcript"))
    cds_locus: Locus | None = field(default=None, init=False, on_setattr=invalidates("hgvs_transcript", "hgvs_protein"))
    protein_locus: Locus | None = field(default=None, init=False, on_setattr=invalidates("hgvs_protein"))
    distance_to_transcript: int | None = field(default=None, init=False)
    exon_number: int | None = field(default=None, init=False)
    intron_number: int | None = field(default=None, init=False)
    _lock: Lock = field(factory=Lock, init=False, repr=False)

    @property
    def shifted(self) -> bool:
        """
        Public accessor for __shifted.

        Returns:
            bool: True if the allele has been shifted, False otherwise.
        """
        return self.__shifted
    
    @cached_property
    def __shifted(self) -> bool:
        """
        Returns whether the allele has been shifted from its original position.

        Returns:
            bool: True if the allele has been shifted, False otherwise.
        """
        return False
    
    @property
    def shifted_cdna_locus(self) -> Locus | None:
        """
        Public accessor for __shifted_cdna_locus.

        Returns:
            Locus | None: The cDNA locus of the variation allele after shifting, or None if not applicable.
        """
        return self.__shifted_cdna_locus
    
    @cached_property
    def __shifted_cdna_locus(self) -> Locus | None:
        """
        Returns the cDNA locus of the variation allele after shifting, or None if not applicable.

        Returns:
            Locus | None: The cDNA locus of the variation allele after shifting, or None if not applicable.
        """
        return None
    
    @property
    def shifted_cds_locus(self) -> Locus | None:
        """
        Public accessor for __shifted_cds_locus.

        Returns:
            Locus | None: The CDS locus of the variation allele after shifting, or None if not applicable.
        """
        return self.__shifted_cds_locus
    
    @cached_property
    def __shifted_cds_locus(self) -> Locus | None:
        """
        Returns the CDS locus of the variation allele after shifting, or None if not applicable.

        Returns:
            Locus | None: The CDS locus of the variation allele after shifting, or None if not applicable.
        """
        return None
    
    @property
    def shifted_protein_locus(self) -> Locus | None:
        """
        Public accessor for __shifted_protein_locus.

        Returns:
            Locus | None: The protein locus of the variation allele after shifting, or None if not applicable.
        """
        return self.__shifted_protein_locus
    
    @cached_property
    def __shifted_protein_locus(self) -> Locus | None:
        """
        Returns the protein locus of the variation allele after shifting, or None if not applicable.

        Returns:
            Locus | None: The protein locus of the variation allele after shifting, or None if not applicable.
        """
        return self.__shifted_protein_locus
    
    @property
    def original_cdna_locus(self) -> Locus | None:
        """Public accessor for __original_cdna_locus.

        Returns:
            Locus | None: The original cDNA locus of the variation allele before shifting, or None if not applicable.
        """
        return self.__original_cdna_locus
    
    @cached_property
    def __original_cdna_locus(self) -> Locus | None:
        """Returns the original cDNA locus of the variation allele before shifting, or None if not applicable.

        Returns:
            Locus | None: The original cDNA locus of the variation allele before shifting, or None if not applicable.
        """
        return None
    
    @property
    def original_cds_locus(self) -> Locus | None:
        """Public accessor for __original_cds_locus.

        Returns:
            Locus | None: The original CDS locus of the variation allele before shifting, or None if not applicable.
        """
        return self.__original_cds_locus
    
    @cached_property
    def __original_cds_locus(self) -> Locus | None:
        """Returns the original CDS locus of the variation allele before shifting, or None if not applicable.

        Returns:
            Locus | None: The original CDS locus of the variation allele before shifting, or None if not applicable.
        """
        return None
    
    @property
    def original_protein_locus(self) -> Locus | None:
        """Public accessor for __original_protein_locus.

        Returns:
            Locus | None: The original protein locus of the variation allele before shifting, or None if not applicable.
        """
        return self.__original_protein_locus
    
    @cached_property
    def __original_protein_locus(self) -> Locus | None:
        """Returns the original protein locus of the variation allele before shifting, or None if not applicable.

        Returns:
            Locus | None: The original protein locus of the variation allele before shifting, or None if not applicable.
        """
        return None
    
    def shift(self) -> None:
        pass

    def hgvs_genomic(self) -> str | None:
        pass

    def hgvs_transcript(self) -> str | None:
        pass

    def hgvs_protein(self) -> str | None:
        pass