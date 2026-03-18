from __future__ import annotations

from attrs import define, field
from functools import cached_property

from ensembl.vep.common import VariationLocus, trim_sequences, shift_sequence
from ensembl.vep.models import Feature
# from ensembl.vep.models.overlap import OverlapConsequence, TranscriptVariation, RegulatoryVariation, MotifVariation, IntergenicVariation

MAX_SEARCH_AREA = 1000


@define(kw_only=True, auto_attribs=True)
class VariationAllele(Feature):
    reference_allele: str | None = None
    alternate_allele: str | None = None
    original_locus: VariationLocus | None = field(default=None, init=False)
    minimised_locus: VariationLocus | None = field(default=None, init=False)
    shifted_locus: VariationLocus | None = field(default=None, init=False)
    _minimised: bool = field(default=False, init=False, repr=False)
    _shifted: bool = field(default=False, init=False, repr=False)
    # overlap_consequences: list[OverlapConsequence] | None = None
    # _transcript_variations: list[TranscriptVariation] | None = field(default=None, init=False, repr=False)
    # _regulatory_variations: list[RegulatoryVariation] | None = field(default=None, init=False, repr=False)
    # _motif_variations: list[MotifVariation] | None = field(default=None, init=False, repr=False)
    # _intergenic_variation: IntergenicVariation | None = field(default=None, init=False, repr=False)
    _DERIVED = frozenset({
        "consequence_types",
        "minimised",
        "shifted"
    })
    _INVALIDATES = {
        "overlap_consequences":["consequence_types"],
    }

    @property
    def minimised(self) -> bool:
        """Returns True if the variation allele has been minimised, False otherwise."""
        return self._minimised
    
    @property
    def shifted(self) -> bool:
        """Returns True if the variation allele has been shifted, False otherwise."""
        return self._shifted

    @cached_property
    def consequence_types(self) -> list[str] | None:
        """Returns the list of consequence types for the variation allele.
        
        Returns:
            list[str] | None: List of consequence types or None if not available."""
        for overlap_consequence in self.overlap_consequences or []:
            if overlap_consequence.SO_term not in self.consequence_types:
                self.consequence_types.append(overlap_consequence.SO_term)

        return self.consequence_types

    def minimise(self) -> None:
        """Minimises the variation allele by computing the minimal representation of the variant.
        
        This method also stores the original locus information before minimisation.
        
        Raises:
            ValueError: If reference or alternate allele is None.
        """
        if self.original_locus is None:
            self.original_locus = VariationLocus(
                start=self.start,
                end=self.end,
                allele_string=f"{self.reference_allele}/{self.alternate_allele}"
            )
        
        (self.reference_allele, self.alternate_allele, self.start, self.end, _) = trim_sequences(
            ref=self.reference_allele, 
            alt=self.alternate_allele, 
            start=self.start, 
            end=self.end
        )
        
        self.minimised_locus = VariationLocus(
            start=self.start,
            end=self.end,
            allele_string=f"{self.reference_allele}/{self.alternate_allele}"
        )

        if 'minimised' in self.__slots__:
            delattr(self, 'minimised')
        self._minimised = True

    def shift(self) -> None:
        """Shifts the variation allele to the leftmost position in the sequence.
        
        This method also stores the original locus information before shifting.
        
        Raises:
            ValueError: If reference or alternate allele is None.
        """
        if self.reference_allele is None or self.alternate_allele is None:
            raise ValueError("Missing reference or alternate allele for shifting")

        if self.original_locus is None:
            self.original_locus = VariationLocus(
                start=self.start,
                end=self.end,
                allele_string=f"{self.reference_allele}/{self.alternate_allele}"
            )

        if len(self.reference_allele) > len(self.alternate_allele):
            seq_to_check = self.reference_allele
            using_ref = True
        else:
            seq_to_check = self.alternate_allele
            using_ref = False
        search_start = max(0, self.start - MAX_SEARCH_AREA - 1)
        search_end = min(self.slice.end, self.end + MAX_SEARCH_AREA)
        pre_seq = self.slice.seq[search_start:(self.start - 1)]
        post_seq = self.slice.seq[self.end:search_end]
        (_, shifted_seq, self.start, self.end) = shift_sequence(
            seq=seq_to_check, 
            post_seq=post_seq, 
            pre_seq=pre_seq, 
            start=self.start, 
            end=self.end
        )
        if using_ref:
            self.reference_allele = shifted_seq
        else:
            self.alternate_allele = shifted_seq

        self.shifted_locus = VariationLocus(
            start=self.start,
            end=self.end,
            allele_string=f"{self.reference_allele}/{self.alternate_allele}"
        )

        if 'shifted' in self.__slots__:
            delattr(self, 'shifted')
        self._shifted = True