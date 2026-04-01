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

from ensembl.vep.common import Locus, trim_sequences, shift_sequence, invalidates
from ensembl.vep.models import Feature
if TYPE_CHECKING:
    from ensembl.vep.models.overlap import OverlapConsequence, TranscriptVariationAllele

logger = logging.getLogger(__name__)
MAX_SEARCH_AREA = 1000


@define(kw_only=True, auto_attribs=True)
class VariationAllele(Feature):
    start: int = field(on_setattr=invalidates("__original_locus", "__minimised_locus", "__shifted_locus", "__minimised", "__shifted"))
    end: int = field(on_setattr=invalidates("__original_locus", "__minimised_locus", "__shifted_locus", "__minimised", "__shifted"))
    reference_allele: str | None = field(default=None, on_setattr=invalidates("__original_locus", "__minimised_locus", "__shifted_locus", "__minimised", "__shifted"))
    alternate_allele: str | None = field(default=None, on_setattr=invalidates("__original_locus", "__minimised_locus", "__shifted_locus", "__minimised", "__shifted"))
    _lock: Lock = field(factory=Lock, init=False, repr=False)
    transcript_variation_alleles: list[TranscriptVariationAllele] = field(default=[], init=False)
    # regulatory_variation_alleles: list[RegulatoryVariationAllele] = field(default=[], init=False)
    # motif_variation_alleles: list[MotifVariationAllele] = field(default=[], init=False)
    # intergenic_variation_alleles: list[IntergenicVariationAllele] = field(default=[], init=False)

    @property
    def original_locus(self) -> Locus | None:
        """Public accessor for __original_locus.
        
        Returns:
            Locus | None: The original locus if available, None otherwise.
        """
        return self.__original_locus
    
    @cached_property
    def __original_locus(self) -> Locus | None:
        """Returns the original locus of the variation allele before any minimisation or shifting.
        
        Returns:
            Locus | None: The original locus if available, None otherwise.
        """
        return None
    
    @property
    def minimised_locus(self) -> Locus | None:
        """Public accessor for __minimised_locus.
        
        Returns:
            Locus | None: The minimised locus if available, None otherwise.
        """
        return self.__minimised_locus
    
    @cached_property
    def __minimised_locus(self) -> Locus | None:
        """Returns the locus of the variation allele after minimisation.
        
        Returns:
            Locus | None: The minimised locus if available, None otherwise.
        """
        return None
    
    @property
    def shifted_locus(self) -> Locus | None:
        """Public accessor for __shifted_locus.
        
        Returns:
            Locus | None: The shifted locus if available, None otherwise.
        """
        return self.__shifted_locus
    
    @cached_property
    def __shifted_locus(self) -> Locus | None:
        """Returns the locus of the variation allele after shifting.
        
        Returns:
            Locus | None: The shifted locus if available, None otherwise.
        """
        return None

    @property
    def minimised(self) -> bool:
        """Public accessor for __minimised.
        
        Returns:
            bool: True if the variation allele has been minimised, False otherwise.
        """
        return self.__minimised
    
    @cached_property
    def __minimised(self) -> bool:
        """Returns whether the variation allele has been minimised.
        
        Returns:
            bool: True if the variation allele has been minimised, False otherwise.
        """       
        return False

    @property
    def shifted(self) -> bool:
        """Public accessor for __shifted.
        
        Returns:
            bool: True if the variation allele has been shifted, False otherwise.
        """
        return self.__shifted
    
    @cached_property
    def __shifted(self) -> bool:
        """Returns whether the variation allele has been shifted.
        
        Returns:
            bool: True if the variation allele has been shifted, False otherwise.
        """
        return False

    def minimise(self, force=False) -> None:
        """Minimises the variation allele by computing the minimal representation of the variant.
        
        This method also stores the original locus information before minimisation.

        Args:
            force (bool): If True, force minimisation even if already minimised.
        
        Raises:
            ValueError: If reference or alternate allele is None.
        """
        if self.minimised and not force:
            logger.debug("Variation allele is already minimised. Use force=True to re-minimise.")
            return
        
        original_locus = self.__original_locus or Locus(
            start=self.start,
            end=self.end,
            allele_string=f"{self.reference_allele}/{self.alternate_allele}"
        )
        shifted = self.shifted
        shifted_locus = self.shifted_locus
        
        (self.reference_allele, self.alternate_allele, self.start, self.end, changed) = trim_sequences(
            ref=self.reference_allele, 
            alt=self.alternate_allele, 
            start=self.start, 
            end=self.end
        )
        
        if self.__original_locus is None:
            self.__original_locus = original_locus
        self.__minimised_locus = Locus(
            start=self.start,
            end=self.end,
            allele_string=f"{self.reference_allele}/{self.alternate_allele}"
        )
        self.__minimised = True
        if not changed:
            self.__shifted = shifted
            self.__shifted_locus = shifted_locus

    def shift(self, force=False) -> None:
        """Shifts the variation allele to the leftmost position in the sequence.
        
        This method also stores the original locus information before shifting.

        Args:
            force (bool): If True, force shifting even if already shifted.
        
        Raises:
            ValueError: If reference or alternate allele is None.
        """
        if self.shifted and not force:
            logger.debug("Variation allele is already shifted. Use force=True to re-shift.")
            return

        if self.reference_allele is None or self.alternate_allele is None:
            raise ValueError("Missing reference or alternate allele for shifting")

        original_locus = self.__original_locus or Locus(
            start=self.start,
            end=self.end,
            allele_string=f"{self.reference_allele}/{self.alternate_allele}"
        )
        minimised = self.minimised
        minimised_locus = self.minimised_locus

        search_start = max(0, self.start - MAX_SEARCH_AREA - 1)
        search_end = min(self.slice.end, self.end + MAX_SEARCH_AREA)
        pre_seq = self.slice.seq[search_start:(self.start - 1)]
        post_seq = self.slice.seq[self.end:search_end]
        if len(self.reference_allele) < len(self.alternate_allele):
            # for insertions, we check alternate allele
            seq_to_check = self.alternate_allele
            using_ref = False
        else:
            # for SNV and deletions, we check reference allele
            seq_to_check = self.reference_allele
            using_ref = True

        (_, shifted_seq, self.start, self.end, changed) = shift_sequence(
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

        if self.__original_locus is None:
            self.__original_locus = original_locus
        self.__shifted_locus = Locus(
            start=self.start,
            end=self.end,
            allele_string=f"{self.reference_allele}/{self.alternate_allele}"
        )
        self.__shifted = True
        if not changed:
            self.__minimised = minimised
            self.__minimised_locus = minimised_locus