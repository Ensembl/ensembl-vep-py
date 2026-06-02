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
"""VariationAllele class."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING
from functools import cached_property
from threading import RLock
from collections.abc import Sequence
from attrs import define, field

from ensembl.vep.common import Locus, trim_sequences, shift_sequence, invalidates, resets_attrs
from ensembl.vep.models import Feature

if TYPE_CHECKING:
    from ensembl.vep.models.overlap import OverlapConsequence, FeatureOverlap, TranscriptVariationAllele

logger = logging.getLogger(__name__)
MAX_SEARCH_AREA = 1000


@define(kw_only=True, auto_attribs=True)
class VariationAllele(Feature):
    """Represents a specific allele of a variation.

    Attributes:
        start (int):
            The start position in relative to slice.
        end (int):
            The end position in relative to slice.
        reference_allele (str | None):
            The reference allele.
        alternate_allele (str | None):
            The alternate allele.
        transcript_variation_alleles (list[TranscriptVariationAllele]):
            List of overlaps between this variation allele and transcripts.
        original_locus (Locus | None):
            The original locus of the variation allele before any minimisation or shifting.
        minimised_locus (Locus | None):
            The locus of the variation allele after minimisation.
        shifted_locus (Locus | None):
            The locus of the variation allele after shifting.
        minimised (bool):
            Whether the variation allele has been minimised.
        shifted (bool):
            Whether the variation allele has been shifted.
    """

    start: int = field(
        on_setattr=[resets_attrs(
            ("_shifted", False), 
            ("_minimised_locus", None), 
            ("_shifted_locus", None)
        )]
    )
    end: int = field(
        on_setattr=[resets_attrs(
            ("_shifted", False),
            ("_minimised_locus", None),
            ("_shifted_locus", None)
        )]
    )
    reference_allele: str | None = field(
        default=None,
        on_setattr=[resets_attrs(
            ("_minimised", False),
            ("_shifted", False),
            ("_minimised_locus", None),
            ("_shifted_locus", None)
        )]
    )
    alternate_allele: str | None = field(
        default=None,
        on_setattr=[resets_attrs(
            ("_minimised", False),
            ("_shifted", False),
            ("_minimised_locus", None),
            ("_shifted_locus", None)
        )]
    )
    _original_locus: Locus | None = field(default=None, init=False)
    """The original locus of the variation allele before any minimisation or shifting (if any occurred)."""
    _minimised: bool = field(default=False, init=False)
    """Whether the variation allele has been minimised."""
    _minimised_locus: Locus | None = field(default=None, init=False)
    _shifted: bool = field(default=False, init=False)
    """Whether the variation allele has been shifted."""
    _shifted_locus: Locus | None = field(default=None, init=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)
    transcript_variation_alleles: list[TranscriptVariationAllele] = field(
        factory=list, init=False, on_setattr=invalidates("_overlap_consequences")
    )
    # regulatory_variation_alleles: list[RegulatoryVariationAllele] = field(factory=list, init=False)
    # motif_variation_alleles: list[MotifVariationAllele] = field(factory=list, init=False)
    # intergenic_variation_alleles: list[IntergenicVariationAllele] = field(factory=list, init=False)

    @property
    def _locus(self) -> Locus:
        return Locus(
            start=self.start, end=self.end, allele_string=f"{self.reference_allele}/{self.alternate_allele}"
        )

    @property
    def original_locus(self) -> Locus | None:
        """Public accessor for _original_locus.

        Returns:
            Locus | None: The original locus if available, None otherwise.
        """
        return self._original_locus

    @property
    def minimised_locus(self) -> Locus | None:
        """Public accessor for _minimised_locus.

        Returns:
            Locus | None: The minimised locus if available, None otherwise.
        """
        return self._minimised_locus

    @property
    def shifted_locus(self) -> Locus | None:
        """Public accessor for _shifted_locus.

        Returns:
            Locus | None: The shifted locus if available, None otherwise.
        """
        return self._shifted_locus

    @property
    def minimised(self) -> bool:
        """Public accessor for _minimised.

        Returns:
            bool: True if the variation allele has been minimised, False otherwise.
        """
        return self._minimised

    @property
    def shifted(self) -> bool:
        """Public accessor for _shifted.

        Returns:
            bool: True if the variation allele has been shifted, False otherwise.
        """
        return self._shifted

    def minimise(self, force: bool = False) -> None:
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

        if self.reference_allele is None or self.alternate_allele is None:
            raise ValueError("Missing reference or alternate allele for minimisation")

        original_locus = self._original_locus or self._locus

        new_ref, new_alt, new_start, new_end, _changed = trim_sequences(
            ref=self.reference_allele, alt=self.alternate_allele, start=self.start, end=self.end
        )

        if _changed:
            self.reference_allele = new_ref
            self.alternate_allele = new_alt
            self.start = new_start
            self.end = new_end

        self._original_locus = self._original_locus or original_locus
        self._minimised_locus = self._locus
        self._minimised = True

    def shift(self, force: bool = False) -> None:
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

        if self.slice is None or self.slice.seq is None:
            raise ValueError("Cannot shift variation allele without attached slice sequence.")

        original_locus = self._original_locus or self._locus

        search_start = max(0, self.start - MAX_SEARCH_AREA - 1)
        search_end = min(self.slice.end, self.end + MAX_SEARCH_AREA)
        pre_seq = self.slice.seq[search_start : (self.start - 1)]
        post_seq = self.slice.seq[self.end : search_end]
        if len(self.reference_allele) < len(self.alternate_allele):
            # for insertions, we check alternate allele
            seq_to_check = self.alternate_allele
            using_ref = False
        else:
            # for SNV and deletions, we check reference allele
            seq_to_check = self.reference_allele
            using_ref = True

        _, shifted_seq, new_start, new_end, _changed = shift_sequence(
            seq=seq_to_check, post_seq=post_seq, pre_seq=pre_seq, start=self.start, end=self.end
        )

        if _changed:
            if using_ref:
                self.reference_allele = shifted_seq
            else:
                self.alternate_allele = shifted_seq

            self.start = new_start
            self.end = new_end

        self._original_locus = self._original_locus or original_locus
        self._shifted_locus = self._locus
        self._shifted = True

    @property
    def overlap_consequences(self) -> list[OverlapConsequence]:
        """Public accessor for _overlap_consequences.

        Returns:
            list[OverlapConsequence]: The list of overlap consequences for this variation allele.
        """
        return self._overlap_consequences

    @cached_property
    def _overlap_consequences(self) -> list[OverlapConsequence]:
        """Returns the list of overlap consequences for this variation allele.

        Returns:
            list[OverlapConsequence]: A list of overlap consequences for this variation allele.
        """
        feature_overlaps = self.get_all_featureoverlaps()
        consequences = []
        for feature_overlap in feature_overlaps:
            consequences.extend(feature_overlap.overlap_consequences)
        return consequences

    def get_all_featureoverlaps(self) -> Sequence[FeatureOverlap]:
        """Returns all feature overlaps for this variation allele.

        Returns:
            list[FeatureOverlap]: A list of all feature overlaps for this variation allele.
        """
        feature_overlaps = []
        for tva in self.transcript_variation_alleles:
            feature_overlaps.append(tva)
        # Similarly, we would extend with regulatory_variation_alleles,
        #  motif_variation_alleles, etc. when they are implemented.
        return feature_overlaps
