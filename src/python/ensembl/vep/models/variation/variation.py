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
"""Variation class."""

from __future__ import annotations

from functools import cached_property
from threading import RLock
from typing import Sequence
from attrs import define, field

from ensembl.vep.common import ClassSOTerm, Locus, invalidate, invalidates, infer_variant_so_term
from ensembl.vep.models import Feature
from ensembl.vep.models.overlap import OverlapConsequence
from ensembl.vep.models.variation.variation_allele import VariationAllele


@define(kw_only=True, auto_attribs=True)
class Variation(Feature):
    """Represents a variation to genomic sequences in a specific locus.

    Attributes:
        stable_id (str | None):
            Stable identifier for the variation.
        allele_string (str | None):
            String representation of the alleles.
        uploaded_locus (Locus | None):
            Original position and allele string supplied by the variation source.
        variation_alleles (list[VariationAllele] | None):
            List of VariationAllele objects describing the alleles of the variation.
        is_somatic (bool | None):
            Indicates whether the variation is somatic.
        class_so_term (str | None):
            Sequence Ontology term describing the class of the variation.
        most_severe_consequence (OverlapConsequence | None):
            Most severe consequence type for the variation.
    """

    stable_id: str | None = field(default=None)
    allele_string: str | None = field(
        default=None,
        on_setattr=invalidates("_variation_alleles", "_class_so_term", "_most_severe_consequence"),
    )
    uploaded_locus: Locus | None = field(default=None)
    is_somatic: bool = field(default=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    @property
    def variation_alleles(self) -> Sequence[VariationAllele]:
        """Public accessor for the _variation_alleles.

        Returns:
            Sequence[VariationAllele]: The list of VariationAllele objects for the variation.
        """
        return self._variation_alleles

    @cached_property
    def _variation_alleles(self) -> Sequence[VariationAllele]:
        """Returns the list of VariationAllele objects for the variation.

        Returns:
            Sequence[VariationAllele]: The list of VariationAllele objects for the variation.
        """
        invalidate(self, "_most_severe_consequence")

        with self._lock:
            variation_alleles = []
            alleles = (self.allele_string or "").split("/")

            ref = alleles[0]
            alts = alleles[1:]
            for alt in alts:
                # ref is still same for all alleles, so we can use the start and end for all alleles
                variation_allele = VariationAllele(
                    start=self.start,
                    end=self.end,
                    strand=self.strand,
                    slice=self.slice,
                    reference_allele=ref,
                    alternate_allele=alt,
                )

                variation_alleles.append(variation_allele)

            return variation_alleles

    @property
    def most_severe_consequence(self) -> OverlapConsequence | None:
        """Public accessor for the _most_severe_consequence.

        Returns:
            OverlapConsequence | None: The most severe consequence type for the variation
                or None if not deducible.
        """
        return self._most_severe_consequence

    @cached_property
    def _most_severe_consequence(self) -> OverlapConsequence | None:
        """Returns the most severe consequence type for the variation.

        Returns:
            OverlapConsequence | None: The most severe consequence type or None if not deducible.
        """
        highest = None
        for variation_allele in self.variation_alleles or []:
            for overlap_consequence in variation_allele.overlap_consequences or []:
                highest = highest or overlap_consequence
                if overlap_consequence.rank < highest.rank:
                    highest = overlap_consequence
                elif overlap_consequence.rank == highest.rank and overlap_consequence.term < highest.term:
                    highest = overlap_consequence

        return highest

    @property
    def class_so_term(self) -> ClassSOTerm | None:
        """Public accessor for the _class_so_term.

        Returns:
            ClassSOTerm | None: The SO term for the class of the variation.
        """
        return self._class_so_term

    @cached_property
    def _class_so_term(self) -> ClassSOTerm | None:
        """
        Return the SO term for the class of the variation.

        Returns:
            ClassSOTerm | None: The SO term string for the variation.
        """
        return infer_variant_so_term(self.allele_string or "", ref_correct=True)

    def minimise_alleles(self) -> None:
        """Minimises all alleles."""
        for variation_allele in self.variation_alleles or []:
            variation_allele.minimise()

    def shift_alleles(self) -> None:
        """Shifts all alleles."""
        for variation_allele in self.variation_alleles or []:
            variation_allele.shift()
