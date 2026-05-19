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
"""StructuralVariation class."""

from functools import cached_property
from threading import RLock
from attrs import define, field

from ensembl.vep.common import ClassSOTerm, invalidate, invalidates, infer_variant_so_term
from ensembl.vep.models.overlap import OverlapConsequence
from ensembl.vep.models.variation.variation import Variation
from ensembl.vep.models.variation.structural_variation_allele import StructuralVariationAllele


@define(kw_only=True, auto_attribs=True)
class StructuralVariation(Variation):
    """Represents a structural variation to genomic sequences in a specific locus.

    Attributes:
            Stable identifier for the variation.
        allele_string (str | None):
            String representation of the alleles.
        structural_variation_alleles (list[StructuralVariationAllele] | None):
            List of StructuralVariationAllele objects describing the alleles of the variation.
        is_somatic (bool | None):
            Indicates whether the variation is somatic.
        copy_number (float | None):
            Number of sequence copies of a structural variant classified as "copy number variant" (CNV) when
            available.
        inner_start (int | None):
            The inner bound on the start position when the variant is imprecise
        outer_start (int | None):
            The outer bound on the start position when the variant is imprecise
        inner_end (int | None):
            The inner bound on the end position when the variant is imprecise
        outer_end (int | None):
            The outer bound on the end position when the variant is imprecise
            Sequence Ontology term describing the class of the variation.
        most_severe_consequence (OverlapConsequence | None):
            Most severe consequence type for the variation.
    """

    allele_string: str | None = field(
        default=None,
        on_setattr=invalidates(
            "__structural_variation_alleles", "__class_so_term", "__most_severe_consequence"
        ),
    )
    is_somatic: bool = field(default=False)
    copy_number: int | None = field(default=None)
    inner_start: int | None = field(default=None)
    outer_start: int | None = field(default=None)
    inner_end: int | None = field(default=None)
    outer_end: int | None = field(default=None)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    @property
    def structural_variation_alleles(self) -> list[StructuralVariationAllele]:
        """Public accessor for the __structural_variation_alleles.

        Returns:
            list[StructuralVariationAllele]: The list of StructuralVariationAllele objects for the variation.
        """
        return self.__structural_variation_alleles

    @cached_property
    def __structural_variation_alleles(self) -> list[StructuralVariationAllele]:
        """Returns the list of StructuralVariationAllele objects for the variation.

        Returns:
            list[StructuralVariationAllele]: The list of StructuralVariationAllele objects for the variation.
        """
        invalidate(self, "__most_severe_consequence")

        with self._lock:
            structural_variation_alleles = []
            alleles = (self.allele_string or "").split("/")

            ref = alleles[0]
            alts = alleles[1:]
            for alt in alts:
                # ref is still same for all alleles, so we can use the start and end for all alleles
                structural_variation_allele = StructuralVariationAllele(
                    start=self.start,
                    end=self.end,
                    strand=self.strand,
                    slice=self.slice,
                    reference_allele=ref,
                    alternate_allele=alt,
                )

                structural_variation_alleles.append(structural_variation_allele)

            return structural_variation_alleles

    @property
    def most_severe_consequence(self) -> OverlapConsequence | None:
        """Public accessor for the __most_severe_consequence.

        Returns:
            OverlapConsequence | None: The most severe consequence type for the variation
                or None if not deducible.
        """
        return self.__most_severe_consequence

    @cached_property
    def __most_severe_consequence(self) -> OverlapConsequence | None:
        """Returns the most severe consequence type for the variation.

        Returns:
            OverlapConsequence | None: The most severe consequence type or None if not deducible.
        """
        highest = None
        for structural_variation_allele in self.structural_variation_alleles or []:
            for overlap_consequence in structural_variation_allele.overlap_consequences or []:
                highest = highest or overlap_consequence
                if overlap_consequence.rank < highest.rank:
                    highest = overlap_consequence
                elif overlap_consequence.rank == highest.rank and overlap_consequence.term < highest.term:
                    highest = overlap_consequence

        return highest

