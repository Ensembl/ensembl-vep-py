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

import re
from functools import cached_property
from threading import RLock
from attrs import define, field

from ensembl.vep.common import ClassSOTerm, invalidate, invalidates
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
        on_setattr=invalidates(
            "__variation_alleles", "__class_so_term_ref_correct", "__most_severe_consequence"
        ),
    )
    is_somatic: bool = field(default=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    @property
    def variation_alleles(self) -> list[VariationAllele]:
        """Public accessor for the __variation_alleles.

        Returns:
            list[VariationAllele]: The list of VariationAllele objects for the variation.
        """
        return self.__variation_alleles

    @cached_property
    def __variation_alleles(self) -> list[VariationAllele]:
        """Returns the list of VariationAllele objects for the variation.

        Returns:
            list[VariationAllele]: The list of VariationAllele objects for the variation.
        """
        invalidate(self, "__most_severe_consequence")

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
        """Public accessor for the __class_so_term.

        Returns:
            ClassSOTerm | None: The SO term for the class of the variation.
        """
        return self.__class_so_term_ref_correct

    @cached_property
    def __class_so_term_ref_correct(self) -> ClassSOTerm | None:
        """
        Return the SO term for the class of the variation.

        Args:
            ref_correct: Whether the first allele is known to be the reference.

        Returns:
            ClassSOTerm | None: The SO term string for the variation.
        """
        allele_string = self.allele_string or ""

        # default to sequence_alteration
        so_term = ClassSOTerm.SEQUENCE_ALTERATION

        # A/T, A/T/G
        if re.match(r"^[A-Z](/[A-Z])+$", allele_string):
            so_term = ClassSOTerm.SNV

        # (CAG)8/(CAG)9
        elif re.search(r"\)\d+", allele_string):
            so_term = ClassSOTerm.TANDEM_REPEAT

        else:
            alleles = allele_string.split("/")

            if len(alleles) > 1:
                ref = alleles[0]
                alts = alleles[1:]

                if ref == "-":
                    if len(alts) == 1 and re.search(r"DEL", alts[0]):
                        # -/(LARGEDELETION)
                        so_term = ClassSOTerm.DELETION

                    if not any(not re.match(r"^[A-Z]+$|INS", alt) for alt in alts):
                        # -/ATT, -/(LARGEINSERTION)
                        so_term = ClassSOTerm.INSERTION

                elif re.match(r"^[A-Z]+$", ref):
                    if not any(not re.match(r"^-$|DEL", alt) for alt in alts):
                        # A/-, A/(LARGEDELETION)
                        so_term = ClassSOTerm.DELETION

                    elif re.match(r"^[A-Z]+(/[A-Z]+)+$", allele_string):
                        # AA/TT => substitution, AA/TTTT => indel
                        same_size = all(len(alt) == len(ref) for alt in alts)
                        so_term = ClassSOTerm.SUBSTITUTION if same_size else ClassSOTerm.INDEL

                elif re.search(r"DEL", ref):
                    if not any(alt != "-" for alt in alts):
                        # (LARGEDELETION)/-, (2345 BP DELETION)/-
                        so_term = ClassSOTerm.DELETION

            elif len(alleles) == 1:
                if re.search(r"INS", alleles[0]):
                    # (LARGEINSERTION)
                    so_term = ClassSOTerm.INSERTION
                elif re.search(r"DEL", alleles[0]):
                    # (308 BP DELETION)
                    so_term = ClassSOTerm.DELETION

        return so_term

    # TBD when and how to call __class_so_term_ref_incorrect
    # instead of __class_so_term_ref_correct
    @cached_property
    def __class_so_term_ref_incorrect(self) -> ClassSOTerm | None:  # pylint: disable=unused-private-member
        """
        Return the SO term for the class of the variation
        with the first allele kown not to be the reference.

        Returns:
            ClassSOTerm | None: The SO term string for the variation.
        """
        allele_string = self.allele_string or ""

        # default to sequence_alteration
        so_term = ClassSOTerm.SEQUENCE_ALTERATION

        # A/T, A/T/G
        if re.match(r"^[A-Z](/[A-Z])+$", allele_string):
            so_term = ClassSOTerm.SNV

        # (CAG)8/(CAG)9
        elif re.search(r"\)\d+", allele_string):
            so_term = ClassSOTerm.TANDEM_REPEAT

        else:
            alleles = allele_string.split("/")

            if len(alleles) > 1:
                ref = alleles[0]
                alts = alleles[1:]

                if ref == "-":
                    if len(alts) == 1 and re.search(r"DEL", alts[0]):
                        # -/(LARGEDELETION)
                        so_term = ClassSOTerm.INDEL

                    if not any(not re.match(r"^[A-Z]+$|INS", alt) for alt in alts):
                        # -/ATT, -/(LARGEINSERTION)
                        so_term = ClassSOTerm.INDEL

                elif re.match(r"^[A-Z]+$", ref):
                    if not any(not re.match(r"^-$|DEL", alt) for alt in alts):
                        # A/-, A/(LARGEDELETION)
                        so_term = ClassSOTerm.INDEL

                    elif re.match(r"^[A-Z]+(/[A-Z]+)+$", allele_string):
                        # AA/TT => substitution, AA/TTTT => indel
                        same_size = all(len(alt) == len(ref) for alt in alts)
                        so_term = ClassSOTerm.SUBSTITUTION if same_size else ClassSOTerm.INDEL

                elif re.search(r"DEL", ref):
                    if not any(alt != "-" for alt in alts):
                        # (LARGEDELETION)/-, (2345 BP DELETION)/-
                        so_term = ClassSOTerm.INDEL

            elif len(alleles) == 1:
                if re.search(r"INS", alleles[0]):
                    # (LARGEINSERTION)
                    so_term = ClassSOTerm.INDEL
                elif re.search(r"DEL", alleles[0]):
                    # (308 BP DELETION)
                    so_term = ClassSOTerm.INDEL

        return so_term

    def minimise_alleles(self) -> None:
        """Minimises all alleles."""
        for variation_allele in self.variation_alleles or []:
            variation_allele.minimise()

    def shift_alleles(self) -> None:
        """Shifts all alleles."""
        for variation_allele in self.variation_alleles or []:
            variation_allele.shift()
