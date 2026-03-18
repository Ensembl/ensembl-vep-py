from __future__ import annotations

import re
from attrs import define, field
from typing import TYPE_CHECKING
from functools import cached_property
from threading import Lock

from ensembl.vep.common import ClassSOTerm
from ensembl.vep.models import Feature
from ensembl.vep.models.overlap import OverlapConsequence
if TYPE_CHECKING:
    from ensembl.vep.models.variation import VariationAllele


@define(kw_only=True, auto_attribs=True)
class Variation(Feature):
    """Represents an Ensembl-specific sequence variation describing a location on a sequence in an arbitrary coordinate system.
    
    Args:
        stable_id (str | None): Stable identifier for the variation.
        allele_string (str | None): String representation of the alleles.
        variation_alleles (VariationAllele | None): VariationAllele object describing the alleles of the variation.
        is_somatic (bool | None): Indicates whether the variation is somatic.
        class_so_term (str | None): Sequence Ontology term describing the class of the variation.
        
    Attributes:
        stable_id (str | None): Stable identifier for the variation.
        allele_string (str | None): String representation of the alleles.
        variation_alleles (VariationAllele | None): VariationAllele object describing the alleles of the variation.
        is_somatic (bool | None): Indicates whether the variation is somatic.
        class_so_term (str | None): Sequence Ontology term describing the class of the variation.
        most_severe_consequence (OverlapConsequence | None): Most severe consequence type for the variation.
    """
    stable_id: str | None = None
    allele_string: str | None = None
    is_somatic: bool | None = False
    minimise_alleles: bool = False
    shift_alleles: bool = False
    _lock: Lock = field(factory=Lock, init=False, repr=False)
    _DERIVED = frozenset({
        "variation_alleles",
        "most_severe_consequence",
        "class_so_term"
    })
    _INVALIDATES = {
        "variation_alleles":["most_severe_consequence"],
        "allele_string":["variation_alleles", "class_so_term"]
    }

    def __setattr__(self, name, value):
        if name == "minimise_alleles":
            if getattr(self, "minimise_alleles", False) and not value:
                raise AttributeError(f"Cannot reset '{name}' after it has been set.")
            elif value:
                for variation_allele in self.variation_alleles or []:
                    variation_allele.minimise()
        if name == "shift_alleles":
            if getattr(self, "shift_alleles", False) and not value:
                raise AttributeError(f"Cannot reset '{name}' after it has been set.")
            elif value:
                for variation_allele in self.variation_alleles or []:
                    variation_allele.shift()
        return super().__setattr__(name, value)

    @cached_property
    def variation_alleles(self) -> list[VariationAllele]:
        """Returns the list of VariationAllele objects for the variation.
        
        Returns:
            list[VariationAllele]: List of VariationAllele objects.
        """
        from ensembl.vep.models.variation import VariationAllele # local import to avoid circular import

        with self._lock:
            variation_alleles = []
            alleles = (self.allele_string or "").split("/")

            ref = alleles[0]
            alts = alleles[1:]
            for alt in alts:
                # ref is still same for all alleles, so we can use the start and end for all alleles
                variation_allele = VariationAllele(
                    start = self.start,
                    end = self.end,
                    strand = self.strand,
                    slice = self.slice,
                    reference_allele=ref, 
                    alternate_allele=alt
                )

                if self.minimise_alleles:
                    variation_allele.minimise()

                if self.shift_alleles:
                    variation_allele.shift()
                
                variation_alleles.append(variation_allele)

            return variation_alleles

    @cached_property
    def most_severe_consequence(self) -> OverlapConsequence | None:
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
                elif overlap_consequence.rank == highest.rank:
                    if overlap_consequence.SO_term < highest.SO_term:
                        highest = overlap_consequence
        
        return highest
    
    @cached_property
    def class_so_term(self: str, ref_correct: bool = True) -> ClassSOTerm | None:
        """
        Return the SO term for the class of the variation.

        Args:
            ref_correct: Whether the first allele is known to be the reference.

        Returns:
            ClassSOTerm | None: The SO term string for the variation.
        """
        allele_string = self.allele_string or ""

        # default to sequence_alteration
        cls = ClassSOTerm.SEQUENCE_ALTERATION

        # A/T, A/T/G
        if re.match(r"^[A-Z](/[A-Z])+$", allele_string):
            cls = ClassSOTerm.SNV

        # (CAG)8/(CAG)9
        elif re.search(r"\)\d+", allele_string):
            cls = ClassSOTerm.TANDEM_REPEAT

        else:
            alleles = allele_string.split("/")

            if len(alleles) > 1:
                ref = alleles[0]
                alts = alleles[1:]

                if ref == "-":
                    if len(alts) == 1 and re.search(r"DEL", alts[0]):
                        # -/(LARGEDELETION)
                        cls = ClassSOTerm.DELETION if ref_correct else ClassSOTerm.INDEL

                    if not any(
                        not re.match(r"^[A-Z]+$|INS", alt) for alt in alts
                    ):
                        # -/ATT, -/(LARGEINSERTION)
                        cls = ClassSOTerm.INSERTION if ref_correct else ClassSOTerm.INDEL

                elif re.match(r"^[A-Z]+$", ref):
                    if not any(not re.match(r"^-$|DEL", alt) for alt in alts):
                        # A/-, A/(LARGEDELETION)
                        cls = ClassSOTerm.DELETION if ref_correct else ClassSOTerm.INDEL

                    elif re.match(
                        r"^[A-Z]+(/[A-Z]+)+$", allele_string
                    ):
                        # AA/TT => substitution, AA/TTTT => indel
                        same_size = all(len(alt) == len(ref) for alt in alts)
                        cls = ClassSOTerm.SUBSTITUTION if same_size else ClassSOTerm.INDEL

                elif re.search(r"DEL", ref):
                    if not any(alt != "-" for alt in alts):
                        # (LARGEDELETION)/-, (2345 BP DELETION)/-
                        cls = ClassSOTerm.DELETION if ref_correct else ClassSOTerm.INDEL

            elif len(alleles) == 1:
                if re.search(r"INS", alleles[0]):
                    # (LARGEINSERTION)
                    cls = ClassSOTerm.INSERTION if ref_correct else ClassSOTerm.INDEL
                elif re.search(r"DEL", alleles[0]):
                    # (308 BP DELETION)
                    cls = ClassSOTerm.DELETION if ref_correct else ClassSOTerm.INDEL

        return cls