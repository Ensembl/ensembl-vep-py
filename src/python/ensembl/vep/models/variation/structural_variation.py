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
from typing import Sequence
from attrs import define, field

from ensembl.vep.common import ClassSOTerm, infer_symbolic_variant_so_term, invalidates
from ensembl.vep.models.variation.variation import Variation
from ensembl.vep.models.variation.structural_variation_allele import StructuralVariationAllele


@define(kw_only=True, auto_attribs=True)
class StructuralVariation(Variation):
    """Represents a structural variation to genomic sequences in a specific locus.

    Attributes:
        sv_type (str | None):
            Type of the structural variant.
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
        variation_alleles (Sequence[StructuralVariationAllele]):
            List of VariationAllele objects describing the alleles of the variation.
    """

    sv_type: str | None = field(default=None, on_setattr=invalidates("_class_so_term"))
    copy_number: int | None = field(default=None)
    inner_start: int | None = field(default=None)
    outer_start: int | None = field(default=None)
    inner_end: int | None = field(default=None)
    outer_end: int | None = field(default=None)
    variation_alleles: Sequence[StructuralVariationAllele] = field(
        factory=list,
        init=False,
        on_setattr=invalidates("_most_severe_consequence"),
    )

    @cached_property
    def _class_so_term(self) -> ClassSOTerm | None:
        """Return the inferred SO term for this structural variation."""
        allele_string = self.allele_string or ""
        _, separator, alternate_alleles = allele_string.partition("/")
        raw_type = alternate_alleles if separator else allele_string
        return infer_symbolic_variant_so_term(raw_type, svtype=self.sv_type)
