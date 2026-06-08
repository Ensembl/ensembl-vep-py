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

from typing import Sequence
from attrs import define, field

from ensembl.vep.models.variation.variation import Variation
from ensembl.vep.models.variation.structural_variation_allele import StructuralVariationAllele


@define(kw_only=True, auto_attribs=True)
class StructuralVariation(Variation):
    """Represents a structural variation to genomic sequences in a specific locus.

    Attributes:
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
        variation_alleles (Sequence[StructuralVariationAllele]):
            List of VariationAllele objects describing the alleles of the variation.
    """

    copy_number: int | None = field(default=None)
    inner_start: int | None = field(default=None)
    outer_start: int | None = field(default=None)
    inner_end: int | None = field(default=None)
    outer_end: int | None = field(default=None)
    variation_alleles: Sequence[StructuralVariationAllele] = field(factory=list, init=False)
