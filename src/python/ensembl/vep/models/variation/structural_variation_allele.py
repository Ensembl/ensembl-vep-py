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
"""StructuralVariationAllele class."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING
from attrs import define, field

from ensembl.vep.models import Feature

if TYPE_CHECKING:
    from ensembl.vep.models.overlap import TranscriptStructuralVariationAllele

logger = logging.getLogger(__name__)


@define(kw_only=True, auto_attribs=True)
class StructuralVariationAllele(VariationAllele):
    """Represents a specific allele of a structural variation."""

    start: int
    end: int
    reference_allele: str | None
    alternate_allele: str | None
    transcript_structural_variation_alleles: list[TranscriptStructuralVariationAllele] = field(
        factory=list, init=False
    )
    # regulatory_variation_alleles: list[RegulatoryVariationAllele] = field(factory=list, init=False)
    # motif_variation_alleles: list[MotifVariationAllele] = field(factory=list, init=False)
    # intergenic_variation_alleles: list[IntergenicVariationAllele] = field(factory=list, init=False)
