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
"""Colocated-variant overlap finding."""

from __future__ import annotations

from attrs import define

from ensembl.vep.models.overlap import FeatureOverlap
from ensembl.vep.models.variation import Variation, VariationAllele
from ensembl.vep.overlaps.base import CandidateLookup, FeatureMatcher
from ensembl.vep.overlaps.feature import FeatureOverlapFinder

__all__ = ["ColocatedVariantOverlapFinder", "colocated_variant_overlap"]


def colocated_variant_overlap(variation: Variation, variation_allele: VariationAllele) -> FeatureOverlap:
    """Construct an overlap between a known variation and an input allele."""
    return FeatureOverlap(feature=variation, overlapping_feature=variation_allele)


@define(kw_only=True, auto_attribs=True)
class ColocatedVariantOverlapFinder(FeatureOverlapFinder[Variation, FeatureOverlap]):
    """Find known variants colocated with input variation alleles."""

    candidates: CandidateLookup[Variation]
    matcher: FeatureMatcher

    def __attrs_post_init__(self) -> None:
        """Use the colocated-variant overlap factory."""
        self.overlap_factory = colocated_variant_overlap
