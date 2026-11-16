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
"""Custom-feature overlap finding."""

from __future__ import annotations

from attrs import define

from ensembl.vep.models.custom_feature import CustomFeature
from ensembl.vep.models.overlap import CustomVariationAllele
from ensembl.vep.models.variation import VariationAllele
from ensembl.vep.overlaps.base import CandidateLookup, FeatureMatcher
from ensembl.vep.overlaps.feature import FeatureOverlapFinder

__all__ = ["CustomOverlapFinder", "custom_overlap"]


def custom_overlap(custom_feature: CustomFeature, variation_allele: VariationAllele) -> CustomVariationAllele:
    """Construct a typed custom-feature overlap."""
    return CustomVariationAllele(feature=custom_feature, overlapping_feature=variation_allele)


@define(kw_only=True, auto_attribs=True)
class CustomOverlapFinder(FeatureOverlapFinder[CustomFeature, CustomVariationAllele]):
    """Find custom feature overlaps for a variation-allele buffer."""

    candidates: CandidateLookup[CustomFeature]
    matcher: FeatureMatcher

    def __attrs_post_init__(self) -> None:
        """Use the typed custom-overlap factory."""
        self.overlap_factory = custom_overlap
