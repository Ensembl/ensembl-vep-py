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
"""Regulatory-overlap finding."""

from __future__ import annotations

from attrs import define

from ensembl.vep.models.overlap import FeatureOverlap
from ensembl.vep.models.regulatory import Regulatory
from ensembl.vep.models.variation import VariationAllele
from ensembl.vep.overlaps.base import CandidateLookup, FeatureMatcher
from ensembl.vep.overlaps.feature import FeatureOverlapFinder

__all__ = ["RegulatoryOverlapFinder", "regulatory_overlap"]


def regulatory_overlap(regulatory: Regulatory, variation_allele: VariationAllele) -> FeatureOverlap:
    """Construct a generic regulatory overlap."""
    return FeatureOverlap(feature=regulatory, overlapping_feature=variation_allele)


@define(kw_only=True, auto_attribs=True)
class RegulatoryOverlapFinder(FeatureOverlapFinder[Regulatory, FeatureOverlap]):
    """Find regulatory overlaps for a variation-allele buffer."""

    candidates: CandidateLookup[Regulatory]
    matcher: FeatureMatcher

    def __attrs_post_init__(self) -> None:
        """Use the regulatory-overlap factory."""
        self.overlap_factory = regulatory_overlap
