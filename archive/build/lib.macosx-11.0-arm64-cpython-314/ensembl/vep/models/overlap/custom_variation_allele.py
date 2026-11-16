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
"""Custom-feature overlap model."""

from __future__ import annotations

from typing import TYPE_CHECKING

from attrs import define, field

from ensembl.vep.models.custom_feature import CustomFeature
from ensembl.vep.models.overlap.feature_overlap import FeatureOverlap

if TYPE_CHECKING:
    from ensembl.vep.models.variation.variation_allele import VariationAllele

__all__ = ["CustomVariationAllele"]


@define(kw_only=True, auto_attribs=True)
class CustomVariationAllele(FeatureOverlap):
    """Represent an association between a custom feature and variation allele.

    Attributes:
        feature (CustomFeature): The custom record represented as a feature.
        overlapping_feature (VariationAllele): The associated input allele.
        summary_statistics (dict[str, object]): Track- or annotator-derived summary values.
    """

    feature: CustomFeature
    overlapping_feature: VariationAllele
    summary_statistics: dict[str, object] = field(factory=dict)

    @property
    def overlap_length(self) -> int:
        """Return the number of shared reference bases."""
        start = max(self.feature.start, self.overlapping_feature.start)
        end = min(self.feature.end, self.overlapping_feature.end)
        return max(0, end - start + 1)

    @property
    def feature_overlap_fraction(self) -> float:
        """Return the fraction of the custom feature sharing reference bases."""
        feature_length = max(0, self.feature.end - self.feature.start + 1)
        if feature_length == 0:
            return 0.0
        return self.overlap_length / feature_length

    @property
    def variation_overlap_fraction(self) -> float | None:
        """Return the fraction of the allele's reference span shared with the feature.

        Insertions have no reference span, so their fraction is undefined even when the selected
        matcher associates the insertion point with this record.
        """
        variation_length = max(0, self.overlapping_feature.end - self.overlapping_feature.start + 1)
        if variation_length == 0:
            return None
        return self.overlap_length / variation_length
