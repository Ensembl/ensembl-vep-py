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
"""Tests for custom variation overlaps."""

from __future__ import annotations

from ensembl.vep.common import Strand
from ensembl.vep.models.custom_feature import CustomFeature
from ensembl.vep.models.overlap import CustomVariationAllele
from ensembl.vep.models.variation import VariationAllele


class TestCustomFeature:
    """Custom source-record payload ownership."""

    def test_uses_an_independent_record_payload_per_feature(self) -> None:
        """Avoid sharing mutable source records between custom features."""
        first = CustomFeature(start=10, end=20, strand=Strand.FORWARD)
        second = CustomFeature(start=10, end=20, strand=Strand.FORWARD)

        first.record["score"] = 5

        assert second.record == {}


class TestCustomVariationAllele:
    """Reference-span geometry and summary payloads."""

    def test_reports_overlap_fractions_for_normal_intervals(self) -> None:
        """Calculate each side's fraction from shared reference bases."""
        overlap = CustomVariationAllele(
            feature=CustomFeature(start=10, end=19, strand=Strand.FORWARD, record={"score": 4}),
            overlapping_feature=VariationAllele(
                start=15,
                end=24,
                strand=Strand.FORWARD,
                reference_allele="AAAAAAAAAA",
                alternate_allele="T",
            ),
            summary_statistics={"mean": 3.5},
        )

        assert overlap.overlap_length == 5
        assert overlap.feature_overlap_fraction == 0.5
        assert overlap.variation_overlap_fraction == 0.5
        assert overlap.feature.record == {"score": 4}
        assert overlap.summary_statistics == {"mean": 3.5}

    def test_marks_insertion_fraction_as_undefined(self) -> None:
        """Avoid inventing a reference-span percentage for an insertion."""
        overlap = CustomVariationAllele(
            feature=CustomFeature(start=10, end=20, strand=Strand.FORWARD),
            overlapping_feature=VariationAllele(
                start=15,
                end=14,
                strand=Strand.FORWARD,
                reference_allele="",
                alternate_allele="T",
            ),
        )

        assert overlap.overlap_length == 0
        assert overlap.feature_overlap_fraction == 0.0
        assert overlap.variation_overlap_fraction is None
