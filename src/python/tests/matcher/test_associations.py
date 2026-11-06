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
"""Tests for exact and reciprocal association matchers."""

from __future__ import annotations

from attrs import define
import pytest

from ensembl.vep.matcher import ExactCoordinateMatcher, ReciprocalOverlapMatcher, get_matched_variant_alleles


@define(kw_only=True, auto_attribs=True)
class _Interval:
    """Minimal interval used to test source-independent matchers."""

    start: int
    end: int


@define(kw_only=True, auto_attribs=True)
class _Allele(_Interval):
    """Minimal allele used to test plugin-facing lookup."""

    reference_allele: str | None = None
    alternate_allele: str | None = None


class TestExactCoordinateMatcher:
    """Exact coordinate association."""

    def test_matches_identical_normal_intervals(self) -> None:
        """Match equal start and end coordinates."""
        matcher = ExactCoordinateMatcher()

        assert matcher.matches(_Interval(start=10, end=20), _Interval(start=10, end=20))

    def test_rejects_partially_overlapping_intervals(self) -> None:
        """Require identical coordinates rather than ordinary overlap."""
        matcher = ExactCoordinateMatcher()

        assert not matcher.matches(_Interval(start=10, end=20), _Interval(start=20, end=30))

    def test_matches_identical_insertion_positions(self) -> None:
        """Treat an interbase insertion position as an exact coordinate pair."""
        matcher = ExactCoordinateMatcher()

        assert matcher.matches(_Interval(start=10, end=9), _Interval(start=10, end=9))


class TestReciprocalOverlapMatcher:
    """Reference-span fraction association."""

    def test_applies_independent_fraction_thresholds(self) -> None:
        """Require each side to meet its configured coverage fraction."""
        matcher = ReciprocalOverlapMatcher(left_minimum=0.5, right_minimum=0.75)

        assert matcher.matches(_Interval(start=10, end=19), _Interval(start=15, end=19))
        assert not matcher.matches(_Interval(start=10, end=19), _Interval(start=15, end=24))

    def test_rejects_zero_width_insertions(self) -> None:
        """Keep reciprocal fractions undefined for insertion coordinates."""
        matcher = ReciprocalOverlapMatcher(left_minimum=0.0, right_minimum=0.0)

        assert not matcher.matches(_Interval(start=10, end=9), _Interval(start=10, end=9))

    @pytest.mark.parametrize("fraction", [-0.1, 1.1])
    def test_rejects_invalid_fraction_thresholds(self, fraction: float) -> None:
        """Require configured fractions to be probabilities."""
        with pytest.raises(ValueError, match="between zero and one"):
            ReciprocalOverlapMatcher(left_minimum=fraction)


class TestMatchedVariantAlleles:
    """Plugin-facing exact allele lookup."""

    def test_retains_candidates_matching_coordinates_and_alleles(self) -> None:
        """Return matching alleles in their input order."""
        query = _Allele(start=10, end=10, reference_allele="A", alternate_allele="T")
        candidates = [
            _Allele(start=10, end=10, reference_allele="A", alternate_allele="T"),
            _Allele(start=10, end=10, reference_allele="A", alternate_allele="G"),
            _Allele(start=11, end=11, reference_allele="A", alternate_allele="T"),
        ]

        matched = get_matched_variant_alleles(candidates, query)

        assert matched == [candidates[0]]

    def test_can_ignore_alternate_allele_identity(self) -> None:
        """Allow coordinate and reference-only matching when requested."""
        query = _Allele(start=10, end=10, reference_allele="A", alternate_allele="T")
        candidate = _Allele(start=10, end=10, reference_allele="A", alternate_allele="G")

        matched = get_matched_variant_alleles([candidate], query, match_alternate=False)

        assert matched == [candidate]
