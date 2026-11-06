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
"""Tests for general overlap matching."""

from __future__ import annotations

from attrs import define
import pytest

from ensembl.vep.common.enums import Strand
from ensembl.vep.matcher import InsertionMatchPolicy, OverlapMatcher


@define(kw_only=True, auto_attribs=True)
class _Feature:
    """Minimal interval-bearing feature."""

    start: int
    end: int
    strand: Strand = Strand.FORWARD


class TestOverlapMatcher:
    """General overlap relationships."""

    @pytest.mark.parametrize(
        ("left", "right"),
        [
            pytest.param((10, 30), (15, 20), id="contain"),
            pytest.param((15, 20), (10, 30), id="within"),
            pytest.param((10, 20), (15, 25), id="intersect"),
        ],
    )
    def test_matches_every_overlap_relationship(
        self,
        left: tuple[int, int],
        right: tuple[int, int],
    ) -> None:
        """Match contain, within, and partial intersection relationships."""

        matcher = OverlapMatcher()

        assert matcher.matches(
            _Feature(start=left[0], end=left[1]),
            _Feature(start=right[0], end=right[1]),
        )

    def test_rejects_disjoint_normal_intervals(self) -> None:
        """Reject normal intervals without a shared base."""

        matcher = OverlapMatcher()

        assert not matcher.matches(_Feature(start=10, end=20), _Feature(start=21, end=30))

    @pytest.mark.parametrize(
        ("policy", "insertion_start", "expected"),
        [
            pytest.param(InsertionMatchPolicy.CONTAINED, 10, False, id="contained_start_boundary"),
            pytest.param(InsertionMatchPolicy.CONTAINED, 11, True, id="contained_interior"),
            pytest.param(InsertionMatchPolicy.CONTAINED, 21, False, id="contained_end_boundary"),
            pytest.param(InsertionMatchPolicy.TOUCHING, 10, True, id="touching_start_boundary"),
            pytest.param(InsertionMatchPolicy.TOUCHING, 21, True, id="touching_end_boundary"),
        ],
    )
    def test_applies_insertion_boundary_policy(
        self,
        policy: InsertionMatchPolicy,
        insertion_start: int,
        expected: bool,
    ) -> None:
        """Apply the selected relationship between an insertion and a feature."""

        matcher = OverlapMatcher(insertion_policy=policy)
        insertion = _Feature(start=insertion_start, end=insertion_start - 1)
        feature = _Feature(start=10, end=20)

        assert matcher.matches(insertion, feature) is expected

    @pytest.mark.parametrize(
        ("strand", "expected"),
        [
            pytest.param(Strand.FORWARD, True, id="forward"),
            pytest.param(Strand.REVERSE, False, id="reverse"),
            pytest.param(Strand.UNKNOWN, False, id="unknown"),
        ],
    )
    def test_uses_feature_strand_for_three_prime_policy(self, strand: Strand, expected: bool) -> None:
        """Project an insertion onto the feature-relative 3-prime flank."""

        matcher = OverlapMatcher(insertion_policy=InsertionMatchPolicy.THREE_PRIME)
        insertion = _Feature(start=10, end=9)
        feature = _Feature(start=10, end=20, strand=strand)

        assert matcher.matches(insertion, feature) is expected

    def test_matches_insertions_at_the_same_interbase_position(self) -> None:
        """Match two insertions sharing an interbase coordinate."""

        matcher = OverlapMatcher()

        assert matcher.matches(_Feature(start=10, end=9), _Feature(start=10, end=9))
