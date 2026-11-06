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
"""Tests for the binned feature cache."""

from __future__ import annotations

from ensembl.vep.caches import BinnedCache, BinnedFeatureCache, Interval


class TestBinnedFeatureCache:
    """Sequence-region candidate retrieval."""

    def test_satisfies_the_binned_cache_protocol(self) -> None:
        """Serve the interface its consumers depend on."""

        assert isinstance(BinnedFeatureCache[str](), BinnedCache)

    def test_returns_features_from_matching_sequence_region(self) -> None:
        """Restrict cached candidates to their sequence region."""

        cache = BinnedFeatureCache[str]()
        cache.add("1", "first", Interval(start=100, end=200))
        cache.add("2", "second", Interval(start=100, end=200))

        objects = cache.get_objects("1", [Interval(start=150, end=150)])

        assert objects == ["first"]

    def test_finds_large_feature_across_bin_boundaries(self) -> None:
        """Find a large feature from either side of a bin boundary."""

        cache = BinnedFeatureCache[str]()
        cache.add("1", "feature", Interval(start=100_000, end=300_000))

        objects = cache.get_objects("1", [Interval(start=120_000, end=280_000)])

        assert objects == ["feature"]
        assert len(cache) == 1

    def test_finds_a_feature_wider_than_every_bin_level(self) -> None:
        """Retain a feature too wide for any bin level in the top-level bin."""

        cache = BinnedFeatureCache[str]()
        cache.add("1", "chromosome_wide", Interval(start=1, end=200_000_000))

        objects = cache.get_objects("1", [Interval(start=150_000_000, end=150_000_000)])

        assert objects == ["chromosome_wide"]

    def test_ignores_repeated_logical_features(self) -> None:
        """Avoid duplicate annotations after overlapping source fetches."""

        cache = BinnedFeatureCache[str]()
        cache.add("1", "first", Interval(start=100, end=200), key="feature-1")
        cache.add("1", "reloaded", Interval(start=100, end=200), key="feature-1")

        objects = cache.get_objects("1", [Interval(start=100, end=200), Interval(start=150, end=150)])

        assert objects == ["first"]

    def test_expands_insertions_for_candidate_retrieval(self) -> None:
        """Return boundary features for an insertion point."""

        cache = BinnedFeatureCache[str]()
        cache.add("1", "left", Interval(start=1, end=9))
        cache.add("1", "right", Interval(start=10, end=20))

        objects = cache.get_objects("1", [Interval(start=10, end=9)])

        assert set(objects) == {"left", "right"}

    def test_clamps_insertion_before_first_base_for_candidate_retrieval(self) -> None:
        """Find a feature starting at the first base."""

        cache = BinnedFeatureCache[str]()
        cache.add("1", "first", Interval(start=1, end=1))

        objects = cache.get_objects("1", [Interval(start=1, end=0)])

        assert objects == ["first"]

    def test_clear_removes_cached_features(self) -> None:
        """Remove all cached feature regions."""

        cache = BinnedFeatureCache[str]()
        cache.add("1", "feature", Interval(start=1, end=10))
        cache.clear()

        assert not cache.get_objects("1", [Interval(start=1, end=10)])
