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
"""Unit tests for interval tree caches."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import Iterator

import pytest
from attrs import define, field

from ensembl.vep.caches.interval_tree_cache import FeatureIntervalTreeCache, IntervalTreeCache
from ensembl.vep.common import SourceType
from ensembl.vep.models import Feature
from ensembl.vep.sources.feature_source import FeatureSource


@define(kw_only=True, auto_attribs=True)
class DummyFeatureSource(FeatureSource[Feature]):
    """Minimal feature source for cache identity tests."""

    features: list[Feature] = field(factory=list)

    def get_all_features(self) -> list[Feature]:
        """Return all stored features."""

        return self.features

    def iter_all_features(self) -> Iterator[Feature]:
        """Return an iterator over stored features."""

        return iter(self.features)

    def get_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> list[Feature]:
        """Return all stored features."""

        return self.features

    def iter_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> Iterator[Feature]:
        """Return an iterator over stored features."""

        return iter(self.features)

    def close(self) -> None:
        """Clear stored features."""

        self.features = []


class TestFeatureIntervalTreeCache:
    """Tests for the shared interval tree cache."""

    def test_cache_implements_interval_tree_cache_interface(self) -> None:
        """Test the concrete cache satisfies the cache protocol."""

        assert isinstance(FeatureIntervalTreeCache(), IntervalTreeCache)

    def test_reuses_tree_for_same_identity(self) -> None:
        """Test identical requests return the same tree."""

        source = DummyFeatureSource(source_type=SourceType.FILE)
        cache = FeatureIntervalTreeCache()

        first = cache.get_interval_tree("chr1", 1, 100, source)
        second = cache.get_interval_tree("chr1", 1, 100, source)

        assert first is second

    def test_key_distinguishes_region_coordinates(self) -> None:
        """Test region name and coordinates are part of tree identity."""

        source = DummyFeatureSource(source_type=SourceType.FILE)
        cache = FeatureIntervalTreeCache()

        original = cache.get_interval_tree("chr1", 1, 100, source)
        different_region = cache.get_interval_tree("chr2", 1, 100, source)
        shifted_start = cache.get_interval_tree("chr1", 2, 100, source)
        shifted_end = cache.get_interval_tree("chr1", 1, 101, source)

        assert len({id(original), id(different_region), id(shifted_start), id(shifted_end)}) == 4

    def test_key_distinguishes_sources(self) -> None:
        """Test equal regions from different feature sources do not collide."""

        first_source = DummyFeatureSource(source_type=SourceType.FILE)
        second_source = DummyFeatureSource(source_type=SourceType.FILE)
        cache = FeatureIntervalTreeCache()

        first = cache.get_interval_tree("chr1", 1, 100, first_source)
        second = cache.get_interval_tree("chr1", 1, 100, second_source)

        assert first is not second

    def test_clear_stops_reusing_cached_tree(self) -> None:
        """Test clear removes entries even when consumers retain trees."""

        source = DummyFeatureSource(source_type=SourceType.FILE)
        cache = FeatureIntervalTreeCache()
        first = cache.get_interval_tree("chr1", 1, 100, source)

        cache.clear()
        second = cache.get_interval_tree("chr1", 1, 100, source)

        assert first is not second

    @pytest.mark.parametrize(("start", "end"), [(0, 1), (2, 1)])
    def test_rejects_invalid_coordinates(self, start: int, end: int) -> None:
        """Test invalid slice coordinates are rejected."""

        source = DummyFeatureSource(source_type=SourceType.FILE)

        with pytest.raises(ValueError, match="slice start"):
            FeatureIntervalTreeCache().get_interval_tree("chr1", start, end, source)

    def test_concurrent_requests_share_tree(self) -> None:
        """Test concurrent misses converge on one canonical interval tree."""

        source = DummyFeatureSource(source_type=SourceType.FILE)
        cache = FeatureIntervalTreeCache()

        with ThreadPoolExecutor(max_workers=4) as executor:
            trees = list(
                executor.map(lambda _index: cache.get_interval_tree("chr1", 1, 100, source), range(20))
            )

        assert all(tree is trees[0] for tree in trees)
