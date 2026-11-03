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

# pylint: disable=redefined-outer-name
"""Tests for feature sources."""

from __future__ import annotations

from collections.abc import Iterator
import pytest
from attrs import define, field

from ensembl.vep.caches import SliceCache, WeakSliceCache
from ensembl.vep.common import SourceType, Strand
from ensembl.vep.sources.sequence_source import SequenceSource
from ensembl.vep.models import Feature
from ensembl.vep.sources.feature_source import FeatureSource


@define(kw_only=True, auto_attribs=True)
class ConcreteFeatureSource(FeatureSource[int, Feature]):
    """Concrete feature source holding one record per position, used to test the interface.

    Only the record half is implemented, which is the point: the feature getters come from the
    interface composing these with ``records_to_features``.
    """

    records: list[int] = field(factory=list)
    region_records: list[int] = field(factory=list)
    closed: bool = field(default=False)
    sequence_source: SequenceSource | None = field(default=None)  # ignore[assignment]

    def get_all_records(self) -> list[int]:
        """Return every record."""
        return list(self.records)

    def iter_all_records(self) -> Iterator[int]:
        """Yield every record."""
        return iter(self.records)

    def get_records_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
        **kwargs: object,
    ) -> list[int]:
        """Return the records this double was given for a region."""
        del region, seq_region_name, start, end, kwargs
        return list(self.region_records)

    def iter_records_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
        **kwargs: object,
    ) -> Iterator[int]:
        """Yield the records this double was given for a region."""
        del region, seq_region_name, start, end, kwargs
        return iter(self.region_records)

    def records_to_features(self, records: list[int]) -> list[Feature]:
        """Return one single-base feature per record."""
        return [Feature(start=value, end=value, strand=Strand.FORWARD, slice=None) for value in records]

    def close(self) -> None:
        self.closed = True


class TestFeatureSourceInitialization:
    """FeatureSource initialization."""

    def test_abstract_feature_source_cannot_be_instantiated(self) -> None:
        """FeatureSource cannot be instantiated directly."""

        with pytest.raises(TypeError, match="abstract class"):
            FeatureSource(  # type: ignore[abstract]  # pylint: disable=abstract-class-instantiated
                source_type=SourceType.FILE, sequence_source=None
            )

    def test_concrete_feature_source_sets_source_type(self) -> None:
        """Concrete feature source sets the source_type correctly."""

        feature_source = ConcreteFeatureSource(source_type=SourceType.FILE)

        assert feature_source.source_type == SourceType.FILE

    def test_concrete_feature_source_has_default_slice_cache(self) -> None:
        """Feature sources receive a slice cache by default."""

        feature_source = ConcreteFeatureSource(source_type=SourceType.FILE)

        assert isinstance(feature_source.slice_cache, SliceCache)

    def test_concrete_feature_source_accepts_shared_slice_cache(self) -> None:
        """A slice cache can be injected and shared between feature sources."""

        slice_cache = WeakSliceCache()
        first = ConcreteFeatureSource(source_type=SourceType.FILE, slice_cache=slice_cache)
        second = ConcreteFeatureSource(source_type=SourceType.FILE, slice_cache=slice_cache)

        assert first.slice_cache is slice_cache
        assert second.slice_cache is slice_cache

    def test_a_source_declares_no_record_count_by_default(self) -> None:
        """A count is optional metadata a source answers only when its format carries one."""

        feature_source = ConcreteFeatureSource(source_type=SourceType.FILE)

        assert feature_source.record_count() is None
        assert feature_source.record_count(scan=True) is None


class TestFeatureSourceContextManager:
    """FeatureSource context manager."""

    @pytest.fixture
    def source(self) -> ConcreteFeatureSource:
        """Create a minimal FeatureSource instance for testing."""

        return ConcreteFeatureSource(source_type=SourceType.FILE)

    def test_enter_returns_self(self, source: ConcreteFeatureSource) -> None:
        """__enter__ returns the source instance itself."""

        with source as s:
            assert s is source

    def test_exit_calls_close(self, source: ConcreteFeatureSource) -> None:
        """__exit__ calls close on normal exit."""

        with source:
            pass

        assert source.closed is True


class TestFeatureGettersComposeTheRecordHalf:
    """What a source gets for free once it can answer in records."""

    def test_every_record_becomes_a_feature(self) -> None:
        """``get_all_features`` is every record converted."""

        feature_source = ConcreteFeatureSource(source_type=SourceType.FILE, records=[1, 2, 3])

        assert [feature.start for feature in feature_source.get_all_features()] == [1, 2, 3]

    def test_streaming_converts_one_record_at_a_time(self) -> None:
        """``iter_all_features`` yields without reading the whole file first."""

        feature_source = ConcreteFeatureSource(source_type=SourceType.FILE, records=[4, 5])

        assert [feature.start for feature in feature_source.iter_all_features()] == [4, 5]

    def test_a_region_query_converts_only_its_records(self) -> None:
        """Both region getters answer from ``get_records_by_region`` and nothing else."""

        feature_source = ConcreteFeatureSource(
            source_type=SourceType.FILE,
            records=[1, 2, 3],
            region_records=[7],
        )

        found = feature_source.get_features_by_region("1:7-7")
        streamed = feature_source.iter_features_by_region("1:7-7")

        assert [feature.start for feature in found] == [7]
        assert [feature.start for feature in streamed] == [7]
