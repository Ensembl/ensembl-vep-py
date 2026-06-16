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
"""Unit tests for FeatureSource."""

from typing import Iterator
<<<<<<< HEAD
import pytest
from attrs import define, field

from ensembl.vep.caches import SliceCache, WeakSliceCache
from ensembl.vep.common import SourceType
from ensembl.vep.sources.sequence_source import SequenceSource
=======

import pytest
from attrs import define, field

from ensembl.vep.common import SourceType
>>>>>>> 0069c1b (adding docstrings)
from ensembl.vep.models import Feature
from ensembl.vep.sources.feature_source import FeatureSource


@define(kw_only=True, auto_attribs=True)
<<<<<<< HEAD
class ConcreteFeatureSource(FeatureSource):
    """Concrete feature source used to test interface helpers."""

    closed: bool = field(default=False)
    sequence_source: SequenceSource | None = field(default=None)  # ignore[assignment]
=======
class _MinimalFeatureSource(FeatureSource):
    """Smallest possible concrete FeatureSource, for exercising the abstract base."""

    source_type: SourceType = SourceType.FILE
    closed: bool = field(default=False)
>>>>>>> 0069c1b (adding docstrings)

    def get_all_features(self) -> list[Feature]:
        return []

    def iter_all_features(self) -> Iterator[Feature]:
        return iter([])

    def get_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> list[Feature]:
        return []

    def iter_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> Iterator[Feature]:
        return iter([])

    def close(self) -> None:
        self.closed = True


<<<<<<< HEAD
class TestFeatureSourceInitialization:
    """Test FeatureSource initialization."""

    def test_abstract_feature_source_cannot_be_instantiated(self) -> None:
        """Test FeatureSource cannot be instantiated directly."""

        with pytest.raises(TypeError, match="abstract class"):
            FeatureSource(  # type: ignore[abstract]  # pylint: disable=abstract-class-instantiated
                source_type=SourceType.FILE, sequence_source=None
            )

    def test_concrete_feature_source_sets_source_type(self) -> None:
        """Test concrete feature source sets the source_type correctly."""

        feature_source = ConcreteFeatureSource(source_type=SourceType.FILE)

        assert feature_source.source_type == SourceType.FILE

    def test_concrete_feature_source_has_default_slice_cache(self) -> None:
        """Test feature sources receive a slice cache by default."""

        feature_source = ConcreteFeatureSource(source_type=SourceType.FILE)

        assert isinstance(feature_source.slice_cache, SliceCache)

    def test_concrete_feature_source_accepts_shared_slice_cache(self) -> None:
        """Test a slice cache can be injected and shared between feature sources."""

        slice_cache = WeakSliceCache()
        first = ConcreteFeatureSource(source_type=SourceType.FILE, slice_cache=slice_cache)
        second = ConcreteFeatureSource(source_type=SourceType.FILE, slice_cache=slice_cache)

        assert first.slice_cache is slice_cache
        assert second.slice_cache is slice_cache


class TestFeatureSourceContextManager:
    """Test for FeatureSource context manager."""

    @pytest.fixture
    def source(self) -> ConcreteFeatureSource:
        """Create a minimal FeatureSource instance for testing."""

        return ConcreteFeatureSource(source_type=SourceType.FILE)

    def test_enter_returns_self(self, source: ConcreteFeatureSource) -> None:
        """Test __enter__ returns the source instance itself."""

        with source as s:
            assert s is source

    def test_exit_calls_close(self, source: ConcreteFeatureSource) -> None:
        """Test __exit__ calls close on normal exit."""

        with source:
            pass

        assert source.closed is True
=======
class TestFeatureSourceAbstractness:
    """Test that FeatureSource cannot be instantiated directly."""

    def test_cannot_instantiate_directly(self) -> None:
        """Test FeatureSource raises TypeError when instantiated directly."""
        with pytest.raises(TypeError, match="abstract"):
            FeatureSource()  # type: ignore[abstract]

    def test_concrete_subclass_can_be_instantiated(self) -> None:
        """Test a subclass implementing all abstract methods can be instantiated."""
        source = _MinimalFeatureSource()
        assert isinstance(source, FeatureSource)


class TestFeatureSourceContextManager:
    """Test the default __enter__/__exit__ implementation."""

    @pytest.fixture
    def source(self) -> _MinimalFeatureSource:
        """Create a minimal FeatureSource instance for testing."""
        return _MinimalFeatureSource()

    def test_enter_returns_self(self, source: _MinimalFeatureSource) -> None:
        """Test __enter__ returns the source instance itself."""
        with source as entered:
            assert entered is source

    def test_exit_calls_close(self, source: _MinimalFeatureSource) -> None:
        """Test __exit__ calls close on normal exit."""
        with source:
            pass
        assert source.closed is True

    def test_close_called_even_on_exception(self, source: _MinimalFeatureSource) -> None:
        """Test __exit__ calls close even when the with-block raises."""
        with pytest.raises(ValueError):
            with source:
                raise ValueError("boom")
        assert source.closed is True

    def test_close_not_called_without_context_manager(self, source: _MinimalFeatureSource) -> None:
        """Test close is not called just by constructing the source."""
        assert source.closed is False
>>>>>>> 0069c1b (adding docstrings)
