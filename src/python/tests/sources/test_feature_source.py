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

import pytest
from attrs import define, field

from ensembl.vep.common import SourceType
from ensembl.vep.models import Feature
from ensembl.vep.sources.feature_source import FeatureSource


@define(kw_only=True, auto_attribs=True)
class _MinimalFeatureSource(FeatureSource):
    """Smallest possible concrete FeatureSource, for exercising the abstract base."""

    source_type: SourceType = SourceType.FILE
    closed: bool = field(default=False)

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
