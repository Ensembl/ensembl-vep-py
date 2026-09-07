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
"""Unit tests for the buffer interface."""

import pytest
from typing import Iterator
from attrs import define

from ensembl.vep.buffers.buffer import Buffer
from ensembl.vep.common import Strand, SourceType
from ensembl.vep.models import Feature
from ensembl.vep.sources.feature_source import FeatureSource

def create_feature_objects(length: int) -> list[Feature]:
    return([Feature(start=value, end=value, strand=Strand.FORWARD) for value in range(length)])

_FEATURES: list[Feature] = create_feature_objects(11)

@define(kw_only=True, auto_attribs=True)
class DummyFeatureSource(FeatureSource):
    """Dummy FeatureSource for the Concrete Buffer"""

    def get_all_features(self) -> list[Feature]:
        """Return all features from this source."""

        return _FEATURES

    def iter_all_features(self) -> Iterator[Feature]:
        """Lazily yield all features from this source."""

        return iter(_FEATURES)

    def get_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> list[Feature]:
        """Return a subset of the features."""

        return _FEATURES[0:3]

    def iter_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> Iterator[Feature]:
        """Lazily yield a subset of the features."""

        return iter(_FEATURES[5:8])

    def close(self) -> None:
        """Empty _FEATURES list."""

        _FEATURES = []


@define(kw_only=True, auto_attribs=True)
class ConcreteBuffer(Buffer):
    """Concrete buffer used to test Buffer interface helpers."""
    buffer_size: int = 5

    def next(self) -> None:
        """Fill the buffer."""
        self.buffer = []
        self.buffer.extend(self.pre_buffer)
        self.pre_buffer = []
        for feat in self.feature_iterator:
            if len(self.buffer) == self.buffer_size:
                self.pre_buffer.append(feat)
                break
            else:
                self.buffer.append(feat)
        if len(self.buffer) < self.buffer_size:
            self._is_empty = True


class TestBufferInitialization:
    """Test Buffer initialization."""

    def test_abstract_buffer_cannot_be_instantiated(self) -> None:
        """Test Buffer cannot be instantiated directly."""

        with pytest.raises(TypeError, match="abstract class"):
            _ = Buffer()  # type: ignore[abstract]  # pylint: disable=abstract-class-instantiated

    def test_concrete_buffer_initial_attributes(self) -> None:
        """Test concrete buffer stores the buffer format."""

        feat_source = DummyFeatureSource(source_type=SourceType.FILE)
        buffer = ConcreteBuffer(source=feat_source)

        assert buffer.source is feat_source
        assert buffer.pre_buffer == []
        assert buffer.buffer == []
        assert not buffer.is_empty


class TestBufferProperties:
    """Test Buffer properties method."""

    def test_concrete_buffer_properties(self) -> None:
        """Test concrete buffer stores the buffer format."""

        feat_source = DummyFeatureSource(source_type=SourceType.FILE)
        buffer = ConcreteBuffer(source=feat_source)

        assert isinstance(buffer.feature_iterator, Iterator)
        buffer.next()
        assert [feat.start for feat in buffer.buffer] == [0, 1, 2, 3, 4]
        buffer.next()
        assert [feat.start for feat in buffer.buffer] == [5, 6, 7, 8, 9]
        assert not buffer.is_empty
        buffer.next()
        assert [feat.start for feat in buffer.buffer] == [10]
        assert buffer.is_empty
