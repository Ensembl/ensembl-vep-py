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
"""Unit tests for the FixedSizeInputBuffer class."""

from collections import deque
from typing import Iterator
import pytest

from ensembl.vep.buffers.fixed_size_input_buffer import FixedSizeInputBuffer
from ensembl.vep.common import SourceType
from ensembl.vep.models.variation.variation import Variation
from tests.buffers.conftest import DummyVariationSource


class TestBufferInitialization:
    """Test Buffer initialization."""

    def test_initialization_without_keywords(self, variation_list: list[Variation]) -> None:
        """Test FixedSizeInputBuffer initialization without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            var_source = DummyVariationSource(source_type=SourceType.FILE, features=variation_list)
            _ = FixedSizeInputBuffer(var_source)  # type: ignore[misc]

    def test_fixed_size_input_buffer_initial_attributes(self, variation_list: list[Variation]) -> None:
        """Test FixedSizeInputBuffer initialises attributes correctly."""

        var_source = DummyVariationSource(source_type=SourceType.FILE, features=variation_list)
        buffer = FixedSizeInputBuffer(source=var_source)

        assert buffer.source is var_source
        assert buffer.pre_buffer == deque([])
        assert buffer.buffer == []
        assert not buffer.is_empty


class TestBufferMethods:
    """Test Buffer methods."""

    def test_fixed_size_input_buffer_empty_source(self) -> None:
        """Test fixed_size_input_buffer ``next`` method."""

        var_source = DummyVariationSource(source_type=SourceType.FILE, features=[])
        buffer = FixedSizeInputBuffer(source=var_source)

        assert isinstance(buffer.feature_iterator, Iterator)
        buffer.next()
        assert not buffer.buffer
        assert buffer.is_empty

    def test_fixed_size_input_buffer_next(self, variation_list: list[Variation]) -> None:
        """Test fixed_size_input_buffer ``next`` method and ``is_empty`` property."""

        var_source = DummyVariationSource(source_type=SourceType.FILE, features=variation_list)
        buffer = FixedSizeInputBuffer(source=var_source)

        assert isinstance(buffer.feature_iterator, Iterator)
        buffer.next()
        assert len(buffer.buffer) == 4
        assert buffer.buffer[0].start == 0
        assert buffer.buffer[0].slice.seq_region_name == "1"
        assert len(buffer.pre_buffer) == 1
        assert buffer.pre_buffer[0].start == 4
        assert buffer.pre_buffer[0].slice.seq_region_name == "2"
        buffer.next()
        assert len(buffer.buffer) == 4
        assert buffer.buffer[0].start == 4
        assert buffer.buffer[0].slice.seq_region_name == "2"
        assert len(buffer.pre_buffer) == 1
        assert buffer.pre_buffer[0].start == 8
        assert buffer.pre_buffer[0].slice.seq_region_name == "3"
        assert not buffer.is_empty
        buffer.next()
        assert len(buffer.buffer) == 3
        assert buffer.buffer[0].start == 8
        assert buffer.buffer[0].slice.seq_region_name == "3"
        assert buffer.pre_buffer == deque([])
        assert buffer.is_empty
