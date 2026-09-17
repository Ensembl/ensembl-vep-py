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

from ensembl.vep.buffers.buffer import Config
from ensembl.vep.buffers.fixed_size_input_buffer import FixedSizeInputBuffer
from ensembl.vep.common import SourceType
from ensembl.vep.models.variation.variation import Variation
from tests.buffers.conftest import DummyVariationSource


@pytest.fixture(name="var_source")
def var_source_fixture(variation_list: list[Variation]) -> DummyVariationSource:
    """Create a dummy VariationSource."""
    return DummyVariationSource(source_type=SourceType.FILE, features=variation_list)


@pytest.fixture(name="fixed_size_buffer")
def fixed_size_buffer_fixture(config: Config, var_source: DummyVariationSource) -> FixedSizeInputBuffer:
    """Create a FixedSizeInputBuffer object."""
    return FixedSizeInputBuffer(config=config, source=var_source)


@pytest.fixture(name="unordered_var_source")
def unordered_var_source_fixture(variation_list_out_of_order: list[Variation]) -> DummyVariationSource:
    """Create a dummy VariationSource with out of order variants."""
    return DummyVariationSource(source_type=SourceType.FILE, features=variation_list_out_of_order)


@pytest.fixture(name="unordered_buffer")
def unordered_buffer_fixture(
    config: Config, unordered_var_source: DummyVariationSource
) -> FixedSizeInputBuffer:
    """Create a FixedSizeInputBuffer object with unordered source."""
    return FixedSizeInputBuffer(config=config, source=unordered_var_source)


class TestBufferInitialization:
    """Test Buffer initialization."""

    def test_initialization_without_keywords(self, variation_list: list[Variation]) -> None:
        """Test FixedSizeInputBuffer initialization without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            var_source = DummyVariationSource(source_type=SourceType.FILE, features=variation_list)
            _ = FixedSizeInputBuffer(var_source)  # type: ignore[call-arg, arg-type, misc]

    def test_fixed_size_input_buffer_initial_attributes(
        self, fixed_size_buffer: FixedSizeInputBuffer, var_source: DummyVariationSource
    ) -> None:
        """Test FixedSizeInputBuffer initialises attributes correctly."""

        assert fixed_size_buffer.source is var_source
        assert fixed_size_buffer._pre_buffer == deque([])  # pylint: disable=protected-access
        assert fixed_size_buffer.buffer == []
        assert not fixed_size_buffer.is_empty


class TestBufferNext:
    """Test Buffer ``next`` method."""

    def test_fixed_size_input_buffer_empty_source(self) -> None:
        """Test fixed_size_input_buffer ``next`` method."""
        var_source = DummyVariationSource(source_type=SourceType.FILE, features=[])
        buffer = FixedSizeInputBuffer(source=var_source)

        assert isinstance(buffer.variation_iterator, Iterator)
        buffer.next()
        assert not buffer.buffer
        assert buffer.is_empty

    def test_fixed_size_input_buffer_next(self, fixed_size_buffer: FixedSizeInputBuffer) -> None:
        """Test fixed_size_input_buffer ``next`` method and ``is_empty`` property."""

        assert isinstance(fixed_size_buffer.variation_iterator, Iterator)
        fixed_size_buffer.next()
        assert len(fixed_size_buffer.buffer) == 5
        assert fixed_size_buffer.buffer[0].start == 1
        assert fixed_size_buffer.buffer[0].slice is not None
        assert fixed_size_buffer.buffer[0].slice.seq_region_name == "1"
        assert len(fixed_size_buffer._pre_buffer) == 1  # pylint: disable=protected-access
        assert fixed_size_buffer._pre_buffer[0].start == 6  # pylint: disable=protected-access
        assert fixed_size_buffer._pre_buffer[0].slice is not None  # pylint: disable=protected-access
        assert (
            fixed_size_buffer._pre_buffer[0].slice.seq_region_name == "1"  # pylint: disable=protected-access
        )
        fixed_size_buffer.next()
        assert len(fixed_size_buffer.buffer) == 3
        assert fixed_size_buffer.buffer[0].start == 6
        assert fixed_size_buffer.buffer[0].slice is not None
        assert fixed_size_buffer.buffer[0].slice.seq_region_name == "1"
        assert len(fixed_size_buffer._pre_buffer) == 1  # pylint: disable=protected-access
        assert fixed_size_buffer._pre_buffer[0].start == 10  # pylint: disable=protected-access
        assert fixed_size_buffer._pre_buffer[0].slice is not None  # pylint: disable=protected-access
        assert (
            fixed_size_buffer._pre_buffer[0].slice.seq_region_name == "2"  # pylint: disable=protected-access
        )
        assert not fixed_size_buffer.is_empty
        fixed_size_buffer.next()
        assert len(fixed_size_buffer.buffer) == 4
        assert fixed_size_buffer.buffer[0].start == 10
        assert fixed_size_buffer.buffer[0].slice is not None
        assert fixed_size_buffer.buffer[0].slice.seq_region_name == "2"
        assert fixed_size_buffer._pre_buffer == deque([])  # pylint: disable=protected-access
        assert fixed_size_buffer.is_empty

    def test_fixed_size_input_buffer_large_pre_buffer(
        self, fixed_size_buffer: FixedSizeInputBuffer, variation_list: list[Variation]
    ) -> None:
        """Test fixed_size_input_buffer ``next`` method stops if buffer is filled by pre_buffer."""
        fixed_size_buffer._pre_buffer = deque(variation_list)  # pylint: disable=protected-access

        fixed_size_buffer.next()
        assert len(fixed_size_buffer.buffer) == 5
        assert len(fixed_size_buffer._pre_buffer) == 7  # pylint: disable=protected-access

    def test_fixed_size_input_buffer_large_pre_buffer_diff_chr(
        self, fixed_size_buffer: FixedSizeInputBuffer, variation_list: list[Variation]
    ) -> None:
        """Test fixed_size_input_buffer ``next`` method stops if
        pre_buffer contains variants on different chromosomes."""
        fixed_size_buffer._pre_buffer = deque(variation_list[0:3])  # pylint: disable=protected-access
        fixed_size_buffer._pre_buffer.extend(variation_list[8:12])  # pylint: disable=protected-access

        fixed_size_buffer.next()
        assert len(fixed_size_buffer.buffer) == 3
        assert len(fixed_size_buffer._pre_buffer) == 4  # pylint: disable=protected-access

    def test_fixed_size_input_buffer_rejects_unsorted_variations(
        self, unordered_buffer: FixedSizeInputBuffer
    ) -> None:
        """Test that ``next`` rejects variations that are out of order."""
        with pytest.raises(RuntimeError, match="The input file appears to be unsorted"):
            unordered_buffer.next()

    def test_fixed_size_input_buffer_max_not_ordered_variants_threshold(
        self, unordered_buffer: FixedSizeInputBuffer
    ) -> None:
        """Test that ``next`` doesn't rejects unsorted variations if the thresholds are high enough."""
        unordered_buffer.config.max_not_ordered_variants = 5
        unordered_buffer.next()

    def test_fixed_size_input_buffer_max_not_ordered_variants_distance_threshold(
        self, unordered_buffer: FixedSizeInputBuffer
    ) -> None:
        """Test that ``next`` doesn't rejects unsorted variations if the thresholds are high enough."""
        unordered_buffer.config.max_not_ordered_variants_distance = 100
        unordered_buffer.next()


class TestBufferGetOverlappingVariation:
    """Test Buffer ``get_overlapping_variation`` method."""

    def test_empty_buffer_returns_none(self, fixed_size_buffer: FixedSizeInputBuffer) -> None:
        """Test that an empty buffer returns None."""
        assert fixed_size_buffer.get_overlapping_variation(start=1, end=10) is None

    def test_overlap_forward_strand(
        self, fixed_size_buffer: FixedSizeInputBuffer, variation_list: list[Variation]
    ) -> None:
        """Test that ``get_overlapping_variation`` returns the correct Variation objects."""
        assert fixed_size_buffer.get_overlapping_variation(start=1, end=10) is None

        fixed_size_buffer.next()
        assert fixed_size_buffer.get_overlapping_variation(start=2, end=3) == variation_list[1:3]
        assert fixed_size_buffer.get_overlapping_variation(start=6, end=6) == []

    def test_overlap_reverse_strand(
        self, fixed_size_buffer: FixedSizeInputBuffer, variation_list: list[Variation]
    ) -> None:
        """Test that ``get_overlapping_variation`` returns the correct Variation objects."""
        assert fixed_size_buffer.get_overlapping_variation(start=1, end=10) is None

        fixed_size_buffer.next()
        # start > end
        assert fixed_size_buffer.get_overlapping_variation(start=3, end=2) == variation_list[1:3]
