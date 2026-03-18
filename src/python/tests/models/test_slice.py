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
"""Unit tests for BaseFeature class."""

import pytest
from copy import deepcopy

from ensembl.vep.models.slice import Slice


class TestSliceInitialization:
    """Test Slice initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test Slice initialization without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = Slice(
                "chromosome", 
                seq_region_name="1", 
                seq_region_length=1000, 
                start=10, 
                end=200, 
                strand=1,
                seq=None
            )

class TestSliceProperties:
    """Test property accessors."""

    @pytest.fixture
    def slice(self) -> Slice:
        """Create a Slice instance for testing."""
        return Slice(
            coord_system="chromosome", 
            seq_region_name="1", 
            seq_region_length=1000, 
            start=10, 
            end=200, 
            strand=1,
            seq=None
        )
    
    def test_core_properties(self, slice: Slice) -> None:
        """Test core properties are set correctly."""
        assert slice.coord_system == "chromosome"
        assert slice.seq_region_name == "1"
        assert slice.seq_region_length == 1000
        assert slice.start == 10
        assert slice.end == 200
        assert slice.strand == 1
        assert slice.seq is None

class TestSliceEquality:
    """Test equality comparison."""

    @pytest.fixture
    def slice(self) -> Slice:
        """Create first Slice instance for testing."""
        return Slice(
            coord_system="chromosome", 
            seq_region_name="1", 
            seq_region_length=1000, 
            start=10, 
            end=200, 
            strand=1,
            seq=None
        )

    def test_equal_slices(self, slice: Slice) -> None:
        """Test two identical Slices are equal."""
        slice2 = deepcopy(slice)
        assert slice == slice2

    def test_different_coord_system(self, slice: Slice) -> None:
        """Test two Slices with different coord_system are not equal."""
        slice2 = deepcopy(slice)
        slice2.coord_system = "contig"
        assert slice != slice2

    def test_different_seq_region_name(self, slice: Slice) -> None:
        """Test two Slices with different seq_region_name are not equal."""
        slice2 = deepcopy(slice)
        slice2.seq_region_name = "2"
        assert slice != slice2

    def test_different_start(self, slice: Slice) -> None:
        """Test two Slices with different start are not equal."""
        slice2 = deepcopy(slice)
        slice2.start = 20
        assert slice != slice2

    def test_different_end(self, slice: Slice) -> None:
        """Test two Slices with different end are not equal."""
        slice2 = deepcopy(slice)
        slice2.end = 300
        assert slice != slice2

    def test_different_strand(self, slice: Slice) -> None:
        """Test two Slices with different strand are not equal."""
        slice2 = deepcopy(slice)
        slice2.strand = -1
        assert slice != slice2

    def test_diferent_seq_region_length(self, slice: Slice) -> None:
        """Test two Slices with different seq_region_length are still equal."""
        slice2 = deepcopy(slice)
        slice2.seq_region_length = 2000
        assert slice == slice2

    def test_different_type(self, slice: Slice) -> None:
        """Test comparing Slice to different type returns False."""
        class AnotherSlice(Slice):
            pass
        slice2 = AnotherSlice(
            coord_system="chromosome", 
            seq_region_name="1", 
            seq_region_length=1000, 
            start=10, 
            end=200, 
            strand=1,
            seq=None
        )
        assert slice != slice2