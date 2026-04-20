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
"""Unit tests for BaseFeature class."""

from copy import deepcopy
import pytest

from ensembl.vep.common import Strand
from ensembl.vep.models.sequence import Slice


class TestSliceInitialization:
    """Test Slice initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test Slice initialization without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = Slice(
                coord_system="chromosome",
                seq_region_name="1",
                seq_region_length=1000,
                start=10,
                end=200,
                strand=Strand.FORWARD,
                seq=None,
            )


class TestSliceProperties:
    """Test property accessors."""

    @pytest.fixture
    def genomic_slice(self) -> Slice:
        """Create a Slice instance for testing."""
        return Slice(
            coord_system="chromosome",
            seq_region_name="1",
            seq_region_length=1000,
            start=10,
            end=200,
            strand=Strand.FORWARD,
            seq=None,
        )

    def test_core_properties(self, genomic_slice: Slice) -> None:
        """Test core properties are set correctly."""
        assert genomic_slice.coord_system == "chromosome"
        assert genomic_slice.seq_region_name == "1"
        assert genomic_slice.seq_region_length == 1000
        assert genomic_slice.start == 10
        assert genomic_slice.end == 200
        assert genomic_slice.strand == 1
        assert genomic_slice.seq is None


class TestSliceEquality:
    """Test equality comparison."""

    @pytest.fixture
    def genomic_slice(self) -> Slice:
        """Create first Slice instance for testing."""
        return Slice(
            coord_system="chromosome",
            seq_region_name="1",
            seq_region_length=1000,
            start=10,
            end=200,
            strand=Strand.FORWARD,
            seq=None,
        )

    def test_equal_slices(self, genomic_slice: Slice) -> None:
        """Test two identical Slices are equal."""
        slice2 = deepcopy(genomic_slice)
        assert genomic_slice == slice2

    def test_different_coord_system(self, genomic_slice: Slice) -> None:
        """Test two Slices with different coord_system are not equal."""
        slice2 = deepcopy(genomic_slice)
        slice2.coord_system = "contig"
        assert genomic_slice != slice2

    def test_different_seq_region_name(self, genomic_slice: Slice) -> None:
        """Test two Slices with different seq_region_name are not equal."""
        slice2 = deepcopy(genomic_slice)
        slice2.seq_region_name = "2"
        assert genomic_slice != slice2

    def test_different_start(self, genomic_slice: Slice) -> None:
        """Test two Slices with different start are not equal."""
        slice2 = deepcopy(genomic_slice)
        slice2.start = 20
        assert genomic_slice != slice2

    def test_different_end(self, genomic_slice: Slice) -> None:
        """Test two Slices with different end are not equal."""
        slice2 = deepcopy(genomic_slice)
        slice2.end = 300
        assert genomic_slice != slice2

    def test_different_strand(self, genomic_slice: Slice) -> None:
        """Test two Slices with different strand are not equal."""
        slice2 = deepcopy(genomic_slice)
        slice2.strand = Strand.REVERSE
        assert genomic_slice != slice2

    def test_diferent_seq_region_length(self, genomic_slice: Slice) -> None:
        """Test two Slices with different seq_region_length are still equal."""
        slice2 = deepcopy(genomic_slice)
        slice2.seq_region_length = 2000
        assert genomic_slice == slice2

    def test_different_type(self, genomic_slice: Slice) -> None:
        """Test comparing Slice to different type returns False."""

        class AnotherSlice(Slice):
            """Slice subclass for testing."""

        slice2 = AnotherSlice(
            coord_system="chromosome",
            seq_region_name="1",
            seq_region_length=1000,
            start=10,
            end=200,
            strand=Strand.FORWARD,
            seq=None,
        )
        assert genomic_slice != slice2
