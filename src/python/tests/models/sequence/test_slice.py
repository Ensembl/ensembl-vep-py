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
"""Unit tests for Slice class."""

from copy import deepcopy
import pytest

from ensembl.vep.common import CoordSystemLevel, Strand
from ensembl.vep.models.sequence import CoordSystem, Slice
from tests.models.conftest import DummySequenceSource


@pytest.fixture
def genomic_slice() -> Slice:
    """Create a standard linear slice for model tests."""
    return Slice(
        coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
        seq_region_name="1",
        start=10,
        end=200,
        strand=Strand.FORWARD,
    )


class TestSliceInitialization:
    """Test Slice initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test Slice initialization without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = Slice(
                CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),  # type: ignore[misc]
                seq_region_name="1",
                start=10,
                end=200,
                strand=Strand.FORWARD,
            )


class TestSliceProperties:
    """Test property accessors."""

    def test_core_properties(self, genomic_slice: Slice) -> None:
        """Test core properties are set correctly."""
        assert genomic_slice.coord_system == CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME)
        assert genomic_slice.seq_region_name == "1"
        assert genomic_slice.seq_region_length is None
        assert genomic_slice.start == 10
        assert genomic_slice.end == 200
        assert genomic_slice.strand == Strand.FORWARD
        assert genomic_slice.seq is None

    def test_seq_region_length_is_lazy_when_missing(self, dummy_sequence_source: DummySequenceSource) -> None:
        """Test seq_region_length can be resolved from the source lazily."""

        source = dummy_sequence_source
        source.data = {"1": "ACGTACGT"}
        genomic_slice = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            start=2,
            end=5,
            strand=Strand.FORWARD,
            sequence_source=source,
        )

        assert genomic_slice.seq_region_length == 8

    def test_lazy_sequence_loading_uses_source(self, dummy_sequence_source: DummySequenceSource) -> None:
        """Test sequence loading is lazy and cached."""
        source = dummy_sequence_source
        source.data = {"1": "ACGTACGT"}
        genomic_slice = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            start=2,
            end=5,
            strand=Strand.FORWARD,
            sequence_source=source,
        )

        assert source.calls == []
        assert genomic_slice.seq == "CGTA"
        assert source.calls == [("1", 2, 5)]

        # Test that subsequent calls to seq do not trigger additional source calls.
        assert genomic_slice.seq == "CGTA"
        assert source.calls == [("1", 2, 5)]

    def test_reverse_strand_returns_reverse_complement(
        self, dummy_sequence_source: DummySequenceSource
    ) -> None:
        """Test reverse-strand slices return reverse-complemented sequence."""
        source = dummy_sequence_source
        source.data = {"1": "AACG"}
        genomic_slice = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            start=1,
            end=4,
            strand=Strand.REVERSE,
            sequence_source=source,
        )

        assert genomic_slice.seq == "CGTT"

    def test_seq_can_be_set_unsafe(self, genomic_slice: Slice) -> None:
        """Test the cached sequence can be supplied through the unsafe context."""
        with pytest.raises(AttributeError, match="property 'seq' of 'Slice' object has no setter"):
            genomic_slice.seq = "ACGT"  # type: ignore[misc]

        with genomic_slice.unsafe() as unsafe_slice:
            setattr(unsafe_slice, "_seq", "ACGT")

        assert genomic_slice.seq == "ACGT"

    def test_seq_region_length_can_be_set_unsafe(self, genomic_slice: Slice) -> None:
        """Test the cached sequence-region length can be supplied through the unsafe context."""
        with genomic_slice.unsafe() as unsafe_slice:
            setattr(unsafe_slice, "_seq_region_length", 4)

        assert genomic_slice.seq_region_length == 4


class TestSliceEquality:
    """Test equality comparison."""

    def test_equal_slices(self, genomic_slice: Slice) -> None:
        """Test two identical Slices are equal."""
        slice2 = deepcopy(genomic_slice)
        assert genomic_slice == slice2

    def test_different_coord_system(self, genomic_slice: Slice) -> None:
        """Test two Slices with different coord_system are not equal."""
        slice2 = deepcopy(genomic_slice)
        slice2.coord_system = CoordSystem(name="contig", level=CoordSystemLevel.CONTIG)
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

    def test_different_seq_region_length(
        self, genomic_slice: Slice, dummy_sequence_source: DummySequenceSource
    ) -> None:
        """Test two Slices with different seq_region_length are still equal."""
        slice2 = deepcopy(genomic_slice)
        dummy_sequence_source.data = {"1": "ACGTACGTACGT"}
        slice2.sequence_source = dummy_sequence_source
        assert genomic_slice == slice2

    def test_different_type(self, genomic_slice: Slice) -> None:
        """Test comparing Slice to different type returns False."""

        class AnotherSlice(Slice):
            """Slice subclass for testing."""

        slice2 = AnotherSlice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            start=10,
            end=200,
            strand=Strand.FORWARD,
        )
        assert genomic_slice != slice2
