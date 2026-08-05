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
"""Unit tests for CircularSlice class."""

from copy import deepcopy
import pytest

from ensembl.vep.common import CoordSystemLevel, Strand
from ensembl.vep.models.sequence import CircularSlice, CoordSystem
from tests.models.conftest import DummySequenceSource


@pytest.fixture
def circular_slice() -> CircularSlice:
    """Create a standard circular slice for model tests."""
    return CircularSlice(
        coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
        seq_region_name="1",
        start=10,
        end=200,
        strand=Strand.FORWARD,
    )


class TestCircularSliceInitialization:
    """Test CircularSlice initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test CircularSlice initialization without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = CircularSlice(
                CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),  # type: ignore[misc]
                seq_region_name="1",
                start=10,
                end=200,
                strand=Strand.FORWARD,
            )


class TestCircularSliceProperties:
    """Test property accessors."""

    def test_core_properties(self, circular_slice: CircularSlice) -> None:
        """Test core properties are set correctly."""
        assert circular_slice.coord_system == CoordSystem(
            name="chromosome", level=CoordSystemLevel.CHROMOSOME
        )
        assert circular_slice.seq_region_name == "1"
        assert circular_slice.seq_region_length is None
        assert circular_slice.start == 10
        assert circular_slice.end == 200
        assert circular_slice.strand == Strand.FORWARD
        assert circular_slice.seq is None

    def test_seq_region_length_is_lazy_when_missing(self, dummy_sequence_source: DummySequenceSource) -> None:
        """Test seq_region_length can be resolved from the source lazily."""

        source = dummy_sequence_source
        source.data = {"1": "ACGTACGT"}
        genomic_slice = CircularSlice(
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
        genomic_slice = CircularSlice(
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
        genomic_slice = CircularSlice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            start=1,
            end=4,
            strand=Strand.REVERSE,
            sequence_source=source,
        )

        assert genomic_slice.seq == "CGTT"

    def test_circular_slice_spanning_origin_concatenates_segments(
        self, dummy_sequence_source: DummySequenceSource
    ) -> None:
        """Test origin-spanning slices concatenate the source segments."""
        source = dummy_sequence_source
        source.data = {"1": "ACGTTGCA"}
        genomic_slice = CircularSlice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            start=7,
            end=10,
            strand=Strand.FORWARD,
            sequence_source=source,
        )

        assert genomic_slice.seq == "CAAC"
        assert source.calls == [("1", 7, 8), ("1", 1, 2)]

    def test_circular_slice_reverse_complements_after_concatenation(
        self, dummy_sequence_source: DummySequenceSource
    ) -> None:
        """Test reverse complementation happens after segment concatenation."""
        source = dummy_sequence_source
        source.data = {"1": "ACGTTGCA"}
        genomic_slice = CircularSlice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            start=7,
            end=10,
            strand=Strand.REVERSE,
            sequence_source=source,
        )

        assert genomic_slice.seq == "GTTG"

    def test_circular_slice_fills_missing_bases_with_n(
        self, dummy_sequence_source: DummySequenceSource
    ) -> None:
        """Test requests longer than the region are padded with Ns."""
        source = dummy_sequence_source
        source.data = {"1": "ACGT"}
        genomic_slice = CircularSlice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            start=7,
            end=12,
            strand=Strand.FORWARD,
            sequence_source=source,
        )

        assert genomic_slice.seq == "GTACNN"


class TestCircularSliceEquality:
    """Test equality comparison."""

    def test_equal_slices(self, circular_slice: CircularSlice) -> None:
        """Test two identical CircularSlices are equal."""
        slice2 = deepcopy(circular_slice)
        assert circular_slice == slice2

    def test_different_coord_system(self, circular_slice: CircularSlice) -> None:
        """Test two CircularSlices with different coord_system are not equal."""
        slice2 = deepcopy(circular_slice)
        slice2.coord_system = CoordSystem(name="contig", level=CoordSystemLevel.CONTIG)
        assert circular_slice != slice2

    def test_different_seq_region_name(self, circular_slice: CircularSlice) -> None:
        """Test two CircularSlices with different seq_region_name are not equal."""
        slice2 = deepcopy(circular_slice)
        slice2.seq_region_name = "2"
        assert circular_slice != slice2

    def test_different_start(self, circular_slice: CircularSlice) -> None:
        """Test two CircularSlices with different start are not equal."""
        slice2 = deepcopy(circular_slice)
        slice2.start = 20
        assert circular_slice != slice2

    def test_different_end(self, circular_slice: CircularSlice) -> None:
        """Test two CircularSlices with different end are not equal."""
        slice2 = deepcopy(circular_slice)
        slice2.end = 300
        assert circular_slice != slice2

    def test_different_strand(self, circular_slice: CircularSlice) -> None:
        """Test two CircularSlices with different strand are not equal."""
        slice2 = deepcopy(circular_slice)
        slice2.strand = Strand.REVERSE
        assert circular_slice != slice2

    def test_diferent_seq_region_length(
        self, circular_slice: CircularSlice, dummy_sequence_source: DummySequenceSource
    ) -> None:
        """Test two CircularSlices with different seq_region_length are still equal."""
        slice2 = deepcopy(circular_slice)
        dummy_sequence_source.data = {"1": "ACGTACGTACGT"}
        slice2.sequence_source = dummy_sequence_source
        assert circular_slice == slice2

    def test_different_type(self, circular_slice: CircularSlice) -> None:
        """Test comparing CircularSlice to different type returns False."""

        class AnotherCircularSlice(CircularSlice):
            """CircularSlice subclass for testing."""

        slice2 = AnotherCircularSlice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            start=10,
            end=200,
            strand=Strand.FORWARD,
        )
        assert circular_slice != slice2
