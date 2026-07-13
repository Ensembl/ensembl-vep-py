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
from attrs import define, field
import pytest

from ensembl.vep.common import CoordSystemLevel, Strand
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.sources.sequence_source import SequenceSource


@define(kw_only=True, auto_attribs=True)
class DummySequenceSource(SequenceSource):
    """Minimal sequence source used for Slice tests."""

    data: dict[str, str] = field(factory=dict)
    calls: list[tuple[str, int, int]] = field(factory=list)
    closed: bool = field(default=False)

    def close(self) -> None:
        self.closed = True

    def _fetch_sequence_chunk(self, seq_region_name: str, start: int, end: int) -> str:
        self.calls.append((seq_region_name, start, end))
        sequence = self.data.get(seq_region_name)
        if sequence is None:
            return "N" * (end - start + 1)

        raw_start = start - 1
        raw_end = end
        chunk = sequence[raw_start:raw_end]
        if len(chunk) < (end - start + 1):
            chunk += "N" * ((end - start + 1) - len(chunk))
        return chunk

    def _get_seq_region_length(self, seq_region_name: str) -> int | None:
        sequence = self.data.get(seq_region_name)
        return len(sequence) if sequence is not None else None


class TestSliceInitialization:
    """Test Slice initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test Slice initialization without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = Slice(
                CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),  # type: ignore[misc]
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
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            seq_region_length=1000,
            start=10,
            end=200,
            strand=Strand.FORWARD,
            seq=None,
        )

    def test_core_properties(self, genomic_slice: Slice) -> None:
        """Test core properties are set correctly."""
        assert genomic_slice.coord_system == CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME)
        assert genomic_slice.seq_region_name == "1"
        assert genomic_slice.seq_region_length == 1000
        assert genomic_slice.start == 10
        assert genomic_slice.end == 200
        assert genomic_slice.strand == Strand.FORWARD
        assert genomic_slice.seq is None

    def test_seq_region_length_is_lazy_when_missing(self) -> None:
        """Test seq_region_length can be resolved from the source lazily."""

        source = DummySequenceSource(data={"1": "ACGTACGT"})
        genomic_slice = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            start=2,
            end=5,
            strand=Strand.FORWARD,
            sequence_source=source,
        )

        assert genomic_slice.seq_region_length == 8

    def test_lazy_sequence_loading_uses_source(self) -> None:
        source = DummySequenceSource(data={"1": "ACGTACGT"})
        genomic_slice = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            seq_region_length=8,
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

    def test_reverse_strand_returns_reverse_complement(self) -> None:
        source = DummySequenceSource(data={"1": "AACG"})
        genomic_slice = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            seq_region_length=4,
            start=1,
            end=4,
            strand=Strand.REVERSE,
            sequence_source=source,
        )

        assert genomic_slice.seq == "CGTT"


class TestSliceEquality:
    """Test equality comparison."""

    @pytest.fixture
    def genomic_slice(self) -> Slice:
        """Create first Slice instance for testing."""
        return Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
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
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            seq_region_length=1000,
            start=10,
            end=200,
            strand=Strand.FORWARD,
            seq=None,
        )
        assert genomic_slice != slice2
