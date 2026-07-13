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
"""Unit tests for the sequence source interface."""

from collections.abc import Callable

from attrs import define, field
import pytest

from ensembl.vep.common import Strand
from ensembl.vep.sources.sequence_source import SequenceSource


@define(kw_only=True, auto_attribs=True)
class ConcreteSequenceSource(SequenceSource):
    """Concrete sequence source used to test shared helpers."""

    data: dict[str, str] = field(factory=dict)
    closed: bool = False
    fetch_calls: list[tuple[str, int, int]] = field(factory=list)
    length_calls: list[str] = field(factory=list)

    def close(self) -> None:
        """Close the source."""

        self.closed = True

    def _fetch_sequence_chunk(self, seq_region_name: str, start: int, end: int) -> str:
        """Return sequence from the in-memory fixture."""

        self.fetch_calls.append((seq_region_name, start, end))
        sequence = self.data.get(seq_region_name)
        if sequence is None:
            return "N" * (end - start + 1)

        return sequence[start - 1 : end]

    def _get_seq_region_length(self, seq_region_name: str) -> int | None:
        """Return the fixture sequence length."""

        self.length_calls.append(seq_region_name)
        sequence = self.data.get(seq_region_name)
        return len(sequence) if sequence is not None else None


class DirectSequenceCache:
    """Minimal cache that delegates spans directly to the source loader."""

    def __init__(self) -> None:
        """Initialize the cache spy."""

        self.calls: list[tuple[str, int, int]] = []
        self.cleared = False

    def get_sequence(
        self,
        seq_region_name: str,
        start: int,
        end: int,
        loader: Callable[[str, int, int], str | None],
    ) -> str | None:
        """Record the request and delegate loading to the source."""

        self.calls.append((seq_region_name, start, end))
        return loader(seq_region_name, start, end)

    def clear(self) -> None:
        """Record cache clearing."""

        self.cleared = True


class TestSequenceSourceInitialization:
    """Test SequenceSource initialization."""

    def test_abstract_sequence_source_cannot_be_instantiated(self) -> None:
        """Test SequenceSource cannot be instantiated directly."""

        with pytest.raises(TypeError, match="abstract class"):
            _ = SequenceSource()  # type: ignore[abstract,call-arg]

    def test_concrete_source_starts_with_empty_length_cache(self) -> None:
        """Test concrete source initializes derived caches."""

        source = ConcreteSequenceSource(data={"chr1": "ACGT"})

        assert not source._seq_region_lengths  # pylint: disable=protected-access


class TestSequenceSourceSequence:
    """Tests for SequenceSource.get_sequence_by_start_end_strand."""

    def test_get_sequence_returns_forward_sequence(self) -> None:
        """Test get_sequence_by_start_end_strand returns forward sequence."""

        source = ConcreteSequenceSource(data={"chr1": "ACGTTGCA"})

        assert source.get_sequence_by_start_end_strand("chr1", 2, 5) == "CGTT"
        assert source.fetch_calls == [("chr1", 2, 5)]

    def test_get_sequence_pads_missing_bases(self) -> None:
        """Test missing trailing bases are padded with Ns."""

        source = ConcreteSequenceSource(data={"chr1": "ACGT"})

        assert source.get_sequence_by_start_end_strand("chr1", 3, 6) == "GTNN"

    def test_get_sequence_returns_ns_for_missing_region(self) -> None:
        """Test unknown regions return Ns for the requested length."""

        source = ConcreteSequenceSource(data={})

        assert source.get_sequence_by_start_end_strand("chr1", 1, 4) == "NNNN"

    def test_get_sequence_reverse_complements(self) -> None:
        """Test reverse-strand requests are reverse-complemented."""

        source = ConcreteSequenceSource(data={"chr1": "AACG"})

        assert source.get_sequence_by_start_end_strand("chr1", 1, 4, strand=Strand.REVERSE) == "CGTT"

    @pytest.mark.parametrize(
        ("start", "end"),
        [
            (0, 1),
            (-1, 1),
            (3, 2),
        ],
    )
    def test_get_sequence_rejects_invalid_coordinates(self, start: int, end: int) -> None:
        """Test sequence requests must use valid 1-based spans."""

        source = ConcreteSequenceSource(data={"chr1": "ACGT"})

        with pytest.raises(ValueError, match="1-based"):
            source.get_sequence_by_start_end_strand("chr1", start, end)


class TestSequenceSourceCache:
    """Tests for SequenceSource cache delegation."""

    def test_get_sequence_delegates_to_configured_cache(self) -> None:
        """Test SequenceSource uses the generic cache interface."""

        cache = DirectSequenceCache()
        source = ConcreteSequenceSource(data={"chr1": "ABCDEFGHIJKLMNOP"}, cache=cache)

        assert source.get_sequence_by_start_end_strand("chr1", 2, 6) == "BCDEF"
        assert cache.calls == [("chr1", 2, 6)]
        assert source.fetch_calls == [("chr1", 2, 6)]

    def test_clear_cache_clears_length_cache_and_configured_cache(self) -> None:
        """Test clear_cache removes derived source and configured cache state."""

        cache = DirectSequenceCache()
        source = ConcreteSequenceSource(data={"chr1": "ACGT"}, cache=cache)

        assert source.get_seq_region_length("chr1") == 4
        source.clear_cache()

        assert not source._seq_region_lengths  # pylint: disable=protected-access
        assert cache.cleared is True


class TestSequenceSourceLength:
    """Tests for SequenceSource.get_seq_region_length."""

    def test_get_seq_region_length_returns_source_length(self) -> None:
        """Test sequence region length is read from the concrete source."""

        source = ConcreteSequenceSource(data={"chr1": "ACGT"})

        assert source.get_seq_region_length("chr1") == 4

    def test_get_seq_region_length_is_cached(self) -> None:
        """Test sequence region length is cached after first lookup."""

        source = ConcreteSequenceSource(data={"chr1": "ACGT"})

        assert source.get_seq_region_length("chr1") == 4
        assert source.get_seq_region_length("chr1") == 4
        assert source.length_calls == ["chr1"]

    def test_get_seq_region_length_returns_none_for_unknown_region(self) -> None:
        """Test unknown sequence region length returns None."""

        source = ConcreteSequenceSource(data={})

        assert source.get_seq_region_length("chr1") is None


class TestSequenceSourceContextManager:
    """Tests for SequenceSource context manager."""

    def test_enter_returns_self(self) -> None:
        """Test __enter__ returns the source instance itself."""

        source = ConcreteSequenceSource(data={"chr1": "ACGT"})

        with source as opened:
            assert opened is source

    def test_exit_calls_close(self) -> None:
        """Test __exit__ closes the source."""

        source = ConcreteSequenceSource(data={"chr1": "ACGT"})

        with source:
            pass

        assert source.closed is True
