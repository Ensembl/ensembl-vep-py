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
"""Unit tests for sequence sources."""

from pathlib import Path

from attrs import define, field
import pytest

from ensembl.vep.caches import SequenceChunkCache
from ensembl.vep.common import Strand
from ensembl.vep.sources.file.fasta.sequence_source import FastaSequenceSource
from ensembl.vep.sources.sequence_source import SequenceSource


@define(kw_only=True, auto_attribs=True)
class ConcreteSequenceSource(SequenceSource):
    """Concrete sequence source used to test shared helpers."""

    data: dict[str, str] = field(factory=dict)
    calls: list[tuple[str, int, int]] = field(factory=list)

    def close(self) -> None:
        return None

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


class TestSequenceSource:
    """Tests for the shared SequenceSource implementation."""

    def test_sequence_source_is_abstract(self) -> None:
        with pytest.raises(TypeError, match="abstract class"):
            SequenceSource()  # type: ignore[abstract,call-arg]

    def test_get_sequence_pads_missing_bases(self) -> None:
        source = ConcreteSequenceSource(data={"chr1": "ACGT"})

        assert source.get_sequence_by_start_end_strand("chr1", 3, 6) == "GTNN"

    def test_get_sequence_reverse_complements(self) -> None:
        source = ConcreteSequenceSource(data={"chr1": "AACG"})

        assert source.get_sequence_by_start_end_strand("chr1", 1, 4, strand=Strand.REVERSE) == "CGTT"

    def test_get_sequence_returns_ns_for_missing_region(self) -> None:
        source = ConcreteSequenceSource(data={})

        assert source.get_sequence_by_start_end_strand("chr1", 1, 4) == "NNNN"

    def test_length_is_cached(self) -> None:
        source = ConcreteSequenceSource(data={"chr1": "ACGT"})

        assert source.get_seq_region_length("chr1") == 4
        assert source.get_seq_region_length("chr1") == 4
        assert source.calls == []

    def test_chunk_cache_reuses_chunks(self) -> None:
        source = ConcreteSequenceSource(data={"chr1": "ABCDEFGHIJKLMNOP"}, cache=SequenceChunkCache(chunk_power=2))

        assert source.get_sequence_by_start_end_strand("chr1", 2, 6) == "BCDEF"
        assert len(source.cache._cache) == 2  # pylint: disable=protected-access

        assert source.get_sequence_by_start_end_strand("chr1", 3, 7) == "CDEFG"
        assert len(source.cache._cache) == 2  # pylint: disable=protected-access

    def test_cache_can_be_disabled(self) -> None:
        source = ConcreteSequenceSource(
            data={"chr1": "ABCDEFGHIJKLMNOP"},
            cache=SequenceChunkCache(enabled=False, chunk_power=2),
        )

        assert source.get_sequence_by_start_end_strand("chr1", 2, 6) == "BCDEF"
        assert source.get_sequence_by_start_end_strand("chr1", 2, 6) == "BCDEF"
        assert len(source.calls) == 4


class TestFastaSequenceSource:
    """Tests for FASTA-backed sequence loading."""

    def test_loads_sequences_from_fasta(self, tmp_path: Path) -> None:
        fasta_path = tmp_path / "sample.fa"
        fasta_path.write_text(">chr1\nACGTTGCA\n>chr2\nAAAA\n", encoding="utf-8")

        source = FastaSequenceSource(file_path=fasta_path)

        assert source.get_sequence_by_start_end_strand("chr1", 2, 5) == "CGTT"
        assert source.get_sequence_by_start_end_strand("chr1", 2, 5, strand=Strand.REVERSE) == "AACG"
        assert source.get_sequence_by_start_end_strand("chr2", 1, 4) == "AAAA"
        assert source.get_sequence_by_start_end_strand("chr3", 1, 4) == "NNNN"

    def test_length_is_resolved_from_fasta(self, tmp_path: Path) -> None:
        fasta_path = tmp_path / "sample.fa"
        fasta_path.write_text(">chr1\nACGTTGCA\n", encoding="utf-8")

        source = FastaSequenceSource(file_path=fasta_path)
        assert source.get_seq_region_length("chr1") == 8

    def test_open_resets_loaded_state(self, tmp_path: Path) -> None:
        fasta_path = tmp_path / "sample.fa"
        fasta_path.write_text(">chr1\nACGT\n", encoding="utf-8")

        source = FastaSequenceSource(file_path=fasta_path)
        assert source.get_sequence_by_start_end_strand("chr1", 1, 4) == "ACGT"

        source.close()
        assert source.cache is None
        assert not source._seq_region_lengths  # pylint: disable=protected-access
