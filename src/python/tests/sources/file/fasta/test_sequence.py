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
"""Unit tests for the FASTA sequence source."""

from pathlib import Path

import pytest

from ensembl.vep.common import Strand
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource

SAMPLE = ">chr1\n" "ACGTTGCA\n" ">chr2\n" "AAAA\n"
SAMPLE_FAI = "chr1\t8\t6\t8\t9\nchr2\t4\t21\t4\t5\n"


def _write_fasta(path: Path, data: str = SAMPLE) -> Path:
    """Write a FASTA fixture file."""

    path.write_text(data, encoding="utf-8")
    return path


def _write_indexed_fasta(path: Path, data: str = SAMPLE, index: str = SAMPLE_FAI) -> Path:
    """Write a FASTA fixture with an accompanying FAI index."""

    path.write_text(data, encoding="utf-8")
    path.with_suffix(path.suffix + ".fai").write_text(index, encoding="utf-8")
    return path


class MetadataOnlyParser:
    """Minimal parser stub that only exposes metadata."""

    def __init__(self, metadata: dict[str, object]) -> None:
        """Initialize the parser stub."""

        self._metadata = metadata
        self.reader_mode = "sequential"
        self.closed = False

    @property
    def metadata(self) -> dict[str, object]:
        """Return parser metadata."""

        return self._metadata

    def open(self, *_args: object, **_kwargs: object) -> None:
        """Open the parser stub."""

        return None

    def close(self) -> None:
        """Close the parser stub."""

        self.closed = True

    def seek(self, *_args: object, **_kwargs: object) -> bool:
        """Raise if sequence lookup tries to seek."""

        raise AssertionError("seek should not be called")

    def next(self) -> None:
        """Raise if sequence lookup tries to read one record."""

        raise AssertionError("next should not be called")

    def records(self, *_args: object, **_kwargs: object) -> list[object]:
        """Raise if length lookup tries to load records."""

        raise AssertionError("records should not be called")


class TestFastaSequenceSourceInitialization:
    """Test FastaSequenceSource initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test FastaSequenceSource initialization without keywords raises error."""

        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = FastaSequenceSource("sample.fa")  # type: ignore[misc]

    def test_initialization_sets_core_attributes(self) -> None:
        """Test FastaSequenceSource stores source and parser attributes."""

        source = FastaSequenceSource()

        assert source.file_path is None
        assert source.mode == "auto"
        assert source.adaptor.format == "fasta"  # pylint: disable=protected-access
        assert source.adaptor.reader_mode == "auto"  # pylint: disable=protected-access
        assert source._sequences is None  # pylint: disable=protected-access


class TestFastaSequenceSourceOpenClose:
    """Tests for FastaSequenceSource open and close."""

    def test_close_clears_loaded_sequences_and_lengths(self, tmp_path: Path) -> None:
        """Test close clears FASTA-specific and derived source caches."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = FastaSequenceSource(file_path=path, mode="sequential")

        assert source.get_sequence_by_start_end_strand("chr1", 1, 4) == "ACGT"
        assert source.get_seq_region_length("chr1") == 8
        source.close()

        assert source._sequences is None  # pylint: disable=protected-access
        assert not source._seq_region_lengths  # pylint: disable=protected-access

    def test_context_manager_closes_source(self, tmp_path: Path) -> None:
        """Test context manager returns and closes the source."""

        path = _write_fasta(tmp_path / "sample.fa")

        with FastaSequenceSource(file_path=path, mode="sequential") as source:
            assert source.get_sequence_by_start_end_strand("chr1", 1, 4) == "ACGT"

        assert not source.adaptor.metadata  # pylint: disable=protected-access


class TestFastaSequenceSourceSequence:
    """Tests for FASTA-backed sequence loading."""

    def test_sequential_mode_loads_requested_span(self, tmp_path: Path) -> None:
        """Test sequential mode returns the requested sequence span."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = FastaSequenceSource(file_path=path, mode="sequential")

        assert source.get_sequence_by_start_end_strand("chr1", 2, 5) == "CGTT"

    def test_sequential_mode_reverse_complements_requested_span(self, tmp_path: Path) -> None:
        """Test sequential mode returns sequence on the requested strand."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = FastaSequenceSource(file_path=path, mode="sequential")

        assert source.get_sequence_by_start_end_strand("chr1", 2, 5, strand=Strand.REVERSE) == "AACG"

    def test_sequential_mode_pads_bases_after_sequence_end(self, tmp_path: Path) -> None:
        """Test sequential mode pads missing trailing bases with Ns."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = FastaSequenceSource(file_path=path, mode="sequential")

        assert source.get_sequence_by_start_end_strand("chr1", 7, 10) == "CANN"

    def test_sequential_mode_returns_ns_for_missing_region(self, tmp_path: Path) -> None:
        """Test sequential mode returns Ns for an unknown sequence region."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = FastaSequenceSource(file_path=path, mode="sequential")

        assert source.get_sequence_by_start_end_strand("chr3", 1, 4) == "NNNN"

    def test_query_mode_loads_requested_span(self, tmp_path: Path) -> None:
        """Test query mode returns the requested sequence span."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        source = FastaSequenceSource(file_path=path, mode="query")

        assert source.get_sequence_by_start_end_strand("chr1", 2, 5) == "CGTT"

    def test_query_mode_reverse_complements_requested_span(self, tmp_path: Path) -> None:
        """Test query mode returns sequence on the requested strand."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        source = FastaSequenceSource(file_path=path, mode="query")

        assert source.get_sequence_by_start_end_strand("chr1", 2, 5, strand=Strand.REVERSE) == "AACG"

    def test_query_mode_pads_bases_after_sequence_end(self, tmp_path: Path) -> None:
        """Test query mode pads missing trailing bases with Ns."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        source = FastaSequenceSource(file_path=path, mode="query")

        assert source.get_sequence_by_start_end_strand("chr1", 7, 10) == "CANN"

    def test_query_mode_returns_ns_for_missing_region(self, tmp_path: Path) -> None:
        """Test query mode returns Ns for an unknown sequence region."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        source = FastaSequenceSource(file_path=path, mode="query")

        assert source.get_sequence_by_start_end_strand("chr3", 1, 4) == "NNNN"

    def test_sequential_sequences_are_loaded_once(self, tmp_path: Path) -> None:
        """Test sequential mode keeps loaded FASTA records in memory."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = FastaSequenceSource(file_path=path, mode="sequential")

        assert source.get_sequence_by_start_end_strand("chr1", 1, 4) == "ACGT"
        loaded_sequences = source._sequences  # pylint: disable=protected-access
        assert source.get_sequence_by_start_end_strand("chr2", 1, 4) == "AAAA"
        assert source._sequences is loaded_sequences  # pylint: disable=protected-access


class TestFastaSequenceSourceLength:
    """Tests for FASTA sequence region length lookup."""

    def test_length_is_resolved_from_fasta_metadata(self, tmp_path: Path) -> None:
        """Test sequence region length is read from parser metadata."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = FastaSequenceSource(file_path=path, mode="sequential")

        assert source.get_seq_region_length("chr1") == 8
        assert source.get_seq_region_length("chr2") == 4

    def test_length_is_cached_by_sequence_source(self, tmp_path: Path) -> None:
        """Test sequence region length lookup uses the shared source cache."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = FastaSequenceSource(file_path=path, mode="sequential")

        assert source.get_seq_region_length("chr1") == 8
        assert source._seq_region_lengths == {"chr1": 8}  # pylint: disable=protected-access

    def test_length_returns_none_for_missing_region(self, tmp_path: Path) -> None:
        """Test missing sequence region length returns None."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = FastaSequenceSource(file_path=path, mode="sequential")

        assert source.get_seq_region_length("chr3") is None

    def test_length_comes_from_metadata_without_loading_sequences(self, tmp_path: Path) -> None:
        """Test length lookup does not load FASTA records."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = FastaSequenceSource(file_path=path, mode="sequential")
        source._parser = MetadataOnlyParser(  # pylint: disable=protected-access
            {"sequence_lengths": {"chr1": 8}}
        )

        assert source.get_seq_region_length("chr1") == 8
        assert source._sequences is None  # pylint: disable=protected-access
