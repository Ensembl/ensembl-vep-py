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

import os
from pathlib import Path

from attrs.exceptions import FrozenAttributeError
import pytest

from ensembl.vep.common import Strand
from ensembl.vep.parsers.fasta import FastaParser
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource

SAMPLE = """>chr1
ACGTTGCA
>chr2
AAAA
"""
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


def _make_fasta_source(path: Path, mode: str) -> FastaSequenceSource:
    """Create a sequence source with a parser opened for the test FASTA."""
    return FastaSequenceSource(
        file_path=path,
        adaptor=FastaParser(file_path=path, mode=mode),
    )


class MetadataOnlyParser(FastaParser):
    """Minimal parser stub that only exposes metadata."""

    def __init__(self, metadata: dict[str, object]) -> None:  # pylint: disable=super-init-not-called
        """Initialize the parser stub."""

        self._test_metadata = metadata
        self.closed = False

    @property
    def metadata(self) -> dict[str, object]:
        """Return parser metadata."""

        return self._test_metadata

    @property
    def reader_mode(self) -> str:
        """Return the parser mode used by this stub."""
        return "sequential"

    def open(self, file_path: str | Path, mode: str = "auto") -> FastaParser:
        """Open the parser stub."""

        return self

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
            _ = FastaSequenceSource("sample.fa")  # type: ignore[arg-type,misc]

    def test_initialization_sets_core_attributes(self) -> None:
        """Test FastaSequenceSource stores source and parser attributes."""

        source = FastaSequenceSource()

        assert source.file_path is None
        assert source.adaptor.format == "fasta"  # pylint: disable=protected-access
        assert source.adaptor.reader_mode == "auto"  # pylint: disable=protected-access
        assert source._metadata_lengths is None  # pylint: disable=protected-access


class TestFastaSequenceSourceOpenClose:
    """Tests for FastaSequenceSource open and close."""

    def test_close_clears_metadata_lengths(self, tmp_path: Path) -> None:
        """Test close clears FASTA-specific and derived source caches."""

        path = _write_indexed_fasta(tmp_path / "sample.fa", SAMPLE)
        source = _make_fasta_source(path, "query")

        assert source.get_sequence_by_start_end_strand("chr1", 1, 4) == "ACGT"
        assert source.get_seq_region_length("chr1") == 8
        source.close()

        assert source._metadata_lengths is None  # pylint: disable=protected-access
        assert not source._seq_region_lengths  # pylint: disable=protected-access

    def test_context_manager_closes_source(self, tmp_path: Path) -> None:
        """Test context manager returns and closes the source."""

        path = _write_fasta(tmp_path / "sample.fa")

        source = _make_fasta_source(path, "sequential")
        with source:
            assert source.get_sequence_by_start_end_strand("chr1", 1, 4) == "ACGT"

        assert not source.adaptor.metadata  # pylint: disable=protected-access

    def test_changing_file_path_is_error(self, tmp_path: Path) -> None:
        """Test that the file_path attribute is frozen."""

        first_path = _write_indexed_fasta(tmp_path / "sample.fa", SAMPLE)
        source = _make_fasta_source(first_path, "query")
        assert source.get_sequence_by_start_end_strand("chr1", 1, 4) == "ACGT"

        second_path = _write_indexed_fasta(
            tmp_path / "sample-2.fa", ">chr1\nAAAAAAAA\n", "chr1\t8\t6\t8\t9\n"
        )
        with pytest.raises(FrozenAttributeError, match="can't set attribute"):
            source.file_path = second_path


class TestFastaSequenceSourceSequence:
    """Tests for FASTA-backed sequence loading."""

    def test_query_mode_reuses_open_reader_across_spans(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test consecutive span fetches open the indexed FASTA once."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        source = _make_fasta_source(path, "query")
        opened: list[str] = []
        original_open = FastaParser.open

        def counting_open(parser: FastaParser, file_path: str, mode: str = "auto") -> FastaParser:
            opened.append(mode)
            return original_open(parser, file_path, mode)

        monkeypatch.setattr(FastaParser, "open", counting_open)

        assert source.get_sequence_by_start_end_strand("chr1", 1, 2) == "AC"
        assert source.get_sequence_by_start_end_strand("chr1", 3, 4) == "GT"
        assert opened == ["auto"]

    def test_sequential_mode_reuses_open_reader_across_spans(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test consecutive span fetches open the unindexed FASTA once."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = _make_fasta_source(path, "sequential")
        opened: list[str] = []
        original_open = FastaParser.open

        def counting_open(parser: FastaParser, file_path: str, mode: str = "auto") -> FastaParser:
            opened.append(mode)
            return original_open(parser, file_path, mode)

        monkeypatch.setattr(FastaParser, "open", counting_open)

        assert source.get_sequence_by_start_end_strand("chr1", 1, 2) == "AC"
        assert source.get_sequence_by_start_end_strand("chr2", 3, 4) == "AA"
        assert opened == ["auto"]
        assert source.adaptor.reader_mode == "sequential"

    def test_sequential_mode_warns_once(self, tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
        """Test an unindexed FASTA logs one warning across spans and length lookups."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = _make_fasta_source(path, "sequential")

        with caplog.at_level("WARNING"):
            assert source.get_sequence_by_start_end_strand("chr1", 1, 2) == "AC"
            assert source.get_seq_region_length("chr2") == 4
            assert source.get_sequence_by_start_end_strand("chr1", 3, 4) == "GT"

        assert [record.message for record in caplog.records if "no usable FASTA index" in record.message] == [
            f"{path} has no usable FASTA index; sequential fallback is bounded but may be slow"
        ]

    def test_query_mode_reopens_reader_in_another_process(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test a reader inherited from another process id is opened again."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        source = _make_fasta_source(path, "query")
        assert source.get_sequence_by_start_end_strand("chr1", 1, 2) == "AC"

        monkeypatch.setattr(os, "getpid", lambda: -1)

        assert source.get_sequence_by_start_end_strand("chr1", 3, 4) == "GT"
        assert source._adaptor_pid == -1  # pylint: disable=protected-access

    def test_sequential_mode_loads_requested_span(self, tmp_path: Path) -> None:
        """Test sequential mode returns the requested sequence span."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = _make_fasta_source(path, "sequential")

        assert source.get_sequence_by_start_end_strand("chr1", 2, 5) == "CGTT"

    def test_sequential_mode_reverse_complements_requested_span(self, tmp_path: Path) -> None:
        """Test sequential mode returns sequence on the requested strand."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = _make_fasta_source(path, "sequential")

        assert source.get_sequence_by_start_end_strand("chr1", 2, 5, strand=Strand.REVERSE) == "AACG"

    def test_sequential_mode_pads_bases_after_sequence_end(self, tmp_path: Path) -> None:
        """Test sequential mode pads missing trailing bases with Ns."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = _make_fasta_source(path, "sequential")

        assert source.get_sequence_by_start_end_strand("chr1", 7, 10) == "CANN"

    def test_sequential_mode_returns_ns_for_missing_region(self, tmp_path: Path) -> None:
        """Test sequential mode returns Ns for an unknown sequence region."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = _make_fasta_source(path, "sequential")

        assert source.get_sequence_by_start_end_strand("chr3", 1, 4) == "NNNN"

    def test_query_mode_loads_requested_span(self, tmp_path: Path) -> None:
        """Test query mode returns the requested sequence span."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        source = _make_fasta_source(path, "query")

        assert source.get_sequence_by_start_end_strand("chr1", 2, 5) == "CGTT"

    def test_query_mode_reverse_complements_requested_span(self, tmp_path: Path) -> None:
        """Test query mode returns sequence on the requested strand."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        source = _make_fasta_source(path, "query")

        assert source.get_sequence_by_start_end_strand("chr1", 2, 5, strand=Strand.REVERSE) == "AACG"

    def test_query_mode_pads_bases_after_sequence_end(self, tmp_path: Path) -> None:
        """Test query mode pads missing trailing bases with Ns."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        source = _make_fasta_source(path, "query")

        assert source.get_sequence_by_start_end_strand("chr1", 7, 10) == "CANN"

    def test_query_mode_returns_ns_for_missing_region(self, tmp_path: Path) -> None:
        """Test query mode returns Ns for an unknown sequence region."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        source = _make_fasta_source(path, "query")

        assert source.get_sequence_by_start_end_strand("chr3", 1, 4) == "NNNN"

    def test_sequential_spans_do_not_retain_complete_reference_sequences(self, tmp_path: Path) -> None:
        """Test sequential lookup retains only its requested result span."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = _make_fasta_source(path, "sequential")

        assert source.get_sequence_by_start_end_strand("chr1", 1, 4) == "ACGT"
        assert source.get_sequence_by_start_end_strand("chr2", 1, 4) == "AAAA"
        assert not hasattr(source, "_sequences")

    def test_no_file_returns_ns(self) -> None:
        """Test if source.file_path is None, _fetch_sequence_chunk returns Ns"""
        source = FastaSequenceSource()
        # pylint: disable-next=protected-access
        assert source._fetch_sequence_chunk(seq_region_name="chr1", start=101, end=110) == "NNNNNNNNNN"

    def test_no_file_metadata_is_empty_frozenset(self) -> None:
        """Test if source.file_path is None, _fetch_sequence_chunk returns Ns"""
        source = FastaSequenceSource()
        # pylint: disable-next=protected-access
        assert source.sequence_region_names == frozenset()


class TestFastaSequenceSourceLength:
    """Tests for FASTA sequence region length lookup."""

    def test_length_is_resolved_from_fasta_metadata(self, tmp_path: Path) -> None:
        """Test sequence region length is read from parser metadata."""

        path = _write_indexed_fasta(tmp_path / "sample.fa", SAMPLE)
        source = _make_fasta_source(path, "query")

        assert source.get_seq_region_length("chr1") == 8
        assert source.get_seq_region_length("chr2") == 4

    def test_length_is_cached_by_sequence_source(self, tmp_path: Path) -> None:
        """Test sequence region length lookup uses the shared source cache."""

        path = _write_indexed_fasta(tmp_path / "sample.fa", SAMPLE)
        source = _make_fasta_source(path, "query")

        assert source.get_seq_region_length("chr1") == 8
        assert source._seq_region_lengths == {"chr1": 8}  # pylint: disable=protected-access

    def test_length_returns_none_for_missing_region(self, tmp_path: Path) -> None:
        """Test missing sequence region length returns None."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        source = _make_fasta_source(path, "query")

        assert source.get_seq_region_length("chr3") is None

    def test_length_is_scanned_from_an_unindexed_fasta(self, tmp_path: Path) -> None:
        """Test sequential fallback derives authoritative lengths without retaining records."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = _make_fasta_source(path, "sequential")

        assert source.get_seq_region_length("chr1") == 8
        assert source.sequence_region_names == frozenset({"chr1", "chr2"})

    def test_length_comes_from_metadata_without_scanning_sequences(self, tmp_path: Path) -> None:
        """Test length lookup does not load FASTA records."""

        path = _write_fasta(tmp_path / "sample.fa")
        source = FastaSequenceSource(
            file_path=path,
            adaptor=MetadataOnlyParser({"sequence_lengths": {"chr1": 8}}),
        )

        assert source.get_seq_region_length("chr1") == 8
        assert source._metadata_lengths == {"chr1": 8}  # pylint: disable=protected-access
