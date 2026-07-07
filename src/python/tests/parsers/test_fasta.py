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
"""Unit tests for the FASTA parser interface."""

from pathlib import Path
import pytest

from ensembl.vep.parsers.fasta import FastaParser, FastaRecord


SAMPLE = (
    ">chr1\n"
    "ACGTTGCA\n"
    ">chr2\n"
    "AAAA\n"
)
DESCRIBED_SAMPLE = (
    ">chr1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF\n"
    "ACGTTGCA\n"
    ">chr2 alternate sequence\n"
    "AAAA\n"
)
SAMPLE_FAI = "chr1\t8\t6\t8\t9\nchr2\t4\t21\t4\t5\n"
DESCRIBED_SAMPLE_FAI = "chr1\t8\t66\t8\t9\nchr2\t4\t106\t4\t5\n"


def _write_fasta(path: Path, data: str = SAMPLE) -> Path:
    """Write a FASTA fixture file."""

    path.write_text(data, encoding="utf-8")
    return path


def _write_indexed_fasta(path: Path, data: str = SAMPLE, index: str = SAMPLE_FAI) -> Path:
    """Write a FASTA fixture with an accompanying FAI index."""

    path.write_text(data, encoding="utf-8")
    path.with_suffix(path.suffix + ".fai").write_text(index, encoding="utf-8")
    return path


def _assert_chr1_record(record: FastaRecord | None) -> None:
    """Assert a parsed record matches the chr1 fixture."""

    assert record is not None
    assert record.seqid == "chr1"
    assert record.start == 1
    assert record.end == 8
    assert record.sequence == "ACGTTGCA"


def _assert_chr2_record(record: FastaRecord | None) -> None:
    """Assert a parsed record matches the chr2 fixture."""

    assert record is not None
    assert record.seqid == "chr2"
    assert record.start == 1
    assert record.end == 4
    assert record.sequence == "AAAA"


class TestFastaParserInitialization:
    """Test FastaParser initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test FastaParser initialization without keywords raises error."""

        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = FastaParser("sample.fa")  # type: ignore[misc]

    def test_initialization_sets_core_attributes(self) -> None:
        """Test FastaParser sets parser attributes."""

        parser = FastaParser()

        assert parser.format == "fasta"
        assert parser.file_path is None
        assert parser.mode == "auto"
        assert parser.reader_mode == "auto"
        assert parser.metadata == {}

    def test_initialization_with_file_path_opens_parser(self, tmp_path: Path) -> None:
        """Test file_path initialization opens the parser."""

        path = _write_fasta(tmp_path / "sample.fa")
        parser = FastaParser(file_path=path, mode="sequential")

        _assert_chr1_record(parser.next())


class TestFastaParserMetadata:
    """Tests for FastaParser.metadata."""

    def test_metadata_returns_headers_lengths_and_descriptions(self, tmp_path: Path) -> None:
        """Test metadata exposes FASTA headers, lengths and descriptions."""

        path = _write_fasta(tmp_path / "sample.fa", DESCRIBED_SAMPLE)
        parser = FastaParser(file_path=path, mode="sequential")

        assert parser.metadata == {
            "header": [
                ">chr1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF",
                ">chr2 alternate sequence",
            ],
            "sequence_lengths": {"chr1": 8, "chr2": 4},
            "sequence_descriptions": {
                "chr1": "dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF",
                "chr2": "alternate sequence",
            },
        }

    def test_metadata_is_cached(self, tmp_path: Path) -> None:
        """Test metadata is stable after first access."""

        path = _write_fasta(tmp_path / "sample.fa", DESCRIBED_SAMPLE)
        parser = FastaParser(file_path=path, mode="sequential")

        metadata = parser.metadata

        assert parser.metadata == metadata

    def test_open_refreshes_metadata_cache(self, tmp_path: Path) -> None:
        """Test open exposes metadata for the newly opened file."""

        first_path = _write_fasta(tmp_path / "first.fa", ">chr1 first description\nACGT\n")
        second_path = _write_fasta(tmp_path / "second.fa", ">chr2 second description\nAAAA\n")
        parser = FastaParser(file_path=first_path, mode="sequential")

        assert parser.metadata["sequence_descriptions"] == {"chr1": "first description"}
        assert parser.open(second_path, mode="sequential") is parser
        assert parser.metadata["sequence_descriptions"] == {"chr2": "second description"}

    def test_metadata_returns_empty_dictionary_after_close(self, tmp_path: Path) -> None:
        """Test metadata is empty when the parser is closed."""

        path = _write_fasta(tmp_path / "sample.fa", DESCRIBED_SAMPLE)
        parser = FastaParser(file_path=path, mode="sequential")

        assert parser.metadata["sequence_lengths"] == {"chr1": 8, "chr2": 4}

        parser.close()

        assert not parser.metadata

    def test_metadata_does_not_consume_sequential_stream(self, tmp_path: Path) -> None:
        """Test metadata reads do not advance the active sequential stream."""

        path = _write_fasta(tmp_path / "sample.fa", DESCRIBED_SAMPLE)
        parser = FastaParser(file_path=path, mode="sequential")

        _assert_chr1_record(parser.next())
        assert parser.metadata["sequence_lengths"] == {"chr1": 8, "chr2": 4}
        _assert_chr2_record(parser.next())

    def test_metadata_in_query_mode_reads_source_headers(self, tmp_path: Path) -> None:
        """Test query-mode metadata is read from the FASTA source file."""

        path = _write_indexed_fasta(tmp_path / "sample.fa", DESCRIBED_SAMPLE, DESCRIBED_SAMPLE_FAI)
        parser = FastaParser(file_path=path, mode="query")

        assert parser.reader_mode == "query"
        assert parser.metadata["header"] == [
            ">chr1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF",
            ">chr2 alternate sequence",
        ]
        assert parser.metadata["sequence_lengths"] == {"chr1": 8, "chr2": 4}


class TestFastaParserOpenClose:
    """Tests for FastaParser open and close."""

    def test_open_returns_parser(self, tmp_path: Path) -> None:
        """Test open returns the parser instance."""

        path = _write_fasta(tmp_path / "sample.fa")
        parser = FastaParser(mode="sequential")

        assert parser.open(path) is parser
        _assert_chr1_record(parser.next())

    def test_open_overrides_mode(self, tmp_path: Path) -> None:
        """Test open overrides the parser mode."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        parser = FastaParser(file_path=None, mode="auto")
        assert parser.reader_mode == "auto"

        assert parser.open(path, mode="sequential") is parser
        assert parser.reader_mode == "sequential"

    def test_open_auto_query_then_auto_seq_changes_mode(self, tmp_path: Path) -> None:
        """Test auto mode is recomputed each time a file is opened."""

        indexed_path = _write_indexed_fasta(tmp_path / "indexed.fa")
        plain_path = _write_fasta(tmp_path / "plain.fa")
        parser = FastaParser(file_path=None, mode="auto")

        assert parser.open(indexed_path) is parser
        assert parser.reader_mode == "query"

        assert parser.open(plain_path) is parser
        assert parser.reader_mode == "sequential"

    def test_close_closes_backend_reader(self, tmp_path: Path) -> None:
        """Test close prevents subsequent backend reads."""

        path = _write_fasta(tmp_path / "sample.fa")
        parser = FastaParser(file_path=path, mode="sequential")

        parser.close()

        with pytest.raises(RuntimeError, match="not open"):
            parser.next()

    def test_context_manager_closes_parser(self, tmp_path: Path) -> None:
        """Test context manager returns and closes the parser."""

        path = _write_fasta(tmp_path / "sample.fa")

        with FastaParser(file_path=path, mode="sequential") as parser:
            _assert_chr1_record(parser.next())

        with pytest.raises(RuntimeError, match="not open"):
            parser.next()


class TestFastaParserSeek:
    """Tests for FastaParser.seek."""

    def test_fasta_parser_seek_primes_next(self, tmp_path: Path) -> None:
        """Test seek returns availability and next returns the requested sequence."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        parser = FastaParser(file_path=path, mode="query")

        assert parser.seek(seq_region_name="chr1", start=2, end=5) is True
        record = parser.next()
        assert record is not None
        assert record.seqid == "chr1"
        assert record.start == 2
        assert record.end == 5
        assert record.sequence == "CGTT"
        assert parser.next() is None

    def test_fasta_parser_seek_accepts_region_string(self, tmp_path: Path) -> None:
        """Test seek accepts an explicit region string."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        parser = FastaParser(file_path=path, mode="query")

        assert parser.seek("chr2:1-3") is True
        record = parser.next()
        assert record is not None
        assert record.seqid == "chr2"
        assert record.sequence == "AAA"

    def test_fasta_parser_seek_accepts_whole_sequence(self, tmp_path: Path) -> None:
        """Test seek accepts a seq_region_name-only region."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        parser = FastaParser(file_path=path, mode="query")

        assert parser.seek(seq_region_name="chr1") is True
        _assert_chr1_record(parser.next())

    def test_fasta_parser_seek_empty_region(self, tmp_path: Path) -> None:
        """Test seek returns false when the sequence is missing."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        parser = FastaParser(file_path=path, mode="query")

        assert parser.seek("chr3") is False
        assert parser.next() is None

    def test_fasta_parser_seek_pads_beyond_sequence_end(self, tmp_path: Path) -> None:
        """Test indexed seek fills missing trailing bases with Ns."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        parser = FastaParser(file_path=path, mode="query")

        assert parser.seek("chr1:7-10") is True
        record = parser.next()
        assert record is not None
        assert record.sequence == "CANN"

    @pytest.mark.parametrize(
        "kwargs",
        [
            {},
            {"region": "chr1:1-4", "seq_region_name": "chr1"},
            {"seq_region_name": "chr1", "end": 4},
            {"seq_region_name": "chr1", "start": 0, "end": 4},
            {"seq_region_name": "chr1", "start": 1, "end": 0},
        ],
    )
    def test_fasta_parser_seek_rejects_invalid_region_arguments(self, kwargs: dict[str, object]) -> None:
        """Test seek rejects invalid region inputs before querying the backend."""

        parser = FastaParser()

        with pytest.raises(ValueError):
            parser.seek(**kwargs)  # type: ignore[arg-type]

    def test_fasta_parser_seek_rejects_invalid_region_string(self, tmp_path: Path) -> None:
        """Test invalid region strings are rejected by the query backend."""

        path = _write_indexed_fasta(tmp_path / "sample.fa")
        parser = FastaParser(file_path=path, mode="query")

        with pytest.raises(ValueError):
            parser.seek("chr1:5-2")

    def test_sequential_mode_rejects_seek(self, tmp_path: Path) -> None:
        """Test sequential reader mode keeps seek disabled."""

        path = _write_fasta(tmp_path / "sample.fa")
        parser = FastaParser(file_path=path, mode="sequential")

        with pytest.raises(RuntimeError, match="sequential"):
            parser.seek(seq_region_name="chr1", start=1, end=4)


class TestFastaParserRecords:
    """Tests for FastaParser.records."""

    def test_records_returns_remaining_records(self, tmp_path: Path) -> None:
        """Test records returns all remaining records from the stream."""

        path = _write_fasta(tmp_path / "sample.fa")
        parser = FastaParser(file_path=path, mode="sequential")

        _assert_chr1_record(parser.next())

        records = parser.records()
        assert len(records) == 1
        _assert_chr2_record(records[0])

    def test_records_respects_limit(self, tmp_path: Path) -> None:
        """Test records stops once the limit is reached."""

        path = _write_fasta(tmp_path / "sample.fa")
        parser = FastaParser(file_path=path, mode="sequential")

        records = parser.records(limit=1)
        assert len(records) == 1
        _assert_chr1_record(records[0])
        _assert_chr2_record(parser.next())

    def test_records_zero_limit_does_not_consume_record(self, tmp_path: Path) -> None:
        """Test zero limit leaves the next record available."""

        path = _write_fasta(tmp_path / "sample.fa")
        parser = FastaParser(file_path=path, mode="sequential")

        assert parser.records(limit=0) == []
        _assert_chr1_record(parser.next())


class TestFastaParserIteration:
    """Tests for FastaParser iteration."""

    def test_parser_is_its_own_iterator(self) -> None:
        """Test iter returns the parser instance."""

        parser = FastaParser()

        assert iter(parser) is parser

    def test_iteration_streams_records(self, tmp_path: Path) -> None:
        """Test iteration streams records from the parser."""

        path = _write_fasta(tmp_path / "sample.fa")
        parser = FastaParser(file_path=path, mode="sequential")

        assert [record.seqid for record in parser] == ["chr1", "chr2"]

    def test_iteration_raises_stop_iteration(self, tmp_path: Path) -> None:
        """Test iteration raises StopIteration when no more records are available."""

        path = _write_fasta(tmp_path / "sample.fa")
        parser = FastaParser(file_path=path, mode="sequential")

        with pytest.raises(StopIteration):
            while True:
                next(parser)


class TestFastaParserParsingTextAndFiles:
    """Tests for FastaParser parse helpers."""

    def test_parse_text_returns_interface_records(self) -> None:
        """Test parse_text exposes FASTA rows as parser-interface records."""

        records = FastaParser.parse_text(SAMPLE)

        assert len(records) == 2
        _assert_chr1_record(records[0])
        _assert_chr2_record(records[1])

    def test_parse_text_respects_limit(self) -> None:
        """Test parse_text stops once the limit is reached."""

        records = FastaParser.parse_text(SAMPLE, limit=1)

        assert len(records) == 1
        _assert_chr1_record(records[0])

    def test_parse_text_zero_limit_returns_empty_list(self) -> None:
        """Test parse_text respects a zero limit."""

        assert FastaParser.parse_text(SAMPLE, limit=0) == []

    def test_parse_file_returns_interface_records(self, tmp_path: Path) -> None:
        """Test parse_file exposes FASTA rows as parser-interface records."""

        path = _write_fasta(tmp_path / "sample.fa")
        records = FastaParser.parse_file(path)

        assert len(records) == 2
        _assert_chr1_record(records[0])
        _assert_chr2_record(records[1])

    def test_parse_file_respects_limit(self, tmp_path: Path) -> None:
        """Test parse_file stops once the limit is reached."""

        path = _write_fasta(tmp_path / "sample.fa")
        records = FastaParser.parse_file(path, limit=1)

        assert len(records) == 1
        _assert_chr1_record(records[0])

    def test_parse_file_missing_path_errors(self) -> None:
        """Test parse_file reports missing files."""

        with pytest.raises(OSError, match="No such file|os error 2"):
            FastaParser.parse_file("/definitely/missing/sample.fa")
