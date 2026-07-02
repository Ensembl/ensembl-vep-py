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
"""Unit tests for the FASTA parser."""

from pathlib import Path

from ensembl.vep.parsers.fasta import FastaParser


SAMPLE = ">chr1\nACGTTGCA\n>chr2\nAAAA\n"


class TestFastaParserInitialization:
    """Test FASTA parser initialization."""

    def test_concrete_parser_sets_format(self) -> None:
        """Test parser stores the parser format."""

        parser = FastaParser()
        assert parser.format == "fasta"


class TestFastaParserRecords:
    """Tests for FastaParser records and parsing helpers."""

    def test_parse_text_returns_records(self) -> None:
        """Test parsing FASTA from text."""

        records = FastaParser.parse_text(SAMPLE)
        assert len(records) == 2
        assert records[0].seqid == "chr1"
        assert records[0].sequence == "ACGTTGCA"

    def test_parse_file_returns_records(self, tmp_path: Path) -> None:
        """Test parsing FASTA from a file."""

        fasta_path = tmp_path / "sample.fa"
        fasta_path.write_text(SAMPLE, encoding="utf-8")

        records = FastaParser.parse_file(fasta_path)
        assert len(records) == 2
        assert records[1].seqid == "chr2"
        assert records[1].sequence == "AAAA"

    def test_records_respects_limit(self, tmp_path: Path) -> None:
        """Test records stops once the limit is reached."""

        fasta_path = tmp_path / "sample.fa"
        fasta_path.write_text(SAMPLE, encoding="utf-8")

        parser = FastaParser(file_path=fasta_path, mode="sequential")
        assert len(parser.records(limit=1)) == 1


class TestFastaParserSeek:
    """Tests for indexed seek."""

    def test_query_seek_returns_requested_slice(self, tmp_path: Path) -> None:
        """Test indexed seek returns the requested subsequence."""

        fasta_path = tmp_path / "sample.fa"
        fasta_path.write_text(SAMPLE, encoding="utf-8")
        (tmp_path / "sample.fa.fai").write_text("chr1\t8\t6\t8\t9\nchr2\t4\t23\t4\t5\n", encoding="utf-8")

        parser = FastaParser(file_path=fasta_path, mode="query")
        assert parser.seek("chr1:2-5")
        record = parser.next()
        assert record is not None
        assert record.seqid == "chr1"
        assert record.sequence == "CGTT"

    def test_query_seek_whole_sequence_resolves_length(self, tmp_path: Path) -> None:
        """Test a whole-sequence query returns the full record."""

        fasta_path = tmp_path / "sample.fa"
        fasta_path.write_text(SAMPLE, encoding="utf-8")
        (tmp_path / "sample.fa.fai").write_text("chr1\t8\t6\t8\t9\nchr2\t4\t23\t4\t5\n", encoding="utf-8")

        parser = FastaParser(file_path=fasta_path, mode="query")
        assert parser.seek("chr1")
        record = parser.next()
        assert record is not None
        assert record.sequence == "ACGTTGCA"
        assert parser.metadata["sequence_lengths"]["chr1"] == 8

    def test_query_seek_missing_sequence_returns_false(self, tmp_path: Path) -> None:
        """Test seek returns False when the sequence is missing."""

        fasta_path = tmp_path / "sample.fa"
        fasta_path.write_text(SAMPLE, encoding="utf-8")
        (tmp_path / "sample.fa.fai").write_text("chr1\t8\t6\t8\t9\nchr2\t4\t23\t4\t5\n", encoding="utf-8")

        parser = FastaParser(file_path=fasta_path, mode="query")
        assert not parser.seek("chr3")
