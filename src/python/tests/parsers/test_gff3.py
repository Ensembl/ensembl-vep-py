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
"""Unit tests for the GFF3 parser interface."""

from pathlib import Path
from attr import exceptions as attr_ex

import pytest

from ensembl.vep.parsers.gff3 import GFF3Parser, GFF3Record

HEADER_LINES = ["##gff-version 3", "# comment 1"]
CHR1_GENE_ROW = "chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1;Name=BRCA2;Alias=a,b"
CHR1_CDS_ROW = "chr1\tensembl\tCDS\t12\t18\t4.5\t-\t2\tID=cds:1;Parent=gene:1"
CHR1_EXON_ROW = "chr1\ttabix\texon\t15\t18\t.\t+\t0\tID=exon:1"
CHR2_GENE_ROW = "chr2\ttabix\tgene\t100\t120\t.\t-\t.\tID=gene:2"
SAMPLE = "\n".join([HEADER_LINES[0], CHR1_GENE_ROW, CHR1_CDS_ROW, "##FASTA", ">chr1", "ACGT"]) + "\n"
SAMPLE_WITH_COMMENT = "\n".join([*HEADER_LINES, CHR1_GENE_ROW, CHR1_CDS_ROW]) + "\n"

# the default SAMPLE string creates a gff file with the chr1 gene row and the chr1 CDS row
# whereas the mock tabix executable returns the chr1 gene row and the chr1 exon row for queries on chr1
def _write_gff3(path: Path, data: str = SAMPLE) -> Path:
    """Write a GFF3 fixture file."""

    path.write_text(data)
    return path


def _write_query_gff3(path: Path, data: str = SAMPLE_WITH_COMMENT) -> Path:
    """Write a query-mode GFF3 fixture with a placeholder tabix index."""

    path.write_text(data)
    path.with_suffix(path.suffix + ".tbi").write_text("")
    return path


def _write_mock_tabix(tmp_path: Path, body: str) -> Path:
    """Write a mock tabix executable."""

    tabix = tmp_path / "tabix"
    tabix.write_text(f"#!/bin/sh\n{body}\n")
    tabix.chmod(0o755)
    return tabix


def _write_query_mock_tabix(tmp_path: Path) -> Path:
    """Write a mock tabix executable that returns deterministic GFF3 records."""

    return _write_mock_tabix(
        tmp_path,
        f"""
if [ "$1" = "--version" ]; then
    exit 0
fi

case "$2" in
    chr1)
        printf '%s\\n' '{CHR1_GENE_ROW}'
        printf '%s\\n' '{CHR1_EXON_ROW}'
        ;;
    chr1:10-20)
        printf '%s\\n' '{CHR1_GENE_ROW}'
        printf '%s\\n' '{CHR1_EXON_ROW}'
        ;;
    chr2:100-120)
        printf '%s\\n' '{CHR2_GENE_ROW}'
        ;;
esac
""",
    )


def _assert_gene_record(record: GFF3Record) -> None:
    """Assert a parsed record matches the chr1 gene fixture."""

    assert record is not None
    assert record.seqid == "chr1"
    assert record.source == "ensembl"
    assert record.feature_type == "gene"
    assert record.start == 10
    assert record.end == 20
    assert record.score is None
    assert record.strand_symbol == "+"
    assert record.phase_value is None
    assert record.raw_attributes == "ID=gene:1;Name=BRCA2;Alias=a,b"
    assert record.attributes == {"ID": ("gene:1",), "Name": ("BRCA2",), "Alias": ("a", "b")}


class TestGFF3ParserInitialization:
    """Test GFF3Parser initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test GFF3Parser initialization without keywords raises error."""

        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = GFF3Parser("annotations.gff3")  # type: ignore[misc]

    def test_initialization_sets_core_attributes(self) -> None:
        """Test GFF3Parser sets parser attributes."""

        parser = GFF3Parser()

        assert parser.format == "gff3"
        assert parser.file_path is None
        assert parser.tabix_path == "tabix"
        assert parser.mode == "auto"
        assert parser.reader_mode == "auto"

    def test_initialization_with_file_path_opens_parser(self, tmp_path: Path) -> None:
        """Test file_path initialization opens the parser."""

        path = _write_gff3(tmp_path / "sample.gff3")
        parser = GFF3Parser(file_path=path)

        _assert_gene_record(parser.next())

    def test_mode_is_frozen(self) -> None:
        """Test attribute mode is settable on object creation and them immutable"""
        parser = GFF3Parser(mode="sequential")

        assert parser.mode == "sequential", "mode attribute set incorrectly"
        assert parser.reader_mode == "sequential", "reader_mode attribute set incorrectly"
        # check that mode is not settable after object creation
        with pytest.raises(attr_ex.FrozenAttributeError, match="can't set attribute"):
            parser.mode = "query"


class TestGFF3ParserMetadata:
    """Tests for GFF3Parser.metadata."""

    def test_metadata_returns_header_lines_and_directives(self, tmp_path: Path) -> None:
        """Test metadata returns raw header lines and directive values."""

        path = _write_gff3(tmp_path / "sample.gff3", SAMPLE_WITH_COMMENT)
        parser = GFF3Parser(file_path=path, mode="sequential")

        assert parser.metadata == {"header": HEADER_LINES, "gff-version": "3"}

    def test_metadata_is_cached(self, tmp_path: Path) -> None:
        """Test metadata is cached after first access."""

        path = _write_gff3(tmp_path / "sample.gff3", SAMPLE_WITH_COMMENT)
        parser = GFF3Parser(file_path=path, mode="sequential")

        metadata = parser.metadata

        assert parser.metadata == metadata

    def test_open_refreshes_metadata_cache(self, tmp_path: Path) -> None:
        """Test open clears cached metadata for the previous file."""

        first_path = _write_gff3(
            tmp_path / "first.gff3",
            "##gff-version 3\nchr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1\n",
        )
        second_path = _write_gff3(
            tmp_path / "second.gff3",
            "##gff-version 4\nchr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1\n",
        )
        parser = GFF3Parser(file_path=first_path, mode="sequential")

        assert parser.metadata == {"header": ["##gff-version 3"], "gff-version": "3"}
        assert parser.open(second_path) is parser
        assert parser.metadata == {"header": ["##gff-version 4"], "gff-version": "4"}

    def test_metadata_returns_empty_dictionary_after_close(self, tmp_path: Path) -> None:
        """Test metadata is empty when the parser is closed."""

        path = _write_gff3(tmp_path / "sample.gff3", SAMPLE_WITH_COMMENT)
        parser = GFF3Parser(file_path=path, mode="sequential")

        assert parser.metadata == {"header": HEADER_LINES, "gff-version": "3"}

        parser.close()

        assert not parser.metadata

    def test_metadata_does_not_consume_sequential_stream(self, tmp_path: Path) -> None:
        """Test metadata reads do not advance the active sequential stream."""

        path = _write_gff3(tmp_path / "sample.gff3", SAMPLE_WITH_COMMENT)
        parser = GFF3Parser(file_path=path, mode="sequential")

        _assert_gene_record(parser.next())
        assert parser.metadata == {"header": HEADER_LINES, "gff-version": "3"}
        record = parser.next()
        assert record is not None
        assert record.feature_type == "CDS"

    def test_metadata_in_query_mode_reads_source_headers(self, tmp_path: Path) -> None:
        """Test query-mode metadata is read from the GFF3 source file."""

        path = _write_query_gff3(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = GFF3Parser(file_path=path, tabix_path=tabix)

        assert parser.reader_mode == "query"
        assert parser.metadata == {"header": HEADER_LINES, "gff-version": "3"}


class TestGFF3ParserOpenClose:
    """Tests for GFF3Parser open and close."""

    def test_open_returns_parser(self, tmp_path: Path) -> None:
        """Test open returns the parser instance."""

        path = _write_gff3(tmp_path / "sample.gff3")
        parser = GFF3Parser(mode="sequential")

        assert parser.open(path) is parser
        _assert_gene_record(parser.next())

    def test_close_closes_backend_reader(self, tmp_path: Path) -> None:
        """Test close prevents subsequent backend reads."""

        path = _write_gff3(tmp_path / "sample.gff3")
        parser = GFF3Parser(file_path=path, mode="sequential")

        parser.close()

        with pytest.raises(RuntimeError, match="not open"):
            parser.next()

    def test_context_manager_closes_parser(self, tmp_path: Path) -> None:
        """Test context manager returns and closes the parser."""

        path = _write_gff3(tmp_path / "sample.gff3")

        with GFF3Parser(file_path=path, mode="sequential") as parser:
            _assert_gene_record(parser.next())

        with pytest.raises(RuntimeError, match="not open"):
            parser.next()


class TestGFF3ParserSeek:
    """Tests for GFF3Parser.seek."""

    def test_gff3_parser_seek_primes_next(self, tmp_path: Path) -> None:
        """Test seek returns availability and next streams parsed records."""

        path = _write_query_gff3(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = GFF3Parser(file_path=path, tabix_path=tabix)

        assert parser.seek(seqid="chr1", start=10, end=20) is True
        _assert_gene_record(parser.next())
        assert [record.feature_type for record in parser.records()] == ["exon"]
        assert parser.next() is None

    def test_gff3_parser_seek_accepts_region_string(self, tmp_path: Path) -> None:
        """Test seek accepts an explicit region string."""

        path = _write_query_gff3(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = GFF3Parser(file_path=path, tabix_path=tabix)

        assert parser.seek("chr2:100-120") is True
        assert [record.seqid for record in parser.records()] == ["chr2"]

    def test_gff3_parser_seek_accepts_whole_sequence(self, tmp_path: Path) -> None:
        """Test seek accepts a seqid-only region."""

        path = _write_query_gff3(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = GFF3Parser(file_path=path, tabix_path=tabix)

        assert parser.seek(seqid="chr1") is True
        assert [record.feature_type for record in parser.records()] == ["gene", "exon"]

    def test_gff3_parser_seek_empty_region(self, tmp_path: Path) -> None:
        """Test seek returns false when no records are available."""

        path = _write_query_gff3(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = GFF3Parser(file_path=path, tabix_path=tabix)

        assert parser.seek("chr1:30-40") is False
        assert parser.next() is None

    @pytest.mark.parametrize(
        "kwargs",
        [
            {},
            {"region": "chr1:10-20", "seqid": "chr1"},
            {"seqid": "chr1", "start": 0, "end": 20},
            {"seqid": "chr1", "start": 10, "end": 0},
        ],
    )
    def test_gff3_parser_seek_rejects_invalid_region_arguments(self, kwargs: dict[str, object]) -> None:
        """Test seek rejects invalid region inputs before querying the backend."""

        parser = GFF3Parser()

        with pytest.raises(ValueError):
            parser.seek(**kwargs)  # type: ignore[arg-type]

    def test_sequential_mode_rejects_seek(self, tmp_path: Path) -> None:
        """Test sequential reader mode keeps seek disabled."""

        path = _write_gff3(tmp_path / "sample.gff3")
        parser = GFF3Parser(file_path=path, mode="sequential")

        with pytest.raises(RuntimeError, match="sequential"):
            parser.seek(seqid="chr1", start=10, end=20)


class TestGFF3ParserRecords:
    """Tests for GFF3Parser.records."""

    def test_records_returns_remaining_records(self, tmp_path: Path) -> None:
        """Test records returns all remaining records from the stream."""

        path = _write_gff3(tmp_path / "sample.gff3")
        parser = GFF3Parser(file_path=path, mode="sequential")

        _assert_gene_record(parser.next())

        assert [record.feature_type for record in parser.records()] == ["CDS"]

    def test_records_respects_limit(self, tmp_path: Path) -> None:
        """Test records stops once the limit is reached."""

        path = _write_gff3(tmp_path / "sample.gff3")
        parser = GFF3Parser(file_path=path, mode="sequential")

        assert [record.feature_type for record in parser.records(limit=1)] == ["gene"]
        assert [record.feature_type for record in parser.records()] == ["CDS"]

    def test_records_zero_limit_does_not_consume_record(self, tmp_path: Path) -> None:
        """Test zero limit leaves the next record available."""

        path = _write_gff3(tmp_path / "sample.gff3")
        parser = GFF3Parser(file_path=path, mode="sequential")

        assert parser.records(limit=0) == []
        _assert_gene_record(parser.next())


class TestGFF3ParserIteration:
    """Tests for GFF3Parser iteration."""

    def test_parser_is_its_own_iterator(self) -> None:
        """Test iter returns the parser instance."""

        parser = GFF3Parser()

        assert iter(parser) is parser

    def test_iteration_streams_records(self, tmp_path: Path) -> None:
        """Test iteration streams records from the parser."""

        path = _write_gff3(tmp_path / "sample.gff3")
        parser = GFF3Parser(file_path=path, mode="sequential")

        assert [record.feature_type for record in parser] == ["gene", "CDS"]

    def test_iteration_raises_stop_iteration(self, tmp_path: Path) -> None:
        """Test iteration raises StopIteration when no more records are available."""

        path = _write_gff3(tmp_path / "sample.gff3")
        parser = GFF3Parser(file_path=path, mode="sequential")

        with pytest.raises(StopIteration):
            while True:
                next(parser)


class TestGFF3ParserParsingTextAndFiles:
    """Tests for GFF3Parser parse helpers."""

    def test_parse_text_returns_interface_records(self) -> None:
        """Test parse_text exposes GFF3 rows as parser-interface records."""

        records = GFF3Parser.parse_text(SAMPLE)

        assert len(records) == 2
        _assert_gene_record(records[0])
        assert records[1].score == 4.5
        assert records[1].phase_value == 2

    def test_parse_text_respects_limit(self) -> None:
        """Test parse_text stops once the limit is reached."""

        records = GFF3Parser.parse_text(SAMPLE, limit=1)

        assert len(records) == 1
        _assert_gene_record(records[0])

    def test_parse_file_returns_interface_records(self, tmp_path: Path) -> None:
        """Test parse_file exposes GFF3 rows as parser-interface records."""

        path = _write_gff3(tmp_path / "sample.gff3")
        records = GFF3Parser.parse_file(path)

        assert len(records) == 2
        _assert_gene_record(records[0])

    def test_parse_file_respects_limit(self, tmp_path: Path) -> None:
        """Test parse_text stops once the limit is reached."""

        path = _write_gff3(tmp_path / "sample.gff3")
        records = GFF3Parser.parse_file(path, limit=1)

        assert len(records) == 1
        _assert_gene_record(records[0])
