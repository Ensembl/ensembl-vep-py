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
"""Unit tests for the GTF parser interface."""

from pathlib import Path
import subprocess
from typing import Protocol

from attr import exceptions as attr_ex

import pytest

from ensembl.vep.parsers.gtf import GTFParser, GTFRecord

HEADER_LINES = ["#!genome-version GRCh38", "# comment 1"]
CHR1_GENE_ROW = 'chr1\tensembl\tgene\t10\t20\t.\t+\t.\tgene_id "gene:1"; gene_name "BRCA2";'
CHR1_TRANSCRIPT_ROW = (
    "chr1\tensembl\ttranscript\t12\t18\t4.5\t+\t.\t"
    'gene_id "gene:1"; transcript_id "transcript:1"; transcript_name "BRCA2-001";'
)
CHR1_EXON_ROW = (
    'chr1\tensembl\texon\t12\t18\t.\t+\t.\tgene_id "gene:1"; transcript_id "transcript:1"; exon_number "1";'
)
CHR1_CDS_ROW = 'chr1\tensembl\tCDS\t15\t18\t.\t+\t0\tgene_id "gene:1"; transcript_id "transcript:1";'
CHR2_GENE_ROW = 'chr2\tensembl\tgene\t100\t120\t.\t-\t.\tgene_id "gene:2"'
SAMPLE = "\n".join([HEADER_LINES[0], CHR1_GENE_ROW, CHR1_TRANSCRIPT_ROW, CHR1_CDS_ROW]) + "\n"
SAMPLE_WITH_COMMENT = "\n".join([*HEADER_LINES, CHR1_GENE_ROW, CHR1_TRANSCRIPT_ROW, CHR1_EXON_ROW]) + "\n"


# the default SAMPLE string creates a gtf file with the chr1 gene row, transcript row and CDS row
# whereas the mock tabix executable returns the chr1 gene row, transcript row and exon row for queries on chr1
def _write_gtf(path: Path, data: str = SAMPLE) -> Path:
    """Write a GTF fixture file."""

    path.write_text(data)
    return path


def _write_query_gtf(path: Path, data: str = SAMPLE_WITH_COMMENT) -> Path:
    """Write a query-mode GTF fixture with a placeholder tabix index."""

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
    """Write a mock tabix executable that returns deterministic GTF records."""

    return _write_mock_tabix(
        tmp_path,
        f"""
if [ "$1" = "--version" ]; then
    exit 0
fi

case "$2" in
    chr1)
        printf '%s\\n' '{CHR1_GENE_ROW}'
        printf '%s\\n' '{CHR1_TRANSCRIPT_ROW}'
        printf '%s\\n' '{CHR1_EXON_ROW}'
        ;;
    chr1:10-20)
        printf '%s\\n' '{CHR1_GENE_ROW}'
        printf '%s\\n' '{CHR1_TRANSCRIPT_ROW}'
        printf '%s\\n' '{CHR1_EXON_ROW}'
        ;;
    chr2:100-120)
        printf '%s\\n' '{CHR2_GENE_ROW}'
        ;;
esac
""",
    )


def _assert_gene_record(record: GTFRecord) -> None:
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
    assert record.raw_attributes == 'gene_id "gene:1"; gene_name "BRCA2";'
    assert record.attributes == {"gene_id": ("gene:1",), "gene_name": ("BRCA2",)}


class TestGTFParserInitialization:
    """Test GTFParser initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test GTFParser initialization without keywords raises error."""

        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = GTFParser("annotations.gtf")  # type: ignore[misc]

    def test_initialization_sets_core_attributes(self) -> None:
        """Test GTFParser sets parser attributes."""

        parser = GTFParser()

        assert parser.format == "gtf"
        assert parser.file_path is None
        assert parser.tabix_path == "tabix"
        assert parser.mode == "auto"
        assert parser.reader_mode == "auto"

    def test_initialization_with_file_path_opens_parser(self, tmp_path: Path) -> None:
        """Test file_path initialization opens the parser."""

        path = _write_gtf(tmp_path / "sample.gtf")
        parser = GTFParser(file_path=path)

        _assert_gene_record(parser.next())

    def test_mode_is_frozen(self) -> None:
        """Test attribute mode is settable on object creation and them immutable"""
        parser = GTFParser(mode="sequential")

        assert parser.mode == "sequential", "mode attribute set incorrectly"
        assert parser.reader_mode == "sequential", "reader_mode attribute set incorrectly"
        # check that mode is not settable after object creation
        with pytest.raises(attr_ex.FrozenAttributeError, match="can't set attribute"):
            parser.mode = "query"


class TestGTFParserMetadata:
    """Tests for GTFParser.metadata."""

    def test_metadata_returns_header_lines(self, tmp_path: Path) -> None:
        """Test metadata returns raw header lines."""

        path = _write_gtf(tmp_path / "sample.gtf", SAMPLE_WITH_COMMENT)
        parser = GTFParser(file_path=path, mode="sequential")

        assert parser.metadata == {"header": HEADER_LINES}

    def test_metadata_is_cached(self, tmp_path: Path) -> None:
        """Test metadata is cached after first access."""

        path = _write_gtf(tmp_path / "sample.gtf", SAMPLE_WITH_COMMENT)
        parser = GTFParser(file_path=path, mode="sequential")

        metadata = parser.metadata

        assert parser.metadata == metadata

    def test_open_refreshes_metadata_cache(self, tmp_path: Path) -> None:
        """Test open clears cached metadata for the previous file."""

        first_path = _write_gtf(
            tmp_path / "first.gtf",
            (
                "#!genome-version GRCh37\nchr1\tensembl\tgene\t10\t20\t.\t+\t.\t"
                'gene_id "gene:1"; gene_name "BRCA2"\n'
            ),
        )
        second_path = _write_gtf(
            tmp_path / "second.gtf",
            (
                "#!genome-version GRCh38\nchr1\tensembl\tgene\t10\t20\t.\t+\t.\t"
                'gene_id "gene:1"; gene_name "BRCA2"\n'
            ),
        )
        parser = GTFParser(file_path=first_path, mode="sequential")

        assert parser.metadata == {"header": ["#!genome-version GRCh37"]}
        assert parser.open(second_path) is parser
        assert parser.metadata == {"header": ["#!genome-version GRCh38"]}

    def test_metadata_returns_empty_dictionary_after_close(self, tmp_path: Path) -> None:
        """Test metadata is empty when the parser is closed."""

        path = _write_gtf(tmp_path / "sample.gtf", SAMPLE_WITH_COMMENT)
        parser = GTFParser(file_path=path, mode="sequential")

        assert parser.metadata == {"header": HEADER_LINES}

        parser.close()

        assert not parser.metadata

    def test_metadata_does_not_consume_sequential_stream(self, tmp_path: Path) -> None:
        """Test metadata reads do not advance the active sequential stream."""

        path = _write_gtf(tmp_path / "sample.gtf", SAMPLE_WITH_COMMENT)
        parser = GTFParser(file_path=path, mode="sequential")

        _assert_gene_record(parser.next())
        assert parser.metadata == {"header": HEADER_LINES}
        record = parser.next()
        assert record is not None
        assert record.feature_type == "transcript"

    def test_metadata_in_query_mode_reads_source_headers(self, tmp_path: Path) -> None:
        """Test query-mode metadata is read from the GTF source file."""

        path = _write_query_gtf(tmp_path / "sample.gtf.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = GTFParser(file_path=path, tabix_path=tabix)

        assert parser.reader_mode == "query"
        assert parser.metadata == {"header": HEADER_LINES}


class TestGTFParserOpenClose:
    """Tests for GTFParser open and close."""

    def test_open_returns_parser(self, tmp_path: Path) -> None:
        """Test open returns the parser instance."""

        path = _write_gtf(tmp_path / "sample.gtf")
        parser = GTFParser(mode="sequential")

        assert parser.open(path) is parser
        _assert_gene_record(parser.next())

    def test_open_overrides_mode(self, tmp_path: Path) -> None:
        """Test open overrides the parser mode."""

        path = _write_query_gtf(tmp_path / "sample.gtf.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = GTFParser(file_path=None, tabix_path=tabix)
        assert parser.reader_mode == "auto"

        assert parser.open(path, mode="sequential") is parser
        assert parser.reader_mode == "sequential"

    def test_open_auto_query_then_auto_seq_changes_mode(self, tmp_path: Path) -> None:
        """Test open overrides the parser mode."""

        path = _write_query_gtf(tmp_path / "sample.gtf.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = GTFParser(file_path=None, tabix_path=tabix)
        assert parser.reader_mode == "auto"

        assert parser.open(path) is parser
        assert parser.reader_mode == "query"

        path = _write_gtf(tmp_path / "sample.gtf")
        assert parser.open(path) is parser
        assert parser.reader_mode == "sequential"

    def test_close_closes_backend_reader(self, tmp_path: Path) -> None:
        """Test close prevents subsequent backend reads."""

        path = _write_gtf(tmp_path / "sample.gtf")
        parser = GTFParser(file_path=path, mode="sequential")

        parser.close()

        with pytest.raises(RuntimeError, match="not open"):
            parser.next()

    def test_context_manager_closes_parser(self, tmp_path: Path) -> None:
        """Test context manager returns and closes the parser."""

        path = _write_gtf(tmp_path / "sample.gtf")

        with GTFParser(file_path=path, mode="sequential") as parser:
            _assert_gene_record(parser.next())

        with pytest.raises(RuntimeError, match="not open"):
            parser.next()


class TestGTFParserSeek:
    """Tests for GTFParser.seek."""

    def test_gtf_parser_seek_primes_next(self, tmp_path: Path) -> None:
        """Test seek returns availability and next streams parsed records."""

        path = _write_query_gtf(tmp_path / "sample.gtf.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = GTFParser(file_path=path, tabix_path=tabix)

        assert parser.seek(seq_region_name="chr1", start=10, end=20) is True
        _assert_gene_record(parser.next())
        assert [record.feature_type for record in parser.records()] == ["transcript", "exon"]
        assert parser.next() is None

    def test_gtf_parser_seek_accepts_region_string(self, tmp_path: Path) -> None:
        """Test seek accepts an explicit region string."""

        path = _write_query_gtf(tmp_path / "sample.gtf.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = GTFParser(file_path=path, tabix_path=tabix)

        assert parser.seek("chr2:100-120") is True
        assert [record.seqid for record in parser.records()] == ["chr2"]

    def test_gtf_parser_seek_accepts_whole_sequence(self, tmp_path: Path) -> None:
        """Test seek accepts a seq_region_name-only region."""

        path = _write_query_gtf(tmp_path / "sample.gtf.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = GTFParser(file_path=path, tabix_path=tabix)

        assert parser.seek(seq_region_name="chr1") is True
        assert [record.feature_type for record in parser.records()] == ["gene", "transcript", "exon"]

    def test_gtf_parser_seek_empty_region(self, tmp_path: Path) -> None:
        """Test seek returns false when no records are available."""

        path = _write_query_gtf(tmp_path / "sample.gtf.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = GTFParser(file_path=path, tabix_path=tabix)

        assert parser.seek("chr1:30-40") is False
        assert parser.next() is None

    @pytest.mark.parametrize(
        "kwargs",
        [
            {},
            {"region": "chr1:10-20", "seq_region_name": "chr1"},
            {"seq_region_name": "chr1", "start": 0, "end": 20},
            {"seq_region_name": "chr1", "start": 10, "end": 0},
        ],
    )
    def test_gtf_parser_seek_rejects_invalid_region_arguments(self, kwargs: dict[str, object]) -> None:
        """Test seek rejects invalid region inputs before querying the backend."""

        parser = GTFParser()

        with pytest.raises(ValueError):
            parser.seek(**kwargs)  # type: ignore[arg-type]

    def test_sequential_mode_rejects_seek(self, tmp_path: Path) -> None:
        """Test sequential reader mode keeps seek disabled."""

        path = _write_gtf(tmp_path / "sample.gtf")
        parser = GTFParser(file_path=path, mode="sequential")

        with pytest.raises(RuntimeError, match="sequential"):
            parser.seek(seq_region_name="chr1", start=10, end=20)


class TestGTFParserRecords:
    """Tests for GTFParser.records."""

    def test_records_returns_remaining_records(self, tmp_path: Path) -> None:
        """Test records returns all remaining records from the stream."""

        path = _write_gtf(tmp_path / "sample.gtf")
        parser = GTFParser(file_path=path, mode="sequential")

        _assert_gene_record(parser.next())

        assert [record.feature_type for record in parser.records()] == ["transcript", "CDS"]

    def test_records_respects_limit(self, tmp_path: Path) -> None:
        """Test records stops once the limit is reached."""

        path = _write_gtf(tmp_path / "sample.gtf")
        parser = GTFParser(file_path=path, mode="sequential")

        assert [record.feature_type for record in parser.records(limit=1)] == ["gene"]
        assert [record.feature_type for record in parser.records()] == ["transcript", "CDS"]

    def test_records_zero_limit_does_not_consume_record(self, tmp_path: Path) -> None:
        """Test zero limit leaves the next record available."""

        path = _write_gtf(tmp_path / "sample.gtf")
        parser = GTFParser(file_path=path, mode="sequential")

        assert parser.records(limit=0) == []
        _assert_gene_record(parser.next())


class TestGTFParserIteration:
    """Tests for GTFParser iteration."""

    def test_parser_is_its_own_iterator(self) -> None:
        """Test iter returns the parser instance."""

        parser = GTFParser()

        assert iter(parser) is parser

    def test_iteration_streams_records(self, tmp_path: Path) -> None:
        """Test iteration streams records from the parser."""

        path = _write_gtf(tmp_path / "sample.gtf")
        parser = GTFParser(file_path=path, mode="sequential")

        assert [record.feature_type for record in parser] == ["gene", "transcript", "CDS"]

    def test_iteration_raises_stop_iteration(self, tmp_path: Path) -> None:
        """Test iteration raises StopIteration when no more records are available."""

        path = _write_gtf(tmp_path / "sample.gtf")
        parser = GTFParser(file_path=path, mode="sequential")

        with pytest.raises(StopIteration):
            while True:
                next(parser)


class TestGTFParserParsingTextAndFiles:
    """Tests for GTFParser parse helpers."""

    def test_parse_text_returns_interface_records(self) -> None:
        """Test parse_text exposes GTF rows as parser-interface records."""

        records = GTFParser.parse_text(SAMPLE)

        assert len(records) == 3
        _assert_gene_record(records[0])
        assert records[1].score == 4.5  # check score of transcript record
        assert records[2].phase_value == 0  # check phase of CDS record

    def test_parse_text_respects_limit(self) -> None:
        """Test parse_text stops once the limit is reached."""

        records = GTFParser.parse_text(SAMPLE, limit=1)

        assert len(records) == 1
        _assert_gene_record(records[0])

    def test_parse_file_returns_interface_records(self, tmp_path: Path) -> None:
        """Test parse_file exposes GTF rows as parser-interface records."""

        path = _write_gtf(tmp_path / "sample.gtf")
        records = GTFParser.parse_file(path)

        assert len(records) == 3
        _assert_gene_record(records[0])
        assert records[1].score == 4.5  # check score of transcript record
        assert records[2].phase_value == 0  # check phase of CDS record

    def test_parse_file_respects_limit(self, tmp_path: Path) -> None:
        """Test parse_text stops once the limit is reached."""

        path = _write_gtf(tmp_path / "sample.gtf")
        records = GTFParser.parse_file(path, limit=1)

        assert len(records) == 1
        _assert_gene_record(records[0])


class TestGTFRecordFirstAttribute:
    """Tests for GTFRecord.first_attribute."""

    def test_returns_first_value_for_known_key(self) -> None:
        """Test the first attribute value is returned for a single-valued key."""

        record = GTFParser.parse_text(SAMPLE, limit=1)[0]

        assert record.first_attribute("gene_id") == "gene:1"

    def test_returns_none_for_missing_key(self) -> None:
        """Test None is returned when the attribute key is absent."""

        record = GTFParser.parse_text(SAMPLE, limit=1)[0]

        assert record.first_attribute("transcript_id") is None


class TestGTFParserIndexSequenceNames:
    """Tests for GTFParser.index_sequence_names."""

    def test_unopened_parser_has_no_index_sequence_names(self) -> None:
        """Test an unopened parser reports no index sequence names."""

        assert GTFParser().index_sequence_names == []

    def test_sequential_reader_has_no_index_sequence_names(self, tmp_path: Path) -> None:
        """Test a sequentially-read file reports no index sequence names."""

        path = _write_gtf(tmp_path / "sample.gtf")
        parser = GTFParser(file_path=path, mode="sequential")

        assert parser.index_sequence_names == []

    def test_system_backend_has_no_index_sequence_names(self, tmp_path: Path) -> None:
        """Test the system tabix executable fallback reports no index sequence names."""

        path = _write_query_gtf(tmp_path / "sample.gtf.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = GTFParser(file_path=path, tabix_path=tabix)

        assert parser.index_sequence_names == []


class HtslibToolchain(Protocol):
    """Protocol for HTSlib toolchain fixtures."""

    bgzip: str
    tabix: str


def _create_tabix_gtf_cli(path: Path, htslib_toolchain: HtslibToolchain) -> Path:
    """Create and tabix-index a GTF file with a real HTSlib toolchain."""

    unzip_path = path.with_suffix(path.suffix.replace(".gz", ""))
    unzip_path.write_text(
        "\n".join([HEADER_LINES[0], CHR1_GENE_ROW, CHR1_TRANSCRIPT_ROW, CHR1_CDS_ROW, CHR2_GENE_ROW]) + "\n"
    )

    subprocess.run([htslib_toolchain.bgzip, "-f", str(unzip_path)], check=True)
    subprocess.run([htslib_toolchain.tabix, "-f", "-p", "gff", str(path)], check=True)

    return path


@pytest.mark.integration
class TestGTFParserIndexIntegration:
    """Integration tests against a real, in-process (noodles) tabix index."""

    def test_index_sequence_names_matches_real_index(
        self, tmp_path: Path, htslib_toolchain: HtslibToolchain
    ) -> None:
        """Test index_sequence_names reports the real index's declared contigs, in order."""

        path = _create_tabix_gtf_cli(tmp_path / "sample.gtf.gz", htslib_toolchain)
        parser = GTFParser(file_path=str(path), tabix_path="tabix-nonexistent-binary-for-test")

        assert parser.index_sequence_names == ["chr1", "chr2"]

    def test_seek_unknown_contig_returns_false_not_raises(
        self, tmp_path: Path, htslib_toolchain: HtslibToolchain
    ) -> None:
        """Test seeking a contig absent from a real index returns False instead of raising."""

        path = _create_tabix_gtf_cli(tmp_path / "sample.gtf.gz", htslib_toolchain)
        parser = GTFParser(file_path=str(path), tabix_path="tabix-nonexistent-binary-for-test")

        assert parser.seek(seq_region_name="chrZ", start=1, end=100) is False
        assert parser.next() is None

        # The parser stays usable for a subsequent, known contig.
        assert parser.seek(seq_region_name="chr1", start=1, end=100) is True
        assert parser.next() is not None
