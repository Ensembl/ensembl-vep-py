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
"""Unit tests for the VCF parser interface."""

from pathlib import Path
from attr import exceptions as attr_ex
import pytest

from ensembl.vep.parsers.vcf import VCFParser, VCFRecord

HEADER_LINES = [
    "##fileformat=VCFv4.3",
    '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">',
    '##FILTER=<ID=q10,Description="Quality below 10">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample0\tsample1",
]
CHR1_VAR_ROW1 = "chr1\t10\trs1\tA\tC\t13\tPASS\tNS=2\tGT\t0/1\t1/1"
CHR1_VAR_ROW2 = "chr1\t20\trs2\tG\tT\t.\tq10\tNS=2\tGT\t0/1\t0/0"
CHR2_VAR_ROW1 = "chr2\t100\trs3\tC\tA\t.\tPASS\tNS=2\tGT\t1/1\t0/1"

SAMPLE = "\n".join([*HEADER_LINES, CHR1_VAR_ROW1, CHR1_VAR_ROW2, CHR2_VAR_ROW1]) + "\n"


def _write_vcf(path: Path, data: str = SAMPLE) -> Path:
    """Write a VCF fixture file."""

    path.write_text(data)
    return path


def _write_query_vcf(path: Path, data: str = SAMPLE) -> Path:
    """Write a query-mode VCF fixture with a placeholder tabix index."""

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
    """Write a mock tabix executable that returns deterministic VCF records."""

    return _write_mock_tabix(
        tmp_path,
        f"""
if [ "$1" = "--version" ]; then
    exit 0
fi

case "$2" in
    chr1)
        printf '%s\\n' '{CHR1_VAR_ROW1}'
        printf '%s\\n' '{CHR1_VAR_ROW2}'
        ;;
    chr1:10-15)
        printf '%s\\n' '{CHR1_VAR_ROW1}'
        ;;
    chr2:100-120)
        printf '%s\\n' '{CHR2_VAR_ROW1}'
        ;;
esac
""",
    )


def _assert_first_vcf_record(record: VCFRecord) -> None:
    """Assert a parsed record matches the first chr1 fixture."""

    assert record is not None
    assert record.chrom == "chr1"
    assert record.pos == 10
    assert record.id == ["rs1"]
    assert record.reference == "A"
    assert record.alt == ["C"]
    assert record.qual == 13.0
    assert record.filter == ["PASS"]
    assert record.info == {"NS": "2"}
    assert record.samples == [{"GT": "0/1"}, {"GT": "1/1"}]


def _assert_second_vcf_record(record: VCFRecord) -> None:
    """Assert a parsed record matches the second chr1 fixture."""

    assert record is not None
    assert record.chrom == "chr1"
    assert record.pos == 20
    assert record.id == ["rs2"]
    assert record.reference == "G"
    assert record.alt == ["T"]
    assert record.qual is None
    assert record.filter == ["q10"]
    assert record.info == {"NS": "2"}
    assert record.samples == [{"GT": "0/1"}, {"GT": "0/0"}]


class TestVCFParserInitialization:
    """Test VCFParser initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test VCFParser initialization without keywords raises error."""

        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = VCFParser("annotations.vcf")  # type: ignore[misc]

    def test_initialization_sets_core_attributes(self) -> None:
        """Test VCFParser sets parser attributes."""

        parser = VCFParser()

        assert parser.format == "vcf"
        assert parser.file_path is None
        assert parser.tabix_path == "tabix"
        assert parser.mode == "auto"
        assert parser.reader_mode == "auto"

    def test_initialization_with_file_path_opens_parser(self, tmp_path: Path) -> None:
        """Test file_path initialization opens the parser."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        parser = VCFParser(file_path=path)
        record = parser.next()
        _assert_first_vcf_record(record)

    def test_mode_is_frozen(self) -> None:
        """Test attribute mode is settable on object creation and them immutable"""
        parser = VCFParser(mode="sequential")

        assert parser.mode == "sequential", "mode attribute set incorrectly"
        assert parser.reader_mode == "sequential", "reader_mode attribute set incorrectly"
        # check that mode is not settable after object creation
        with pytest.raises(attr_ex.FrozenAttributeError, match="can't set attribute"):
            parser.mode = "query"


class TestVCFParserMetadata:
    """Tests for VCFParser.metadata."""

    @pytest.fixture
    def sample_metadata(self) -> dict:
        """Create a metadata dict for comparison"""
        return {
            "file_format": "VCFv4.3",
            "info": {
                "NS": {
                    "type": "Integer",
                    "idx": "",
                    "number": "1",
                    "description": "Number of samples with data",
                    "other_fields": "",
                }
            },
            "filter": {"q10": {"idx": "", "other_fields": "", "description": "Quality below 10"}},
            "format": {
                "GT": {
                    "number": "1",
                    "description": "Genotype",
                    "type": "String",
                    "other_fields": "",
                    "idx": "",
                }
            },
            "alt": {},
            "contig": {},
            "samples": ["sample0", "sample1"],
            "other": {},
        }

    def test_metadata_returns_header_info(self, tmp_path: Path, sample_metadata: dict) -> None:
        """Test metadata returns header lines."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        parser = VCFParser(file_path=path, mode="sequential")

        assert parser.metadata == sample_metadata

    def test_metadata_is_cached(self, tmp_path: Path) -> None:
        """Test metadata is cached after first access."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        parser = VCFParser(file_path=path, mode="sequential")

        metadata = parser.metadata

        assert parser.metadata == metadata

    def test_open_refreshes_metadata_cache(self, tmp_path: Path, sample_metadata: dict) -> None:
        """Test open clears cached metadata for the previous file."""

        first_path = _write_vcf(tmp_path / "sample-1.vcf", SAMPLE)

        HEADER_TMP = HEADER_LINES.copy()
        HEADER_TMP[0] = "##fileformat=VCFv4.5"
        SAMPLE2 = "\n".join([*HEADER_TMP, CHR1_VAR_ROW1]) + "\n"
        second_path = _write_vcf(tmp_path / "sample-2.vcf", SAMPLE2)

        parser = VCFParser(file_path=first_path, mode="sequential")

        assert parser.metadata == sample_metadata
        assert parser.open(second_path) is parser
        assert parser.metadata["file_format"] == "VCFv4.5"

    def test_metadata_returns_empty_dictionary_after_close(self, tmp_path: Path, sample_metadata: dict) -> None:
        """Test metadata is empty when the parser is closed."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        parser = VCFParser(file_path=path, mode="sequential")

        assert parser.metadata == sample_metadata
        parser.close()

        assert not parser.metadata

    def test_metadata_does_not_consume_sequential_stream(self, tmp_path: Path, sample_metadata: dict) -> None:
        """Test metadata reads do not advance the active sequential stream."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        parser = VCFParser(file_path=path, mode="sequential")

        _assert_first_vcf_record(parser.next())
        assert parser.metadata == sample_metadata
        record = parser.next()
        assert record is not None
        _assert_second_vcf_record(record)

    def test_metadata_in_query_mode_reads_source_headers(self, tmp_path: Path, sample_metadata: dict) -> None:
        """Test query-mode metadata is read from the VCF source file."""

        path = _write_query_vcf(tmp_path / "sample.vcf.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = VCFParser(file_path=path, tabix_path=tabix)

        assert parser.reader_mode == "query"
        assert parser.metadata == sample_metadata


class TestVCFParserOpenClose:
    """Tests for VCFParser open and close."""

    def test_open_returns_parser(self, tmp_path: Path) -> None:
        """Test open returns the parser instance."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        parser = VCFParser(mode="sequential")

        assert parser.open(path) is parser
        _assert_first_vcf_record(parser.next())

    def test_open_overrides_mode(self, tmp_path: Path) -> None:
        """Test open overrides the parser mode."""

        path = _write_query_vcf(tmp_path / "sample.vcf.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = VCFParser(file_path=None, tabix_path=tabix)
        assert parser.reader_mode == "auto"

        assert parser.open(path, mode="sequential") is parser
        assert parser.reader_mode == "sequential"

    def test_open_auto_query_then_auto_seq_changes_mode(self, tmp_path: Path) -> None:
        """Test open overrides the parser mode."""

        path = _write_query_vcf(tmp_path / "sample.vcf.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = VCFParser(file_path=None, tabix_path=tabix)
        assert parser.reader_mode == "auto"

        assert parser.open(path) is parser
        assert parser.reader_mode == "query"

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        assert parser.open(path) is parser
        assert parser.reader_mode == "sequential"

    def test_close_closes_backend_reader(self, tmp_path: Path) -> None:
        """Test close prevents subsequent backend reads."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        parser = VCFParser(file_path=path, mode="sequential")

        parser.close()

        with pytest.raises(RuntimeError, match="not open"):
            parser.next()

    def test_context_manager_closes_parser(self, tmp_path: Path) -> None:
        """Test context manager returns and closes the parser."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)

        with VCFParser(file_path=path, mode="sequential") as parser:
            _assert_first_vcf_record(parser.next())

        with pytest.raises(RuntimeError, match="not open"):
            parser.next()


class TestVCFParserSeek:
    """Tests for VCFParser.seek."""

    def test_vcf_parser_seek_primes_next(self, tmp_path: Path) -> None:
        """Test seek returns availability and next streams parsed records."""

        path = _write_query_vcf(tmp_path / "sample.vcf.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = VCFParser(file_path=path, tabix_path=tabix)

        assert parser.seek(seq_region_name="chr1", start=10, end=15) is True
        _assert_first_vcf_record(parser.next())
        assert parser.next() is None

    def test_vcf_parser_seek_accepts_region_string(self, tmp_path: Path) -> None:
        """Test seek accepts an explicit region string."""

        path = _write_query_vcf(tmp_path / "sample.vcf.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = VCFParser(file_path=path, tabix_path=tabix)

        assert parser.seek("chr2:100-120") is True
        assert [record.chrom for record in parser.records()] == ["chr2"]

    def test_vcf_parser_seek_accepts_whole_sequence(self, tmp_path: Path) -> None:
        """Test seek accepts a seq_region_name-only region."""

        path = _write_query_vcf(tmp_path / "sample.vcf.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = VCFParser(file_path=path, tabix_path=tabix)

        assert parser.seek(seq_region_name="chr1") is True
        assert [record.reference for record in parser.records()] == ["A", "G"]

    def test_vcf_parser_seek_empty_region(self, tmp_path: Path) -> None:
        """Test seek returns false when no records are available."""

        path = _write_query_vcf(tmp_path / "sample.vcf.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = VCFParser(file_path=path, tabix_path=tabix)

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
    def test_vcf_parser_seek_rejects_invalid_region_arguments(self, kwargs: dict[str, object]) -> None:
        """Test seek rejects invalid region inputs before querying the backend."""

        parser = VCFParser()

        with pytest.raises(ValueError):
            parser.seek(**kwargs)  # type: ignore[arg-type]

    def test_sequential_mode_rejects_seek(self, tmp_path: Path) -> None:
        """Test sequential reader mode keeps seek disabled."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        parser = VCFParser(file_path=path, mode="sequential")

        with pytest.raises(RuntimeError, match="sequential"):
            parser.seek(seq_region_name="chr1", start=10, end=20)


class TestVCFParserRecords:
    """Tests for VCFParser.records."""

    def test_records_returns_remaining_records(self, tmp_path: Path) -> None:
        """Test records returns all remaining records from the stream."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        parser = VCFParser(file_path=path, mode="sequential")

        _assert_first_vcf_record(parser.next())

        assert [record.chrom for record in parser.records()] == ["chr1", "chr2"]

    def test_records_respects_limit(self, tmp_path: Path) -> None:
        """Test records stops once the limit is reached."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        parser = VCFParser(file_path=path, mode="sequential")

        assert [record.reference for record in parser.records(limit=1)] == ["A"]
        assert [record.reference for record in parser.records()] == ["G", "C"]

    def test_records_zero_limit_does_not_consume_record(self, tmp_path: Path) -> None:
        """Test zero limit leaves the next record available."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        parser = VCFParser(file_path=path, mode="sequential")

        assert parser.records(limit=0) == []
        _assert_first_vcf_record(parser.next())


class TestVCFParserIteration:
    """Tests for VCFParser iteration."""

    def test_parser_is_its_own_iterator(self) -> None:
        """Test iter returns the parser instance."""

        parser = VCFParser()

        assert iter(parser) is parser

    def test_iteration_streams_records(self, tmp_path: Path) -> None:
        """Test iteration streams records from the parser."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        parser = VCFParser(file_path=path, mode="sequential")

        assert [record.reference for record in parser] == ["A", "G", "C"]

    def test_iteration_raises_stop_iteration(self, tmp_path: Path) -> None:
        """Test iteration raises StopIteration when no more records are available."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        parser = VCFParser(file_path=path, mode="sequential")

        with pytest.raises(StopIteration):
            while True:
                next(parser)


class TestVCFParserParsingTextAndFiles:
    """Tests for VCFParser parse helpers."""

    def test_parse_text_returns_interface_records(self) -> None:
        """Test parse_text exposes VCF rows as parser-interface records."""

        records = VCFParser.parse_text(SAMPLE)

        assert len(records) == 3
        _assert_first_vcf_record(records[0])
        _assert_second_vcf_record(records[1])

    def test_parse_text_respects_limit(self) -> None:
        """Test parse_text stops once the limit is reached."""

        records = VCFParser.parse_text(SAMPLE, limit=1)

        assert len(records) == 1
        _assert_first_vcf_record(records[0])

    def test_parse_file_returns_interface_records(self, tmp_path: Path) -> None:
        """Test parse_file exposes VCF rows as parser-interface records."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        records = VCFParser.parse_file(path)

        assert len(records) == 3
        _assert_first_vcf_record(records[0])

    def test_parse_file_respects_limit(self, tmp_path: Path) -> None:
        """Test parse_text stops once the limit is reached."""

        path = _write_vcf(tmp_path / "sample.vcf", SAMPLE)
        records = VCFParser.parse_file(path, limit=1)

        assert len(records) == 1
        _assert_first_vcf_record(records[0])
