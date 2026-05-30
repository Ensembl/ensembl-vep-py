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
"""Unit tests for the tabix parser interface."""

from pathlib import Path
import shutil
import subprocess
import pytest

from ensembl.vep.parsers.tabix import TabixParser


HEADER_LINES = [
    "##gff-version 3",
    "#seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes",
]
CHR1_GENE_ROW = "chr1\ttabix\tgene\t10\t20\t.\t+\t.\tID=gene:1"
CHR1_EXON_ROW = "chr1\ttabix\texon\t12\t14\t.\t+\t.\tID=exon:1"
CHR1_CDS_ROW = "chr1\ttabix\tCDS\t15\t18\t.\t+\t0\tID=cds:1"
CHR2_GENE_ROW = "chr2\ttabix\tgene\t100\t120\t.\t-\t.\tID=gene:2"


def _write_placeholder_file(path: Path) -> Path:
    """Write a placeholder data file for mocked tabix backends."""

    path.write_text("")
    return path


def _write_mock_tabix(tmp_path: Path, body: str) -> Path:
    """Write a mock tabix executable."""

    tabix = tmp_path / "tabix"
    tabix.write_text(f"#!/bin/sh\n{body}\n")
    tabix.chmod(0o755)
    return tabix


def _write_query_mock_tabix(tmp_path: Path) -> Path:
    """Write a mock tabix executable that returns deterministic records."""

    return _write_mock_tabix(
        tmp_path,
        f"""
if [ "$1" = "-H" ]; then
    printf '%s\\n' '{HEADER_LINES[0]}'
    printf '%s\\n' '{HEADER_LINES[1]}'
    exit 0
fi

case "$2" in
    chr1)
        printf '%s\\n' '{CHR1_GENE_ROW}'
        printf '%s\\n' '{CHR1_EXON_ROW}'
        printf '%s\\n' '{CHR1_CDS_ROW}'
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

def _create_tabix_gff3_cli(path: Path) -> Path:
    unzip_path = path.with_suffix(path.suffix.replace(".gz", ""))

    unzip_path.write_text(f"""{HEADER_LINES[0]}
{HEADER_LINES[1]}
{CHR1_GENE_ROW}
{CHR1_EXON_ROW}
{CHR1_CDS_ROW}
{CHR2_GENE_ROW}
""")

    subprocess.run(["bgzip", "-f", str(unzip_path)], check=True)
    subprocess.run(["tabix", "-p", "gff", str(path)], check=True)

    return path


class TestTabixParserInitialization:
    """Test TabixParser initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test TabixParser initialization without keywords raises error."""

        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = TabixParser("annotations.gff3.gz")  # type: ignore[misc]

    def test_initialization_sets_core_attributes(self) -> None:
        """Test TabixParser sets parser attributes."""

        parser = TabixParser()

        assert parser.format == "tabix"
        assert parser.file_path is None
        assert parser.tabix_path == "tabix"

    def test_initialization_with_file_path_opens_parser(self, tmp_path: Path) -> None:
        """Test file_path initialization opens the parser."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)

        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        assert parser.seek(seqid="chr1", start=10, end=20) is True


class TestTabixParserMetadata:
    """Tests for TabixParser.metadata."""

    def test_metadata_returns_header_lines(self, tmp_path: Path) -> None:
        """Test metadata returns header lines from the backend."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        assert parser.metadata == {"header": HEADER_LINES}

    def test_metadata_is_cached(self, tmp_path: Path) -> None:
        """Test metadata is cached after first access."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        metadata = parser.metadata

        assert parser.metadata is metadata

    def test_open_refreshes_metadata_cache(self, tmp_path: Path) -> None:
        """Test open clears cached metadata for the previous file."""

        first_path = _write_placeholder_file(tmp_path / "first.gff3.gz")
        second_path = _write_placeholder_file(tmp_path / "second.gff3.gz")
        tabix = _write_mock_tabix(
            tmp_path,
            """
if [ "$1" = "-H" ]; then
    case "$2" in
        *first.gff3.gz)
            printf '%s\\n' '##first-header'
            ;;
        *second.gff3.gz)
            printf '%s\\n' '##second-header'
            ;;
    esac
    exit 0
fi
""",
        )
        parser = TabixParser(file_path=first_path, tabix_path=str(tabix))

        assert parser.metadata == {"header": ["##first-header"]}
        assert parser.open(second_path) is parser
        assert parser.metadata == {"header": ["##second-header"]}

    def test_metadata_returns_empty_dictionary_after_close(self, tmp_path: Path) -> None:
        """Test metadata is empty when the parser is closed."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        assert parser.metadata == {"header": HEADER_LINES}

        parser.close()

        assert not parser.metadata

    def test_metadata_raises_backend_errors(self, tmp_path: Path) -> None:
        """Test metadata raises backend errors from tabix header retrieval."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_mock_tabix(
            tmp_path,
            """
if [ "$1" = "-H" ]; then
    exit 2
fi
""",
        )
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        with pytest.raises(RuntimeError, match="tabix -H failed"):
            _ = parser.metadata


class TestTabixParserOpenClose:
    """Tests for TabixParser open and close."""

    def test_open_returns_parser(self, tmp_path: Path) -> None:
        """Test open returns the parser instance."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(tabix_path=str(tabix))

        assert parser.open(path) is parser
        assert parser.seek(seqid="chr2", start=100, end=120) is True
        assert parser.next() == CHR2_GENE_ROW

    def test_close_closes_backend_reader(self, tmp_path: Path) -> None:
        """Test close prevents subsequent backend reads."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        parser.close()

        with pytest.raises(RuntimeError, match="not open"):
            parser.next()

    def test_context_manager_closes_parser(self, tmp_path: Path) -> None:
        """Test context manager returns and closes the parser."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)

        with TabixParser(file_path=path, tabix_path=str(tabix)) as parser:
            assert parser.seek(seqid="chr2", start=100, end=120) is True
            assert parser.next() == CHR2_GENE_ROW

        with pytest.raises(RuntimeError, match="not open"):
            parser.next()


class TestTabixParserSeek:
    """Tests for TabixParser.seek."""

    def test_tabix_parser_seek_primes_next(self, tmp_path: Path) -> None:
        """Test seek returns availability and next streams raw rows."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        assert parser.seek(seqid="chr1", start=10, end=20) is True
        assert parser.next() == CHR1_GENE_ROW
        assert parser.records() == [CHR1_EXON_ROW]
        assert parser.next() is None

    def test_tabix_parser_seek_accepts_region_string(self, tmp_path: Path) -> None:
        """Test seek accepts an explicit region string."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        assert parser.seek("chr2:100-120") is True
        assert parser.records() == [CHR2_GENE_ROW]

    def test_tabix_parser_seek_accepts_whole_sequence(self, tmp_path: Path) -> None:
        """Test seek accepts a seqid-only region."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        assert parser.seek(seqid="chr1") is True
        assert parser.records() == [CHR1_GENE_ROW, CHR1_EXON_ROW, CHR1_CDS_ROW]

    def test_tabix_parser_seek_empty_region(self, tmp_path: Path) -> None:
        """Test seek returns false when no rows are available."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        assert parser.seek("chr1:30-40") is False
        assert parser.next() is None

    @pytest.mark.parametrize(
        "kwargs",
        [
            {},
            {"region": "chr1:10-20", "seqid": "chr1"},
            {"seqid": "chr1", "start": 10},
            {"seqid": "chr1", "start": 0, "end": 20},
            {"seqid": "chr1", "start": 20, "end": 10},
        ],
    )
    def test_tabix_parser_seek_rejects_invalid_region_arguments(self, kwargs: dict[str, object]) -> None:
        """Test seek rejects invalid region inputs before querying the backend."""

        parser = TabixParser()

        with pytest.raises(ValueError):
            parser.seek(**kwargs)  # type: ignore[arg-type]


class TestTabixParserRecords:
    """Tests for TabixParser.records."""

    def test_records_returns_remaining_rows(self, tmp_path: Path) -> None:
        """Test records returns all remaining rows from the active query."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        parser.seek(seqid="chr1")

        assert parser.next() == CHR1_GENE_ROW
        assert parser.records() == [CHR1_EXON_ROW, CHR1_CDS_ROW]

    def test_records_respects_limit(self, tmp_path: Path) -> None:
        """Test records stops once the limit is reached."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        parser.seek(seqid="chr1")

        assert parser.records(limit=2) == [CHR1_GENE_ROW, CHR1_EXON_ROW]
        assert parser.records() == [CHR1_CDS_ROW]

    def test_records_zero_limit_does_not_consume_buffered_row(self, tmp_path: Path) -> None:
        """Test zero limit leaves the primed row available."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        parser.seek(seqid="chr1", start=10, end=20)

        assert parser.records(limit=0) == []
        assert parser.next() == CHR1_GENE_ROW


class TestTabixParserIteration:
    """Tests for TabixParser iteration."""

    def test_parser_is_its_own_iterator(self) -> None:
        """Test iter returns the parser instance."""

        parser = TabixParser()

        assert iter(parser) is parser

    def test_iteration_streams_active_query(self, tmp_path: Path) -> None:
        """Test iteration streams rows from the active query."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        parser.seek(seqid="chr1", start=10, end=20)

        assert list(parser) == [CHR1_GENE_ROW, CHR1_EXON_ROW]

        record_list = [CHR1_GENE_ROW, CHR1_EXON_ROW]
        for idx, record in enumerate(parser):
            assert record == record_list[idx]

    def test_iteration_raises_stop_iteration(self, tmp_path: Path) -> None:
        """Test iteration raises StopIteration when no more rows are available."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        parser.seek(seqid="chr1", start=10, end=20)

        with pytest.raises(StopIteration):
            while True:
                next(parser)


class TestTabixParserNoodlesBackend:
    """Tests for TabixParser with the noodles backend."""

    def test_noodles_backend_reads_records(self, tmp_path: Path) -> None:
        """Test TabixParser can read records with the noodles backend."""

        path = _create_tabix_gff3_cli(tmp_path / "sample.gff3.gz")
        parser = TabixParser(file_path=path, tabix_path=None)

        assert parser.metadata == {"header": HEADER_LINES}
        assert parser.seek(seqid="chr1", start=10, end=20) is True
        assert parser.next() == CHR1_GENE_ROW
        assert parser.records() == [CHR1_EXON_ROW, CHR1_CDS_ROW]
        assert parser.next() is None


class TestTabixParserIntegration:
    """Integration tests for TabixParser with a real tabix executable."""

    def test_integration_with_real_tabix(self, tmp_path: Path) -> None:
        """Test TabixParser can read from a real tabix executable."""

        path = _create_tabix_gff3_cli(tmp_path / "sample.gff3.gz")
        tabix = shutil.which("tabix")

        if tabix is None:
            pytest.skip("tabix executable not found in PATH")

        parser = TabixParser(file_path=path, tabix_path=tabix)

        assert parser.metadata == {"header": HEADER_LINES}
        assert parser.seek(seqid="chr1", start=10, end=20) is True
        assert parser.next() == CHR1_GENE_ROW
        assert parser.records() == [CHR1_EXON_ROW, CHR1_CDS_ROW]
        assert parser.next() is None
