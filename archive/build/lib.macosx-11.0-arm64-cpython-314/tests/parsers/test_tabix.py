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
import subprocess
from typing import Protocol
import pytest

from ensembl.vep.parsers.tabix import TabixParser

HEADER_LINES = [
    "##gff-version 3",
    "#seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes",
]
CHR1_GENE1_ROW = "chr1\ttabix\tgene\t10\t20\t.\t+\t.\tID=gene:1"
CHR1_EXON1_ROW = "chr1\ttabix\texon\t12\t14\t.\t+\t.\tID=exon:1"
CHR1_EXON2_ROW = "chr1\ttabix\texon\t15\t18\t.\t+\t0\tID=exon:2"
CHR1_GENE2_ROW = "chr1\ttabix\tgene\t50\t60\t.\t+\t.\tID=gene:2"
CHR2_GENE_ROW = "chr2\ttabix\tgene\t100\t120\t.\t-\t.\tID=gene:2"
BED_HEADER_ROW = "#chrom\tstart\tend\tname"
BED_CHR1_FEATURE1_ROW = "chr1\t9\t20\tfeature1"
BED_CHR1_FEATURE2_ROW = "chr1\t11\t14\tfeature2"
BED_CHR1_FEATURE3_ROW = "chr1\t14\t18\tfeature3"
BED_CHR1_FEATURE4_ROW = "chr1\t49\t60\tfeature4"
BED_CHR2_FEATURE1_ROW = "chr2\t99\t120\tfeature1"
BED_ZERO_POS_ROW = "chr1\t0\t0\tfeature1"
BED_SAME_POS_ROW = "chr1\t99\t99\tfeature1"
GFF3_ZERO_POS_ROW = "chr1\ttabix\tgene\t0\t0\t.\t+\t.\tID=gene:0"
GFF3_SAME_POS_ROW = "chr1\ttabix\tgene\t100\t100\t.\t+\t.\tID=gene:100"


class HtslibToolchain(Protocol):
    """Protocol for HTSlib toolchain fixtures."""

    bgzip: str
    tabix: str


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
        printf '%s\\n' '{CHR1_GENE1_ROW}'
        printf '%s\\n' '{CHR1_EXON1_ROW}'
        printf '%s\\n' '{CHR1_EXON2_ROW}'
        ;;
    chr1:10-20)
        printf '%s\\n' '{CHR1_GENE1_ROW}'
        printf '%s\\n' '{CHR1_EXON1_ROW}'
        ;;
    chr2:100-120)
        printf '%s\\n' '{CHR2_GENE_ROW}'
        ;;
esac
""",
    )


def _create_tabix_gff3_cli(
    path: Path,
    htslib_toolchain: HtslibToolchain,
    lines: list[str] | None = None,
) -> Path:
    """Create and tabix-index a GFF3 file with a real HTSlib toolchain."""

    unzip_path = path.with_suffix(path.suffix.replace(".gz", ""))

    fixture_lines = lines or [
        HEADER_LINES[0],
        HEADER_LINES[1],
        CHR1_GENE1_ROW,
        CHR1_EXON1_ROW,
        CHR1_EXON2_ROW,
        CHR1_GENE2_ROW,
        CHR2_GENE_ROW,
    ]
    unzip_path.write_text("\n".join(fixture_lines) + "\n")

    subprocess.run([htslib_toolchain.bgzip, "-f", str(unzip_path)], check=True)
    subprocess.run([htslib_toolchain.tabix, "-f", "-p", "gff", str(path)], check=True)

    return path


def _create_tabix_bed_cli(
    path: Path,
    htslib_toolchain: HtslibToolchain,
    lines: list[str] | None = None,
) -> Path:
    """Create and tabix-index a BED file with a real HTSlib toolchain."""

    unzip_path = path.with_suffix(path.suffix.replace(".gz", ""))

    fixture_lines = lines or [
        BED_HEADER_ROW,
        BED_CHR1_FEATURE1_ROW,
        BED_CHR1_FEATURE2_ROW,
        BED_CHR1_FEATURE3_ROW,
        BED_CHR1_FEATURE4_ROW,
        BED_CHR2_FEATURE1_ROW,
    ]
    unzip_path.write_text("\n".join(fixture_lines) + "\n")

    subprocess.run([htslib_toolchain.bgzip, "-f", str(unzip_path)], check=True)
    subprocess.run([htslib_toolchain.tabix, "-f", "-p", "bed", str(path)], check=True)

    return path


def _open_integration_parser(
    path: Path,
    backend: str,
    htslib_toolchain: HtslibToolchain,
) -> TabixParser:
    """Open a parser for a selected integration backend."""

    tabix_path = "tabix-nonexistent-binary-for-test" if backend == "noodles" else htslib_toolchain.tabix

    return TabixParser(file_path=path, tabix_path=tabix_path)


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

        assert parser.seek(seq_region_name="chr1", start=10, end=20) is True


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
        assert parser.seek(seq_region_name="chr2", start=100, end=120) is True
        assert parser.next() == CHR2_GENE_ROW

    def test_open_with_mode_logs_warning(self, tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
        """Test assigning an invalid mode logs a warning."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(tabix_path=str(tabix))

        with caplog.at_level("DEBUG"):
            parser.open(path, mode="mode")
            assert "The 'mode' argument is not used by TabixParser and will be ignored." in caplog.text

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
            assert parser.seek(seq_region_name="chr2", start=100, end=120) is True
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

        assert parser.seek(seq_region_name="chr1", start=10, end=20) is True
        assert parser.next() == CHR1_GENE1_ROW
        assert parser.records() == [CHR1_EXON1_ROW]
        assert parser.next() is None

    def test_tabix_parser_seek_accepts_region_string(self, tmp_path: Path) -> None:
        """Test seek accepts an explicit region string."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        assert parser.seek("chr2:100-120") is True
        assert parser.records() == [CHR2_GENE_ROW]

    def test_tabix_parser_seek_accepts_whole_sequence(self, tmp_path: Path) -> None:
        """Test seek accepts a seq_region_name-only region."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        assert parser.seek(seq_region_name="chr1") is True
        assert parser.records() == [CHR1_GENE1_ROW, CHR1_EXON1_ROW, CHR1_EXON2_ROW]

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
            {"region": "chr1:10-20", "seq_region_name": "chr1"},
            {"seq_region_name": "chr1", "start": 0, "end": 20},
            {"seq_region_name": "chr1", "start": 10, "end": 0},
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

        parser.seek(seq_region_name="chr1")

        assert parser.next() == CHR1_GENE1_ROW
        assert parser.records() == [CHR1_EXON1_ROW, CHR1_EXON2_ROW]

    def test_records_respects_limit(self, tmp_path: Path) -> None:
        """Test records stops once the limit is reached."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        parser.seek(seq_region_name="chr1")

        assert parser.records(limit=2) == [CHR1_GENE1_ROW, CHR1_EXON1_ROW]
        assert parser.records() == [CHR1_EXON2_ROW]

    def test_records_zero_limit_does_not_consume_buffered_row(self, tmp_path: Path) -> None:
        """Test zero limit leaves the primed row available."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        parser.seek(seq_region_name="chr1", start=10, end=20)

        assert parser.records(limit=0) == []
        assert parser.next() == CHR1_GENE1_ROW


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

        parser.seek(seq_region_name="chr1", start=10, end=20)

        assert list(parser) == [CHR1_GENE1_ROW, CHR1_EXON1_ROW]

        record_list = [CHR1_GENE1_ROW, CHR1_EXON1_ROW]
        for idx, record in enumerate(parser):
            assert record == record_list[idx]

    def test_iteration_raises_stop_iteration(self, tmp_path: Path) -> None:
        """Test iteration raises StopIteration when no more rows are available."""

        path = _write_placeholder_file(tmp_path / "sample.gff3.gz")
        tabix = _write_query_mock_tabix(tmp_path)
        parser = TabixParser(file_path=path, tabix_path=str(tabix))

        parser.seek(seq_region_name="chr1", start=10, end=20)

        with pytest.raises(StopIteration):
            while True:
                next(parser)


@pytest.mark.integration
class TestTabixParserIntegration:
    """Integration tests for TabixParser across backends."""

    @pytest.mark.parametrize(
        ("dataset", "expected_header", "expected_first", "expected_rest"),
        [
            ("gff3", HEADER_LINES, CHR1_GENE1_ROW, [CHR1_EXON1_ROW, CHR1_EXON2_ROW]),
            (
                "bed",
                [BED_HEADER_ROW],
                BED_CHR1_FEATURE1_ROW,
                [BED_CHR1_FEATURE2_ROW, BED_CHR1_FEATURE3_ROW],
            ),
        ],
        ids=["gff3", "bed"],
    )
    @pytest.mark.parametrize("backend", ["noodles", "system"], ids=["noodles", "system"])
    def test_integration_seek_next_behavior(
        self,
        tmp_path: Path,
        htslib_toolchain: HtslibToolchain,
        dataset: str,
        expected_header: list[str],
        expected_first: str,
        expected_rest: list[str],
        backend: str,
    ) -> None:
        """Test both backends can read default GFF3 and BED fixtures."""

        if dataset == "gff3":
            path = _create_tabix_gff3_cli(tmp_path / f"sample-{dataset}-{backend}.gff3.gz", htslib_toolchain)
        else:
            path = _create_tabix_bed_cli(tmp_path / f"sample-{dataset}-{backend}.bed.gz", htslib_toolchain)

        parser = _open_integration_parser(path, backend, htslib_toolchain)

        assert parser.metadata == {"header": expected_header}
        assert parser.seek(seq_region_name="chr1", start=10, end=20) is True
        assert parser.next() == expected_first
        assert parser.records() == expected_rest
        assert parser.next() is None

    @pytest.mark.parametrize(
        ("dataset", "expected_header", "expected_first", "expected_rest"),
        [
            ("gff3", HEADER_LINES, CHR1_GENE1_ROW, [CHR1_EXON1_ROW, CHR1_EXON2_ROW, CHR1_GENE2_ROW]),
            (
                "bed",
                [BED_HEADER_ROW],
                BED_CHR1_FEATURE1_ROW,
                [BED_CHR1_FEATURE2_ROW, BED_CHR1_FEATURE3_ROW, BED_CHR1_FEATURE4_ROW],
            ),
        ],
        ids=["gff3", "bed"],
    )
    @pytest.mark.parametrize("backend", ["noodles", "system"], ids=["noodles", "system"])
    def test_integration_query_full_chromosome(
        self,
        tmp_path: Path,
        htslib_toolchain: HtslibToolchain,
        dataset: str,
        expected_header: list[str],
        expected_first: str,
        expected_rest: list[str],
        backend: str,
    ) -> None:
        """Test both backends can read default GFF3 and BED fixtures."""

        if dataset == "gff3":
            path = _create_tabix_gff3_cli(tmp_path / f"sample-{dataset}-{backend}.gff3.gz", htslib_toolchain)
        else:
            path = _create_tabix_bed_cli(tmp_path / f"sample-{dataset}-{backend}.bed.gz", htslib_toolchain)

        parser = _open_integration_parser(path, backend, htslib_toolchain)

        assert parser.metadata == {"header": expected_header}
        assert parser.seek(seq_region_name="chr1") is True
        assert parser.next() == expected_first
        assert parser.records() == expected_rest
        assert parser.next() is None

    @pytest.mark.parametrize(
        ("dataset", "expected_header", "expected_first", "expected_rest"),
        [
            ("gff3", HEADER_LINES, CHR1_GENE1_ROW, [CHR1_EXON1_ROW, CHR1_EXON2_ROW]),
            (
                "bed",
                [BED_HEADER_ROW],
                BED_CHR1_FEATURE1_ROW,
                [BED_CHR1_FEATURE2_ROW, BED_CHR1_FEATURE3_ROW],
            ),
        ],
        ids=["gff3", "bed"],
    )
    @pytest.mark.parametrize("backend", ["noodles", "system"], ids=["noodles", "system"])
    def test_integration_query_subregion(
        self,
        tmp_path: Path,
        htslib_toolchain: HtslibToolchain,
        dataset: str,
        expected_header: list[str],
        expected_first: str,
        expected_rest: list[str],
        backend: str,
    ) -> None:
        """Test both backends can read default GFF3 and BED fixtures."""

        if dataset == "gff3":
            path = _create_tabix_gff3_cli(tmp_path / f"sample-{dataset}-{backend}.gff3.gz", htslib_toolchain)
        else:
            path = _create_tabix_bed_cli(tmp_path / f"sample-{dataset}-{backend}.bed.gz", htslib_toolchain)

        parser = _open_integration_parser(path, backend, htslib_toolchain)

        assert parser.metadata == {"header": expected_header}
        assert parser.seek(seq_region_name="chr1", start=12, end=20) is True
        assert parser.next() == expected_first
        assert parser.records() == expected_rest
        assert parser.next() is None

    @pytest.mark.parametrize(
        ("dataset", "expected_header", "expected_first", "expected_rest"),
        [
            ("gff3", HEADER_LINES, CHR1_GENE1_ROW, [CHR1_EXON1_ROW, CHR1_EXON2_ROW, CHR1_GENE2_ROW]),
            (
                "bed",
                [BED_HEADER_ROW],
                BED_CHR1_FEATURE1_ROW,
                [BED_CHR1_FEATURE2_ROW, BED_CHR1_FEATURE3_ROW, BED_CHR1_FEATURE4_ROW],
            ),
        ],
        ids=["gff3", "bed"],
    )
    @pytest.mark.parametrize("backend", ["noodles", "system"], ids=["noodles", "system"])
    def test_integration_query_without_end(
        self,
        tmp_path: Path,
        htslib_toolchain: HtslibToolchain,
        dataset: str,
        expected_header: list[str],
        expected_first: str,
        expected_rest: list[str],
        backend: str,
    ) -> None:
        """Test both backends can read default GFF3 and BED fixtures."""

        if dataset == "gff3":
            path = _create_tabix_gff3_cli(tmp_path / f"sample-{dataset}-{backend}.gff3.gz", htslib_toolchain)
        else:
            path = _create_tabix_bed_cli(tmp_path / f"sample-{dataset}-{backend}.bed.gz", htslib_toolchain)

        parser = _open_integration_parser(path, backend, htslib_toolchain)

        assert parser.metadata == {"header": expected_header}
        assert parser.seek(seq_region_name="chr1", start=12) is True
        assert parser.next() == expected_first
        assert parser.records() == expected_rest
        assert parser.next() is None

    @pytest.mark.parametrize(("dataset"), ["gff3", "bed"], ids=["gff3", "bed"])
    @pytest.mark.parametrize("backend", ["noodles", "system"], ids=["noodles", "system"])
    def test_integration_query_start_greater_than_end(
        self,
        tmp_path: Path,
        htslib_toolchain: HtslibToolchain,
        dataset: str,
        backend: str,
    ) -> None:
        """Test query with start greater than end returns no records."""

        if dataset == "gff3":
            path = _create_tabix_gff3_cli(tmp_path / f"gff3-before-{backend}.gff3.gz", htslib_toolchain)
        else:
            path = _create_tabix_bed_cli(tmp_path / f"bed-before-{backend}.bed.gz", htslib_toolchain)

        parser = _open_integration_parser(path, backend, htslib_toolchain)

        assert parser.seek("chr2:200-100") is False
        assert parser.records() == []

    @pytest.mark.parametrize(("dataset"), ["gff3", "bed"], ids=["gff3", "bed"])
    @pytest.mark.parametrize("backend", ["noodles", "system"], ids=["noodles", "system"])
    def test_integration_query_just_before(
        self,
        tmp_path: Path,
        htslib_toolchain: HtslibToolchain,
        dataset: str,
        backend: str,
    ) -> None:
        """Test query just before the chr2 feature returns no records."""

        if dataset == "gff3":
            path = _create_tabix_gff3_cli(tmp_path / f"gff3-before-{backend}.gff3.gz", htslib_toolchain)
        else:
            path = _create_tabix_bed_cli(tmp_path / f"bed-before-{backend}.bed.gz", htslib_toolchain)

        parser = _open_integration_parser(path, backend, htslib_toolchain)

        assert parser.seek("chr2:99-99") is False
        assert parser.records() == []

    @pytest.mark.parametrize(
        ("dataset", "expected_row"),
        [
            ("gff3", CHR2_GENE_ROW),
            ("bed", BED_CHR2_FEATURE1_ROW),
        ],
        ids=["gff3", "bed"],
    )
    @pytest.mark.parametrize("backend", ["noodles", "system"], ids=["noodles", "system"])
    def test_integration_query_exactly_on_start(
        self,
        tmp_path: Path,
        htslib_toolchain: HtslibToolchain,
        dataset: str,
        expected_row: str,
        backend: str,
    ) -> None:
        """Test query exactly on chr2 feature start returns matching record."""

        if dataset == "gff3":
            path = _create_tabix_gff3_cli(tmp_path / f"gff3-exact-start-{backend}.gff3.gz", htslib_toolchain)
        else:
            path = _create_tabix_bed_cli(tmp_path / f"bed-exact-start-{backend}.bed.gz", htslib_toolchain)

        parser = _open_integration_parser(path, backend, htslib_toolchain)

        assert parser.seek("chr2:100-100") is True
        assert parser.records() == [expected_row]

    @pytest.mark.parametrize(
        ("dataset", "expected_row"),
        [
            ("gff3", CHR2_GENE_ROW),
            ("bed", BED_CHR2_FEATURE1_ROW),
        ],
        ids=["gff3", "bed"],
    )
    @pytest.mark.parametrize("backend", ["noodles", "system"], ids=["noodles", "system"])
    def test_integration_query_exactly_on_end(
        self,
        tmp_path: Path,
        htslib_toolchain: HtslibToolchain,
        dataset: str,
        expected_row: str,
        backend: str,
    ) -> None:
        """Test query exactly on chr2 feature end returns matching record."""

        if dataset == "gff3":
            path = _create_tabix_gff3_cli(tmp_path / f"gff3-exact-end-{backend}.gff3.gz", htslib_toolchain)
        else:
            path = _create_tabix_bed_cli(tmp_path / f"bed-exact-end-{backend}.bed.gz", htslib_toolchain)

        parser = _open_integration_parser(path, backend, htslib_toolchain)

        assert parser.seek("chr2:120-120") is True
        assert parser.records() == [expected_row]

    @pytest.mark.parametrize("dataset", ["gff3", "bed"], ids=["gff3", "bed"])
    @pytest.mark.parametrize("backend", ["noodles", "system"], ids=["noodles", "system"])
    def test_integration_query_just_after(
        self,
        tmp_path: Path,
        htslib_toolchain: HtslibToolchain,
        dataset: str,
        backend: str,
    ) -> None:
        """Test query just after the chr2 feature returns no records."""

        if dataset == "gff3":
            path = _create_tabix_gff3_cli(tmp_path / f"gff3-after-{backend}.gff3.gz", htslib_toolchain)
        else:
            path = _create_tabix_bed_cli(tmp_path / f"bed-after-{backend}.bed.gz", htslib_toolchain)

        parser = _open_integration_parser(path, backend, htslib_toolchain)

        assert parser.seek("chr2:121-121") is False
        assert parser.records() == []

    @pytest.mark.parametrize(
        ("dataset", "lines"),
        [
            (
                "gff3",
                [HEADER_LINES[0], HEADER_LINES[1], GFF3_ZERO_POS_ROW],
            ),
            (
                "bed",
                [BED_ZERO_POS_ROW],
            ),
        ],
        ids=["gff3", "bed"],
    )
    @pytest.mark.parametrize("backend", ["noodles", "system"], ids=["noodles", "system"])
    def test_integration_zero_position_queries(
        self,
        tmp_path: Path,
        htslib_toolchain: HtslibToolchain,
        dataset: str,
        lines: list[str],
        backend: str,
    ) -> None:
        """Test zero-position query behavior for GFF3 and BED fixtures.

        A readable index is always queried in process, so both backends follow noodles here.
        """

        if dataset == "gff3":
            path = _create_tabix_gff3_cli(
                tmp_path / f"{dataset}-{backend}.gff3.gz",
                htslib_toolchain,
                lines,
            )
        else:
            path = _create_tabix_bed_cli(
                tmp_path / f"{dataset}-{backend}.bed.gz",
                htslib_toolchain,
                lines,
            )
        parser = _open_integration_parser(path, backend, htslib_toolchain)

        if dataset == "gff3":
            with pytest.raises(ValueError):
                parser.seek("chr1:1-1")
            return

        assert parser.seek("chr1:1-1") is True
        assert parser.records() == [BED_ZERO_POS_ROW]

    @pytest.mark.parametrize(
        ("dataset", "lines"),
        [
            (
                "gff3",
                [HEADER_LINES[0], HEADER_LINES[1], GFF3_SAME_POS_ROW],
            ),
            (
                "bed",
                [BED_SAME_POS_ROW],
            ),
        ],
        ids=["gff3", "bed"],
    )
    @pytest.mark.parametrize("backend", ["noodles", "system"], ids=["noodles", "system"])
    def test_integration_same_position_queries(
        self,
        tmp_path: Path,
        htslib_toolchain: HtslibToolchain,
        dataset: str,
        lines: list[str],
        backend: str,
    ) -> None:
        """Test same-position query behavior for GFF3 and BED fixtures."""

        if dataset == "gff3":
            path = _create_tabix_gff3_cli(
                tmp_path / f"{dataset}-{backend}.gff3.gz",
                htslib_toolchain,
                lines,
            )
        else:
            path = _create_tabix_bed_cli(
                tmp_path / f"{dataset}-{backend}.bed.gz",
                htslib_toolchain,
                lines,
            )
        parser = _open_integration_parser(path, backend, htslib_toolchain)

        expected_row = [GFF3_SAME_POS_ROW] if dataset == "gff3" else []
        assert parser.seek("chr1:100-100") is bool(expected_row)
        assert parser.records() == expected_row

        expected_row = [GFF3_SAME_POS_ROW] if dataset == "gff3" else [BED_SAME_POS_ROW]
        assert parser.seek("chr1:99-100") is bool(expected_row)
        assert parser.records() == expected_row
