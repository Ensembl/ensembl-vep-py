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
"""Unit tests for the GFF parser interface."""

from pathlib import Path

import pytest

from ensembl.vep.parsers.gff import GFFParser


SAMPLE = """\
##gff-version 3
chr1\tensembl\tgene\t10\t20\t.\t+\t.\tID=gene:1;Name=BRCA2;Alias=a,b
chr1\tensembl\tCDS\t12\t18\t4.5\t-\t2\tID=cds:1;Parent=gene:1
##FASTA
>chr1
ACGT
"""


def test_parse_text_returns_interface_records() -> None:
    """GFF rows are exposed as parser-interface records."""

    records = GFFParser.parse_text(SAMPLE)

    assert len(records) == 2

    record = records[0]
    assert record.seqid == "chr1"
    assert record.source == "ensembl"
    assert record.feature_type == "gene"
    assert record.start == 10
    assert record.end == 20
    assert record.score is None
    assert record.strand == "+"
    assert record.phase is None
    assert record.raw_attributes == "ID=gene:1;Name=BRCA2;Alias=a,b"
    assert record.attributes == {"ID": ("gene:1",), "Name": ("BRCA2",), "Alias": ("a", "b")}

    assert records[1].score == 4.5
    assert records[1].phase == 2


def test_parser_streams_file(tmp_path: Path) -> None:
    """Parser.next streams records until the FASTA directive."""

    path = tmp_path / "sample.gff3"
    path.write_text(SAMPLE)

    parser = GFFParser(path)

    first = parser.next()
    second = parser.next()

    assert first is not None
    assert second is not None
    assert first.feature_type == "gene"
    assert second.feature_type == "CDS"
    assert parser.next() is None


def test_parser_next_does_not_require_tabix(tmp_path: Path) -> None:
    """Parser.next streams the whole file without opening tabix."""

    path = tmp_path / "sample.gff3"
    path.write_text(SAMPLE)

    parser = GFFParser(path, tabix_path=str(tmp_path / "missing-tabix"))

    assert [record.feature_type for record in parser] == ["gene", "CDS"]


def test_parser_can_force_sequential_reader(tmp_path: Path) -> None:
    """Sequential reader mode keeps seek disabled."""

    path = tmp_path / "sample.gff3"
    path.write_text(SAMPLE)

    parser = GFFParser(path, tabix_path=str(tmp_path / "missing-tabix"), reader="sequential")

    assert parser.reader == "sequential"
    assert [record.feature_type for record in parser] == ["gene", "CDS"]

    with pytest.raises(RuntimeError, match="sequential"):
        parser.seek(seqid="chr1", start=10, end=20)


def test_records_honours_limit(tmp_path: Path) -> None:
    """Parser.records returns a bounded number of rows."""

    path = tmp_path / "sample.gff3"
    path.write_text(SAMPLE)

    parser = GFFParser(path)
    records = parser.records(limit=1)

    assert len(records) == 1
    assert records[0].feature_type == "gene"


def test_seek_uses_tabix_and_parses_region_output(tmp_path: Path) -> None:
    """Parser.seek updates the active stream for next()."""

    path = tmp_path / "sample.gff3.gz"
    path.write_text(SAMPLE)

    tabix = tmp_path / "tabix"
    tabix.write_text(
        """#!/bin/sh
printf '%s\\n' 'chr1\ttabix\tgene\t10\t20\t.\t+\t.\tID=gene:1;Name=BRCA2'
"""
    )
    tabix.chmod(0o755)

    parser = GFFParser(path, tabix_path=str(tabix))
    has_records = parser.seek(seqid="chr1", start=10, end=20)
    record = parser.next()

    assert has_records is True
    assert record is not None
    assert record.seqid == "chr1"
    assert record.source == "tabix"
    assert record.feature_type == "gene"
    assert record.start == 10
    assert record.end == 20
    assert record.attributes == {"ID": ("gene:1",), "Name": ("BRCA2",)}
    assert parser.next() is None


def test_parser_can_force_tabix_reader(tmp_path: Path) -> None:
    """Tabix reader mode can be selected explicitly."""

    path = tmp_path / "sample.gff3.gz"
    path.write_text(SAMPLE)

    tabix = tmp_path / "tabix"
    tabix.write_text(
        """#!/bin/sh
printf '%s\\n' 'chr1\ttabix\tgene\t10\t20\t.\t+\t.\tID=gene:1'
"""
    )
    tabix.chmod(0o755)

    parser = GFFParser(path, tabix_path=str(tabix), reader="tabix")

    assert parser.reader == "tabix"
    assert parser.next() is None
    assert parser.seek(seqid="chr1", start=10, end=20) is True
    assert parser.next() is not None


def test_parser_set_reader_switches_reader_mode(tmp_path: Path) -> None:
    """Reader mode can be switched after construction."""

    path = tmp_path / "sample.gff3.gz"
    path.write_text(SAMPLE)

    tabix = tmp_path / "tabix"
    tabix.write_text(
        """#!/bin/sh
printf '%s\\n' 'chr1\ttabix\tgene\t10\t20\t.\t+\t.\tID=gene:1'
"""
    )
    tabix.chmod(0o755)

    parser = GFFParser(path, tabix_path=str(tabix))

    assert parser.reader == "auto"
    assert parser.set_reader("tabix") is parser
    assert parser.reader == "tabix"
    assert parser.seek(seqid="chr1", start=10, end=20) is True
    assert parser.next() is not None


def test_seek_returns_false_for_empty_region(tmp_path: Path) -> None:
    """Parser.seek returns false when tabix yields no rows."""

    path = tmp_path / "sample.gff3.gz"
    path.write_text(SAMPLE)

    tabix = tmp_path / "tabix"
    tabix.write_text("#!/bin/sh\n")
    tabix.chmod(0o755)

    parser = GFFParser(path, tabix_path=str(tabix))

    assert parser.seek(seqid="chr1", start=30, end=40) is False
    assert parser.next() is None
