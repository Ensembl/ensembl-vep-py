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

from ensembl.vep.parsers.tabix import TabixParser


def test_tabix_parser_seek_primes_next(tmp_path: Path) -> None:
    """TabixParser.seek returns availability and next streams raw rows."""

    path = tmp_path / "sample.gff3.gz"
    path.write_text("")

    tabix = tmp_path / "tabix"
    tabix.write_text(
        """#!/bin/sh
printf '%s\\n' 'chr1\ttabix\tgene\t10\t20\t.\t+\t.\tID=gene:1'
printf '%s\\n' 'chr1\ttabix\texon\t12\t14\t.\t+\t.\tID=exon:1'
"""
    )
    tabix.chmod(0o755)

    parser = TabixParser(path, tabix_path=str(tabix))

    assert parser.seek(seqid="chr1", start=10, end=20) is True
    assert parser.next() == "chr1\ttabix\tgene\t10\t20\t.\t+\t.\tID=gene:1"
    assert parser.records() == ["chr1\ttabix\texon\t12\t14\t.\t+\t.\tID=exon:1"]
    assert parser.next() is None


def test_tabix_parser_seek_empty_region(tmp_path: Path) -> None:
    """TabixParser.seek returns false when no rows are available."""

    path = tmp_path / "sample.gff3.gz"
    path.write_text("")

    tabix = tmp_path / "tabix"
    tabix.write_text("#!/bin/sh\n")
    tabix.chmod(0o755)

    parser = TabixParser(path, tabix_path=str(tabix))

    assert parser.seek("chr1:30-40") is False
    assert parser.next() is None
