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
"""Unit tests for the VCFVariationSource class."""

import pytest

from ensembl.vep.ensembl_vep_core import (  # noqa: F401,F403 pylint: disable=E0401,E0611
    VCFRecord
)
from ensembl.vep.parsers.vcf import VCFParser
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource
from ensembl.vep.sources.file.vcf.variation import VCFVariationSource

HEADER_LINES = [
    "##fileformat=VCFv4.3",
    '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">',
    '##FILTER=<ID=q10,Description="Quality below 10">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample0\tsample1",
]
CHR1_VAR_ROW1 = "chr1\t10\trs1;id4\tA\tC\t13\tPASS\tNS=2\tGT\t0/1\t1/1"
CHR1_VAR_ROW2 = "chr1\t20\t.\tG\tT\t.\tq10\tNS=2\tGT\t0/1\t0/0"
CHR2_VAR_ROW1 = "chr2\t100\trs3\tC\tA\t.\tPASS\tNS=2\tGT\t1/1\t0/1"
SAMPLE = "\n".join([*HEADER_LINES, CHR1_VAR_ROW1, CHR1_VAR_ROW2, CHR2_VAR_ROW1]) + "\n"

@pytest.fixture
def vcf_records(text: str = SAMPLE) -> list[VCFRecord]:
    """Create a list VCFRecord objects for a SNP for testing."""

    return VCFParser().parse_text(text)


class TestVCFVariationSourceInitialization:
    """Test VCFVariationSource initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test VCFVariationSource initialization without keywords raises error."""

        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = VCFVariationSource("sample.vcf")  # type: ignore[misc]

    def test_initialization_sets_core_attributes(self) -> None:
        """Test VCFVariationSource stores source and parser attributes."""

        source = VCFVariationSource(sequence_source=FastaSequenceSource())

        assert source.file_path is None
        assert source.adaptor.format == "vcf"  # pylint: disable=protected-access
        assert source.adaptor.file_path is None  # pylint: disable=protected-access
        assert source.adaptor.tabix_path == "tabix"  # pylint: disable=protected-access
        assert source.adaptor.reader_mode == "auto"  # pylint: disable=protected-access


class TestVCFVariationSourceMethods:
    """Test the VCFVariationSource methods."""

    def test_stable_id_from_vcf_record(self, vcf_records: list[VCFRecord]) -> None:
        """Test _stable_id_from_vcf_record method."""
        assert VCFVariationSource._stable_id_from_vcf_record(vcf_records[0]) == "rs1;id4"
        assert VCFVariationSource._stable_id_from_vcf_record(vcf_records[1]) is None
        assert VCFVariationSource._stable_id_from_vcf_record(vcf_records[2]) == "rs3"
