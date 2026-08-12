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
CHR1_VAR_ROW2 = "chr1\t20\t.\tGAT\tG\t.\tq10\tNS=2\tGT\t0/1\t0/0"
CHR2_VAR_ROW1 = "chr2\t100\trs3\tC\tATC\t.\tPASS\tNS=2\tGT\t1/1\t0/1"
SAMPLE = "\n".join([*HEADER_LINES, CHR1_VAR_ROW1, CHR1_VAR_ROW2, CHR2_VAR_ROW1]) + "\n"

CHR1_SV_ROW1 = "chr1\t50\tsv1\tA\t<DEL>\t.\tPASS\tEND=80;SVTYPE=DEL;\tGT\t0/1\t0/0"
CHR1_SV_ROW2 = "chr1\t14\t.\tT\t<INS>\t.\t.\tIMPRECISE;SVLEN=20;CILEN=-10,10;CIPOS=-10,10\tGT\t0/1\t1/1"
CHR1_SV_ROW3 = "chr1\t101\tsv3\tT\t<DUP>\t.\tPASS\tSVLEN=20;CIPOS=-5,2;CIEND=-3,4\tGT\t0/1\t0/0"
CHR1_SV_ROW4 = "chr1\t201\tsv4\tT\t<INV>\t.\tPASS\tSVLEN=30;CIPOS=-5,2;CIEND=-3,4\tGT\t0/1\t0/0"
CHR1_SV_ROW5 = "chr1\t301\tsv5\tT\t<CNV>\t.\tPASS\tSVLEN=40;CIPOS=-5,2;CIEND=-3,4\tGT\t0/1\t0/0"
SV_SAMPLE = "\n".join([*HEADER_LINES, CHR1_SV_ROW1, CHR1_SV_ROW2, CHR1_SV_ROW3, CHR1_SV_ROW4, CHR1_SV_ROW5]) + "\n"

def create_vcf_records_for_testing(text: str) -> list[VCFRecord]:
    """Create a list of VCFRecord objects for testing."""

    return VCFParser().parse_text(text)

@pytest.fixture
def vcf_records(text: str = SAMPLE) -> list[VCFRecord]:
    """Create a list of SNV VCFRecord objects for testing."""

    return create_vcf_records_for_testing(text)

@pytest.fixture
def sv_vcf_records(text: str = SV_SAMPLE) -> list[VCFRecord]:
    """Create a list VCFRecord objects for a SNP for testing."""

    return create_vcf_records_for_testing(text)

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

    def test_end_from_vcf_record_snv(self, vcf_records: list[VCFRecord]) -> None:
        """Test _end_from_vcf_record method."""
        source = VCFVariationSource(sequence_source=FastaSequenceSource())

        assert source._end_from_vcf_record(vcf_records[0]) == 10
        assert source._end_from_vcf_record(vcf_records[1]) == 22
        assert source._end_from_vcf_record(vcf_records[2]) == 100

    def test_end_from_vcf_record_symbolic(self, sv_vcf_records: list[VCFRecord]) -> None:
        """Test _end_from_vcf_record method."""
        source = VCFVariationSource(sequence_source=FastaSequenceSource())

        assert source._end_from_vcf_record(sv_vcf_records[0], is_symbolic=True) == 80
        assert source._end_from_vcf_record(sv_vcf_records[1], is_symbolic=True) == 33
        assert source._end_from_vcf_record(sv_vcf_records[2], is_symbolic=True) == 120
        assert source._end_from_vcf_record(sv_vcf_records[3], is_symbolic=True) == 230 
        assert source._end_from_vcf_record(sv_vcf_records[4], is_symbolic=True) == 340
