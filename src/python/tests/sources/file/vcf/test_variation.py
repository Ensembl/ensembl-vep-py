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

import logging
import pytest

from ensembl.vep.ensembl_vep_core import VCFRecord  # noqa: F401,F403 pylint: disable=E0401,E0611
from ensembl.vep.common.enums import Strand
from ensembl.vep.models.variation.variation import Variation, VariationAllele
from ensembl.vep.models.variation.structural_variation import StructuralVariation, StructuralVariationAllele
from ensembl.vep.parsers.vcf import VCFParser
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource
from ensembl.vep.sources.file.vcf.variation import VCFVariationSource

logger = logging.getLogger(__name__)


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

CHR1_SV_DEL_SVTYPE = "chr1\t50\tsv1\tA\t<DEL>\t.\tPASS\tEND=80;SVTYPE=DEL;\tGT\t0/1\t0/0"
CHR1_SV_INS_IMPRECISE = "chr1\t14\t.\tT\t<INS>\t.\t.\tIMPRECISE;SVLEN=20;CILEN=-10,10;CIPOS=-10,10\tGT\t0/1\t1/1"
CHR1_SV_DUP_SVLEN = "chr1\t101\tsv3\tT\t<DUP>\t.\tPASS\tSVLEN=20;CIPOS=-5,2;CIEND=-3,4\tGT\t0/1\t0/0"
CHR1_SV_INV = "chr1\t201\tsv4\tT\t<INV>\t.\tPASS\tSVLEN=30;CIPOS=-5,2;CIEND=-3,4\tGT\t0/1\t0/0"
CHR1_SV_CNV = "chr1\t301\tsv5\tT\t<CNV>\t.\tPASS\tSVLEN=40;CIPOS=-5,2;CIEND=-3,4\tGT\t0/1\t0/0"
CHR1_SV_BND1 = "chrA\t2\t.\tT\t<DEL>\t.\t.\tSVLEN=2;SVCLAIM=D;EVENT=DEL_split_bp_cn\tGT\t0/1\t1/1"
CHR1_SV_BND2 = "chrA\t2\tdelbp1\tT\tT[chrA:5[\t.\t.\tMATEID=delbp2;EVENT=DEL_split_bp_cn\tGT\t0/1\t1/1"
CHR1_SV_BND3 = "chrA\t2\tdelbp2\tA\t]chrA:2]A\t.\t.\tMATEID=delbp1;EVENT=DEL_split_bp_cn\tGT\t0/1\t1/1"
CHR1_SV_BND4 = "chrA\t14\t.\tG\t.CCCCCCG\t.\t.\tEVENT=single_breakend\tGT\t0/1\t1/1"

SV_SAMPLE = "\n".join([
    *HEADER_LINES, CHR1_SV_DEL_SVTYPE, CHR1_SV_INS_IMPRECISE, CHR1_SV_DUP_SVLEN, CHR1_SV_INV, CHR1_SV_CNV,
    CHR1_SV_BND1, CHR1_SV_BND2, CHR1_SV_BND3, CHR1_SV_BND4]) + "\n"

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

class DummySequenceSource:
    """Minimal sequence source fixture for slice-length assertions."""

    def __init__(self, lengths: dict[str, int], sequences: dict[str, str] | None = None) -> None:
        self.lengths = lengths
        self.sequences = sequences or {}
        self.closed = False

    def get_sequence_by_start_end_strand(  # pragma: no cover - helper
        self,
        seq_region_name: str,
        _start: int,
        _end: int,
        _strand: Strand,
    ) -> str | None:
        """Return the sequence for a requested region."""

        return self.sequences.get(seq_region_name)

    def get_seq_region_length(self, seq_region_name: str) -> int | None:
        """Return the length for a requested sequence region."""

        return self.lengths.get(seq_region_name)

    def close(self) -> None:
        """Close the sequence source."""

        self.closed = True


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

    def test_end_and_svlen_for_sva_from_vcf_record(
        self, sv_vcf_records: list[VCFRecord], caplog: pytest.LogCaptureFixture
    ) -> None:
        """Test _end_and_svlen_for_sva_from_vcf_record method."""
        source = VCFVariationSource(sequence_source=FastaSequenceSource())

        # test for incorrect allele_index
        with caplog.at_level(logging.WARNING):
            results = source._end_and_svlen_for_sva_from_vcf_record(sv_vcf_records[0], allele_index=1)
            assert (
                "Allele index, 1, is too large, sv1" in caplog.text
            ), "Incorrect warning in _end_and_svlen_for_sva_from_vcf_record"
        assert results == (None, None)

        assert source._end_and_svlen_for_sva_from_vcf_record(sv_vcf_records[0], allele_index=0) == (80, 30)
        assert source._end_and_svlen_for_sva_from_vcf_record(sv_vcf_records[1], allele_index=0) == (14, 20)
        assert source._end_and_svlen_for_sva_from_vcf_record(sv_vcf_records[2], allele_index=0) == (121, 20)
        assert source._end_and_svlen_for_sva_from_vcf_record(sv_vcf_records[3], allele_index=0) == (231, 30)
        assert source._end_and_svlen_for_sva_from_vcf_record(sv_vcf_records[3], allele_index=1) == (201, None)
        assert source._end_and_svlen_for_sva_from_vcf_record(sv_vcf_records[4], allele_index=0) == (341, 40)
        assert source._end_and_svlen_for_sva_from_vcf_record(sv_vcf_records[4], allele_index=1) == (351, 50)
        assert source._end_and_svlen_for_sva_from_vcf_record(sv_vcf_records[6], allele_index=0) == (
            None,
            None,
        )

    def test_int_info_value(self, sv_vcf_records: list[VCFRecord]) -> None:
        """Test _int_info_value method."""

        source = VCFVariationSource(sequence_source=FastaSequenceSource())

        assert source._int_info_value(sv_vcf_records[0], "END") == 80
        assert source._int_info_value(sv_vcf_records[5], "END") == None
        assert source._int_info_value(sv_vcf_records[3], "SVLEN") == None

    def test_is_snv(self, vcf_records: list[VCFRecord], sv_vcf_records: list[VCFRecord]) -> None:
        """Test _is_snv method."""
        source = VCFVariationSource(sequence_source=FastaSequenceSource())

        for record in vcf_records:
            assert source._is_snv(record)

        for record in sv_vcf_records:
            assert not source._is_snv(record)

    def test_has_svtype(
        self,
        vcf_records: list[VCFRecord],
        sv_vcf_records: list[VCFRecord]
    ) -> None:
        """Test _has_svtype method."""
        source = VCFVariationSource(sequence_source=FastaSequenceSource())

        for record in vcf_records:
            assert not source._has_svtype(record)

        assert source._has_svtype(sv_vcf_records[0])
        assert not source._has_svtype(sv_vcf_records[1])

    def test_is_symbolic(
        self,
        vcf_records: list[VCFRecord],
        sv_vcf_records: list[VCFRecord]
    ) -> None:
        """Test _is_symbolic method."""
        source = VCFVariationSource(sequence_source=FastaSequenceSource())

        for record in vcf_records:
            assert not source._is_symbolic(record)

        for record in sv_vcf_records[0:6]:
            assert source._is_symbolic(record)
        for record in sv_vcf_records[6:]:
            assert not source._is_symbolic(record)

    def test_is_breakend(
        self,
        vcf_records: list[VCFRecord]
    ) -> None:
        """Test _is_breakend method."""
        SV_LINES = "\n".join([*HEADER_LINES, CHR1_SV_DEL_SVTYPE, CHR1_SV_BND1, CHR1_SV_BND2, CHR1_SV_BND3, CHR1_SV_BND4]) + "\n"
        sv_vcf_records = create_vcf_records_for_testing(SV_LINES)

        source = VCFVariationSource(sequence_source=FastaSequenceSource())

        for record in vcf_records:
            assert not source._is_breakend(record)

        ## NEED TO THINK ABOUT BREAKEND RECORDS COMPOSED OF MULTIPLE LINES
        assert not source._is_breakend(sv_vcf_records[0])
        assert not source._is_breakend(sv_vcf_records[1])
        assert source._is_breakend(sv_vcf_records[2])
        assert source._is_breakend(sv_vcf_records[3])
        assert source._is_breakend(sv_vcf_records[4])

    def test_is_symbolic_insertion(self) -> None:
        """Test _is_symbolic_insertion method."""
        source = VCFVariationSource(sequence_source=FastaSequenceSource())

        assert source._is_symbolic_insertion("INS")
        assert source._is_symbolic_insertion("<INS>")
        assert not source._is_symbolic_insertion("<DEL>")

    def test_is_structural_variation_record(
        self, vcf_records: list[VCFRecord], sv_vcf_records: list[VCFRecord]
    ) -> None:
        """Test _is_structural_variation_record method."""
        source = VCFVariationSource(sequence_source=FastaSequenceSource())

        for record in vcf_records:
            assert not source._is_structural_variation_record(record)

        for record in sv_vcf_records:
            assert source._is_structural_variation_record(record)

    def test_records_to_variation(self, vcf_records: list[VCFRecord]) -> None:
        """Test _records_to_variation method."""
        vcf_var_source = VCFVariationSource(
            sequence_source=DummySequenceSource({"chr1": 32}, {"chr1": "ACGTACGTACGTACGTACGTACGTACGTACGT"})
        )

        var_feats = vcf_var_source._records_to_features(records=vcf_records)

        assert len(var_feats) == 3
        assert isinstance(var_feats[0], Variation)
        assert var_feats[0].start == 10
        assert var_feats[0].end == 10
        assert var_feats[0].strand == Strand.FORWARD
        assert var_feats[0].stable_id == "rs1;id4"
        assert var_feats[0].allele_string == "A/C"
        assert var_feats[0].slice is not None
        assert var_feats[0].slice.seq_region_name == "chr1"
        assert var_feats[0].slice.seq_region_length == 32

        assert var_feats[1].end == 22
        assert var_feats[1].allele_string == "GAT/G"
