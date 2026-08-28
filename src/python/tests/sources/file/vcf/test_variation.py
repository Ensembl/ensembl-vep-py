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
from ensembl.vep.common import Locus
from ensembl.vep.common.enums import Strand
from ensembl.vep.models.variation.variation import Variation
from ensembl.vep.models.variation.structural_variation import StructuralVariation
from ensembl.vep.parsers.vcf import VCFParser
from ensembl.vep.sources.file.vcf.variation import VCFVariationSource
from tests.models.conftest import DummySequenceSource

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
CHR2_VAR_ROW1 = "chr2\t100\trs3\tC\tCAT\t.\tPASS\tNS=2\tGT\t1/1\t0/1"
SAMPLE = "\n".join([*HEADER_LINES, CHR1_VAR_ROW1, CHR1_VAR_ROW2, CHR2_VAR_ROW1]) + "\n"

CHR1_SV_DEL_SVTYPE = "chr1\t50\tsv1\tA\t<DEL>\t.\tPASS\tEND=80;SVTYPE=DEL;\tGT\t0/1\t0/0"
CHR1_SV_DEL_NEGATIVE_SVLEN = "chr1\t100\tsv2\tA\t<DEL>\t.\tPASS\tSVLEN=-50;SVTYPE=DEL\tGT\t0/1\t0/0"
CHR1_SV_INS_IMPRECISE = (
    "chr1\t14\t.\tT\t<INS>\t.\t.\tIMPRECISE;SVLEN=20;CILEN=-10,10;CIPOS=-10,10\tGT\t0/1\t1/1"
)
CHR1_SV_DUP_SVLEN = "chr1\t101\tsv3\tT\t<DUP>\t.\tPASS\tSVLEN=20;CIPOS=-5,2;CIEND=-3,4\tGT\t0/1\t0/0"
CHR1_SV_INV = "chr1\t201\tsv4\tT\t<INV>,<INV>\t.\tPASS\tSVLEN=30,.;CIPOS=-5,2;CIEND=-3,4\tGT\t0/1\t0/0"
CHR1_SV_CNV = (
    "chr1\t301\tsv5\tT\t<CNV>,<CNV>\t.\tPASS\tSVLEN=40,50;CIPOS=-5,2,.,.;CIEND=-3,4,.,.\tGT\t0/1\t0/2"
)
CHR1_SV_BND1 = "chrA\t2\t.\tT\t<DEL>\t.\t.\tSVLEN=.;SVCLAIM=D;EVENT=DEL_split_bp_cn\tGT\t0/1\t1/1"
CHR1_SV_BND2 = (
    "chrA\t2\tdelbp1\tT\tT[chrA:5[\t.\t.\tSVLEN=.;MATEID=delbp2;EVENT=DEL_split_bp_cn\tGT\t0/1\t1/1"
)
CHR1_SV_BND3 = "chrA\t2\tdelbp2\tA\t]chrA:2]A\t.\t.\tMATEID=delbp1;EVENT=DEL_split_bp_cn\tGT\t0/1\t1/1"
CHR1_SV_BND4 = "chrA\t14\t.\tG\t.CCCCCCG\t.\t.\tEVENT=single_breakend\tGT\t0/1\t1/1"

SV_SAMPLE = (
    "\n".join(
        [
            *HEADER_LINES,
            CHR1_SV_DEL_SVTYPE,
            CHR1_SV_INS_IMPRECISE,
            CHR1_SV_DUP_SVLEN,
            CHR1_SV_INV,
            CHR1_SV_CNV,
            CHR1_SV_BND1,
            CHR1_SV_BND2,
            CHR1_SV_BND3,
            CHR1_SV_BND4,
        ]
    )
    + "\n"
)


def create_vcf_records_for_testing(text: str) -> list[VCFRecord]:
    """Create a list of VCFRecord objects for testing."""

    return VCFParser().parse_text(text)


def create_sequence_source_for_testing() -> DummySequenceSource:
    """Create a sequence source covering every region used in the VCF fixtures."""

    return DummySequenceSource(
        data={
            "chr1": "N" * 400,
            "chr2": "N" * 200,
            "chrA": "N" * 100,
        }
    )


@pytest.fixture(name="vcf_records")
def vcf_records_fixture(text: str = SAMPLE) -> list[VCFRecord]:
    """Create a list of SNV VCFRecord objects for testing."""

    return create_vcf_records_for_testing(text)


@pytest.fixture(name="sv_vcf_records")
def sv_vcf_records_fixture(text: str = SV_SAMPLE) -> list[VCFRecord]:
    """Create a list VCFRecord objects for a SNP for testing."""

    return create_vcf_records_for_testing(text)


def assert_first_variation_object(var: Variation) -> None:
    """Test the variation object of CHR1_VAR_ROW1"""
    assert isinstance(var, Variation)
    assert var.start == 10
    assert var.end == 10
    assert var.strand == Strand.FORWARD
    assert var.stable_id == "rs1;id4"
    assert var.allele_string == "A/C"


def assert_first_structural_variation_object(var: StructuralVariation) -> None:
    """Test the variation object of CHR1_VAR_ROW1"""
    assert var.start == 50
    assert var.end == 80
    assert var.allele_string == "A/<DEL>"


class TestVCFVariationSourceInitialization:
    """Test VCFVariationSource initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test VCFVariationSource initialization without keywords raises error."""

        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = VCFVariationSource(DummySequenceSource(), file_path="sample.vcf")  # type: ignore[misc]

    def test_initialization_sets_core_attributes(self) -> None:
        """Test VCFVariationSource stores source and parser attributes."""

        source = VCFVariationSource(sequence_source=DummySequenceSource(), file_path="sample.vcf")

        assert source.file_path == "sample.vcf"
        assert source.adaptor.format == "vcf"  # pylint: disable=protected-access
        assert source.adaptor.file_path is None  # pylint: disable=protected-access
        assert source.adaptor.tabix_path == "tabix"  # pylint: disable=protected-access
        assert source.adaptor.reader_mode == "auto"  # pylint: disable=protected-access


class TestVCFVariationSourceMethods:
    """Test the VCFVariationSource methods."""

    def test_stable_id_from_vcf_record(self, vcf_records: list[VCFRecord]) -> None:
        """Test _stable_id_from_vcf_record method."""
        assert (
            VCFVariationSource._stable_id_from_vcf_record(vcf_records[0])  # pylint: disable=protected-access
            == "rs1;id4"
        )
        assert (
            VCFVariationSource._stable_id_from_vcf_record(vcf_records[1])  # pylint: disable=protected-access
            is None
        )
        assert (
            VCFVariationSource._stable_id_from_vcf_record(vcf_records[2])  # pylint: disable=protected-access
            == "rs3"
        )

    def test_end_and_svlen_for_sva_from_vcf_record(
        self, sv_vcf_records: list[VCFRecord], caplog: pytest.LogCaptureFixture
    ) -> None:
        """Test _end_and_svlen_for_sva_from_vcf_record method."""
        source = VCFVariationSource(sequence_source=DummySequenceSource(), file_path="sample.vcf")

        # test for incorrect allele_index
        with caplog.at_level(logging.WARNING):
            results = source._end_and_svlen_for_sva_from_vcf_record(  # pylint: disable=protected-access
                sv_vcf_records[0], allele_index=1
            )
            assert (
                "Allele index, 1, is too large, sv1" in caplog.text
            ), "Incorrect warning in _end_and_svlen_for_sva_from_vcf_record"
        assert results == (None, None)

        assert source._end_and_svlen_for_sva_from_vcf_record(  # pylint: disable=protected-access
            sv_vcf_records[0], allele_index=0
        ) == (
            80,
            30,
        )
        assert source._end_and_svlen_for_sva_from_vcf_record(  # pylint: disable=protected-access
            sv_vcf_records[1], allele_index=0
        ) == (
            14,
            20,
        )
        assert source._end_and_svlen_for_sva_from_vcf_record(  # pylint: disable=protected-access
            sv_vcf_records[2], allele_index=0
        ) == (
            121,
            20,
        )
        assert source._end_and_svlen_for_sva_from_vcf_record(  # pylint: disable=protected-access
            sv_vcf_records[3], allele_index=0
        ) == (
            231,
            30,
        )
        assert source._end_and_svlen_for_sva_from_vcf_record(  # pylint: disable=protected-access
            sv_vcf_records[3], allele_index=1
        ) == (
            201,
            None,
        )
        assert source._end_and_svlen_for_sva_from_vcf_record(  # pylint: disable=protected-access
            sv_vcf_records[4], allele_index=0
        ) == (
            341,
            40,
        )
        assert source._end_and_svlen_for_sva_from_vcf_record(  # pylint: disable=protected-access
            sv_vcf_records[4], allele_index=1
        ) == (
            351,
            50,
        )
        assert source._end_and_svlen_for_sva_from_vcf_record(  # pylint: disable=protected-access
            sv_vcf_records[6], allele_index=0
        ) == (
            None,
            None,
        )

    def test_negative_svlen_is_treated_as_positive(self) -> None:
        """Test negative SVLEN values produce positive lengths and forward coordinates."""
        record = create_vcf_records_for_testing(
            "\n".join([*HEADER_LINES, CHR1_SV_DEL_NEGATIVE_SVLEN]) + "\n"
        )[0]
        source = VCFVariationSource(sequence_source=DummySequenceSource(), file_path="sample.vcf")

        assert source._end_and_svlen_for_sva_from_vcf_record(  # pylint: disable=protected-access
            record, allele_index=0
        ) == (150, 50)

        allele = source._create_sv_allele_from_record(  # pylint: disable=protected-access
            record, allele_index=0, seq_region_slice=None
        )
        assert allele.start == 101
        assert allele.end == 150
        assert allele.svlen == 50
        assert allele.svtype == "DEL"

    def test_int_info_value(self, sv_vcf_records: list[VCFRecord]) -> None:
        """Test _int_info_value method."""

        source = VCFVariationSource(sequence_source=DummySequenceSource(), file_path="sample.vcf")

        assert source._int_info_value(sv_vcf_records[0], "END") == 80  # pylint: disable=protected-access
        assert source._int_info_value(sv_vcf_records[5], "END") is None  # pylint: disable=protected-access
        assert source._int_info_value(sv_vcf_records[3], "SVLEN") is None  # pylint: disable=protected-access

    def test_is_bases_only(self, vcf_records: list[VCFRecord], sv_vcf_records: list[VCFRecord]) -> None:
        """Test _is_bases_only method."""
        source = VCFVariationSource(sequence_source=DummySequenceSource(), file_path="sample.vcf")

        for record in vcf_records:
            assert source._is_bases_only(record)  # pylint: disable=protected-access

        for record in sv_vcf_records:
            assert not source._is_bases_only(record)  # pylint: disable=protected-access

    def test_has_svtype(self, vcf_records: list[VCFRecord], sv_vcf_records: list[VCFRecord]) -> None:
        """Test _has_svtype method."""
        source = VCFVariationSource(sequence_source=DummySequenceSource(), file_path="sample.vcf")

        for record in vcf_records:
            assert not source._has_svtype(record)  # pylint: disable=protected-access

        assert source._has_svtype(sv_vcf_records[0])  # pylint: disable=protected-access
        assert not source._has_svtype(sv_vcf_records[1])  # pylint: disable=protected-access

    def test_is_symbolic(self, vcf_records: list[VCFRecord], sv_vcf_records: list[VCFRecord]) -> None:
        """Test _is_symbolic method."""
        source = VCFVariationSource(sequence_source=DummySequenceSource(), file_path="sample.vcf")

        for record in vcf_records:
            assert not source._is_symbolic(record)  # pylint: disable=protected-access

        for record in sv_vcf_records[0:6]:
            assert source._is_symbolic(record)  # pylint: disable=protected-access
        for record in sv_vcf_records[6:]:
            assert not source._is_symbolic(record)  # pylint: disable=protected-access

    def test_is_breakend(self, vcf_records: list[VCFRecord]) -> None:
        """Test _is_breakend method."""
        SV_LINES = (
            "\n".join(
                [*HEADER_LINES, CHR1_SV_DEL_SVTYPE, CHR1_SV_BND1, CHR1_SV_BND2, CHR1_SV_BND3, CHR1_SV_BND4]
            )
            + "\n"
        )
        sv_vcf_records = create_vcf_records_for_testing(SV_LINES)

        source = VCFVariationSource(sequence_source=DummySequenceSource(), file_path="sample.vcf")

        for record in vcf_records:
            assert not source._is_breakend(record)  # pylint: disable=protected-access

        # NOTE: NEED TO THINK ABOUT BREAKEND RECORDS COMPOSED OF MULTIPLE LINES
        assert not source._is_breakend(sv_vcf_records[0])  # pylint: disable=protected-access
        assert not source._is_breakend(sv_vcf_records[1])  # pylint: disable=protected-access
        assert source._is_breakend(sv_vcf_records[2])  # pylint: disable=protected-access
        assert source._is_breakend(sv_vcf_records[3])  # pylint: disable=protected-access
        assert source._is_breakend(sv_vcf_records[4])  # pylint: disable=protected-access

    def test_is_symbolic_insertion(self) -> None:
        """Test _is_symbolic_insertion method."""
        source = VCFVariationSource(sequence_source=DummySequenceSource(), file_path="sample.vcf")

        assert source._is_symbolic_insertion("INS")  # pylint: disable=protected-access
        assert source._is_symbolic_insertion("<INS>")  # pylint: disable=protected-access
        assert not source._is_symbolic_insertion("<DEL>")  # pylint: disable=protected-access

    def test_is_structural_variation_record(
        self, vcf_records: list[VCFRecord], sv_vcf_records: list[VCFRecord]
    ) -> None:
        """Test _is_structural_variation_record method."""
        source = VCFVariationSource(sequence_source=DummySequenceSource(), file_path="sample.vcf")

        for record in vcf_records:
            assert not source._is_structural_variation_record(record)  # pylint: disable=protected-access

        for record in sv_vcf_records:
            assert source._is_structural_variation_record(record)  # pylint: disable=protected-access

    def test_remove_anchor_bases(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test _remove_anchor_bases method."""
        source = VCFVariationSource(sequence_source=DummySequenceSource(), file_path="sample.vcf")

        assert source._remove_anchor_bases(  # pylint: disable=protected-access
            "CAT", "C", start=10, end=12
        ) == ("AT", "-", 11, 12, True)
        assert source._remove_anchor_bases(  # pylint: disable=protected-access
            "G", "GAT", start=10, end=10
        ) == ("-", "AT", 11, 10, True)
        assert source._remove_anchor_bases(  # pylint: disable=protected-access
            "ACGT", "ACCT", start=10, end=13
        ) == ("GT", "CT", 12, 13, True)
        assert source._remove_anchor_bases(  # pylint: disable=protected-access
            "A", "T", start=10, end=10
        ) == ("A", "T", 10, 10, False)
        assert source._remove_anchor_bases(  # pylint: disable=protected-access
            "AT", ["A", "ATT"], start=10, end=11
        ) == ("T", ["-", "TT"], 11, 11, True)

        with caplog.at_level(logging.WARNING):
            assert source._remove_anchor_bases(  # pylint: disable=protected-access
                "A", "AT", start=10, end=10, strand=Strand.UNKNOWN
            ) == ("A", "AT", 10, 10, False)
        assert "Strand is unknown, cannot remove anchor bases with confidence" in caplog.text

    def test_remove_anchor_bases_symbolic_allele(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test _remove_anchor_bases_symbolic_allele method."""
        source = VCFVariationSource(sequence_source=DummySequenceSource(), file_path="sample.vcf")

        # test for unknown strand warning
        with caplog.at_level(logging.WARNING):
            source._remove_anchor_bases_symbolic_allele(  # pylint: disable=protected-access
                "A", "<DEL>", start=20, end=50, strand=Strand.UNKNOWN
            )
            assert (
                "Strand is unknown, cannot trim alleles with confidence" in caplog.text
            ), "Incorrect strand unknown warning in _remove_anchor_bases_symbolic_allele"

        # test for allele already trimmed warning
        with caplog.at_level(logging.WARNING):
            source._remove_anchor_bases_symbolic_allele(  # pylint: disable=protected-access
                "-", "<DEL>", start=20, end=50, strand=Strand.FORWARD
            )
            assert (
                "Allele is already trimmed" in caplog.text
            ), "Incorrect allele already trimmed warning in _remove_anchor_bases_symbolic_allele"

        assert source._remove_anchor_bases_symbolic_allele(  # pylint: disable=protected-access
            "A", "<DEL>", start=20, end=50
        ) == (
            "-",
            "<DEL>",
            21,
            50,
            True,
        )
        assert source._remove_anchor_bases_symbolic_allele(  # pylint: disable=protected-access
            "A", "<DEL>", start=20, end=50, strand=Strand.REVERSE
        ) == (
            "-",
            "<DEL>",
            20,
            49,
            True,
        )
        assert source._remove_anchor_bases_symbolic_allele(  # pylint: disable=protected-access
            "A", "<INS>", start=14, end=14
        ) == (
            "-",
            "<INS>",
            15,
            14,
            True,
        )


class TestVCFVariationSourceRecordsToFeatures:
    """Tests functions for creating Variation/StructuralVariation objects."""

    def test_records_to_variation(self, vcf_records: list[VCFRecord]) -> None:
        """Test _records_to_variation method."""

        seq_source = create_sequence_source_for_testing()
        vcf_var_source = VCFVariationSource(
            sequence_source=seq_source,
            file_path="sample.vcf",
        )

        var_feats = vcf_var_source._records_to_variation(  # pylint: disable=protected-access
            records=vcf_records
        )

        assert len(var_feats) == 3
        assert_first_variation_object(var_feats[0])
        assert var_feats[0].slice is not None
        assert var_feats[0].slice.seq_region_name == "chr1"
        assert var_feats[0].slice.seq_region_length == 400

        assert var_feats[1].start == 21
        assert var_feats[1].end == 22
        assert var_feats[1].allele_string == "AT/-"
        assert var_feats[1].uploaded_locus == Locus(start=20, end=22, allele_string="GAT/G")

        assert var_feats[2].start == 101
        assert var_feats[2].end == 100
        assert var_feats[2].allele_string == "-/AT"
        assert var_feats[2].uploaded_locus == Locus(start=100, end=100, allele_string="C/CAT")

    def test_create_sv_allele_from_record(self, sv_vcf_records: list[VCFRecord]) -> None:
        """Test _create_sv_allele_from_record method."""
        vcf_var_source = VCFVariationSource(
            sequence_source=create_sequence_source_for_testing(),
            file_path="sample.vcf",
        )

        sva = vcf_var_source._create_sv_allele_from_record(  # pylint: disable=protected-access
            sv_vcf_records[0], allele_index=0, seq_region_slice=None
        )
        assert sva.start == 51
        assert sva.svtype == "DEL"

    def test_records_to_structural_variation(self, sv_vcf_records: list[VCFRecord]) -> None:
        """Test _records_to_structural_variation method."""
        vcf_var_source = VCFVariationSource(
            sequence_source=create_sequence_source_for_testing(),
            file_path="sample.vcf",
        )

        var_feats = vcf_var_source._records_to_structural_variation(  # pylint: disable=protected-access
            records=sv_vcf_records
        )
        assert_first_structural_variation_object(var_feats[0])

    def test_records_to_features(
        self,
        vcf_records: list[VCFRecord],
        sv_vcf_records: list[VCFRecord],
    ) -> None:
        """Test _records_to_features method."""
        vcf_var_source = VCFVariationSource(
            sequence_source=create_sequence_source_for_testing(),
            file_path="sample.vcf",
        )

        var_feats = vcf_var_source._records_to_features(  # pylint: disable=protected-access
            [*vcf_records, *sv_vcf_records]
        )

        assert len(var_feats) == 12
        assert_first_variation_object(var_feats[0])

        assert isinstance(var_feats[3], StructuralVariation)
        assert_first_structural_variation_object(var_feats[3])
