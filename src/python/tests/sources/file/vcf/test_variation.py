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
"""Unit tests for the VcfVariationSource class."""

import pytest

from ensembl.vep.ensembl_vep_core import (  # noqa: F401,F403 pylint: disable=E0401,E0611
    VCFRecord
)
from ensembl.vep.sources.file.vcf.variation import VcfVariationSource

# @pytest.fixture
# def snp_vcf_record() -> VCFRecord:
#     """Create a VcfRecord object for a SNP for testing."""
#     return VCFRecord(
#         chrom="chr1",
#         pos=10,
#         id="rs1",
#         reference="A",
#         alt="C",
#         qual=40,
#         filter="PASS",
#     )

class TestVcfVariationSourceInitialization:
    """Test VcfVariationSource initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test VcfVariationSource initialization without keywords raises error."""

        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = VcfVariationSource("sample.vcf")  # type: ignore[misc]

    def test_initialization_sets_core_attributes(self) -> None:
        """Test VcfVariationSource stores source and parser attributes."""

        source = VcfVariationSource()

        assert source.file_path is None
        assert source.adaptor.format == "vcf"  # pylint: disable=protected-access
        assert source.adaptor.file_path is None  # pylint: disable=protected-access
        assert source.adaptor.tabix_path == "tabix"  # pylint: disable=protected-access
        assert source.adaptor.reader_mode == "auto"  # pylint: disable=protected-access


class TestVcfVariationProperties:
    """Tests property accessors and default values."""

    def test_initialization(self, snp_vcf_record: VCFRecord) -> None:
        """VcfVariationSource"""
        vcf_var_source = VcfVariationSource()


