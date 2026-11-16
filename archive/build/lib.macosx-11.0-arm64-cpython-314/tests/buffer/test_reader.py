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
"""Tests for VCF input buffering."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from ensembl.vep.buffer.reader import (
    count_vcf_records,
    input_header_lines,
    iter_vcf_buffers,
    iter_vcf_columns,
)

VCF_WITH_SAMPLES = """\
##fileformat=VCFv4.3
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\tSAMPLE_B
chr1\t10\trs1\tA\tC\t10.00\tPASS\tDP=007\tGT\t0|1\t1/1
chr2\t20\trs2\tG\tT\t.\tPASS\t.\tGT\t0/0\t0/1
chr1\t30\trs3\tC\tA\t.\tPASS\t.\tGT\t1/1\t./.
"""


class TestVCFBufferReader:
    """Raw VCF row transport."""

    def test_preserves_format_and_sample_columns(self, tmp_path: Path) -> None:
        """Carry every field without parser normalisation."""
        path = tmp_path / "samples.vcf"
        path.write_text(VCF_WITH_SAMPLES, encoding="utf-8")

        columns = list(iter_vcf_columns(path))

        assert columns[0] == (
            "chr1",
            "10",
            "rs1",
            "A",
            "C",
            "10.00",
            "PASS",
            "DP=007",
            "GT",
            "0|1",
            "1/1",
        )

    def test_preserves_recurrent_region_order_across_buffers(self, tmp_path: Path) -> None:
        """Split region order A, B, A without sorting or combining it."""
        path = tmp_path / "unsorted.vcf"
        path.write_text(VCF_WITH_SAMPLES, encoding="utf-8")

        payloads = list(iter_vcf_buffers(path, batch_size=10))

        assert [payload.index for payload in payloads] == [0, 1, 2]
        assert [payload.sequence_region for payload in payloads] == ["chr1", "chr2", "chr1"]
        assert [payload.columns[0][2] for payload in payloads] == ["rs1", "rs2", "rs3"]

    def test_reads_compressed_headers_and_columns(self, tmp_path: Path) -> None:
        """Treat gzip compression as transport rather than record content."""
        path = tmp_path / "samples.vcf.gz"
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            handle.write(VCF_WITH_SAMPLES)

        headers = input_header_lines(path)
        columns = list(iter_vcf_columns(path))

        assert headers[-1].endswith("FORMAT\tSAMPLE_A\tSAMPLE_B")
        assert columns[-1][8:] == ("GT", "1/1", "./.")
        assert count_vcf_records(path) == 3

    def test_counts_uncompressed_variant_records(self, tmp_path: Path) -> None:
        """Count data rows without treating headers or blank lines as variants."""
        path = tmp_path / "variants.vcf"
        path.write_text(VCF_WITH_SAMPLES + "\n", encoding="utf-8")

        assert count_vcf_records(path) == 3

    def test_rejects_a_non_positive_batch_size(self, tmp_path: Path) -> None:
        """Require each emitted buffer to have capacity for a record."""
        path = tmp_path / "samples.vcf"
        path.write_text(VCF_WITH_SAMPLES, encoding="utf-8")

        with pytest.raises(ValueError, match="at least one"):
            list(iter_vcf_buffers(path, batch_size=0))

    def test_rejects_a_short_vcf_row(self, tmp_path: Path) -> None:
        """Reject data that cannot be reconstructed as a VCF record."""
        path = tmp_path / "short.vcf"
        path.write_text("#CHROM\tPOS\nchr1\t10\n", encoding="utf-8")

        with pytest.raises(ValueError, match="at least eight"):
            list(iter_vcf_columns(path))
