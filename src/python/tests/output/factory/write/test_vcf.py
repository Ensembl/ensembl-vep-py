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
"""Tests for VCF output writing."""

from __future__ import annotations

import pytest

from ensembl.vep.output.factory import AnnotatedVariant
from ensembl.vep.output.factory.write import VCFOutputWriter


class TestVCFOutputWriter:
    """Annotated-variant VCF serialization."""

    def test_line_replaces_the_existing_csq_annotation(self) -> None:
        """Write snapshot annotations through the existing VCF CSQ encoder."""
        variant = AnnotatedVariant(
            id="rs1",
            seq_region_name="1",
            start=10,
            end=10,
            strand=1,
            allele_string="A/G",
            input_columns=("1", "10", "rs1", "A", "G", ".", ".", "CSQ=old;DP=3"),
            annotations=(
                {
                    "Allele": "G",
                    "Consequence": ["missense_variant"],
                    "Feature_type": "Transcript",
                },
            ),
        )

        line = VCFOutputWriter().line(variant)

        assert line == ("1\t10\trs1\tA\tG\t.\t.\tDP=3;" "CSQ=G|missense_variant||||Transcript||||||||||||")

    def test_line_requires_original_columns_for_valid_indel_spelling(self) -> None:
        """Reject a snapshot that lacks the anchor base needed to construct valid VCF."""
        variant = AnnotatedVariant(
            id="insertion",
            seq_region_name="1",
            start=11,
            end=10,
            strand=1,
            allele_string="-/AT",
        )

        with pytest.raises(ValueError, match="requires the original input columns"):
            VCFOutputWriter().line(variant)
