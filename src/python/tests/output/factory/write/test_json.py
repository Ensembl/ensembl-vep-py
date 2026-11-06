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
"""Tests for JSON output writing."""

from __future__ import annotations

from ensembl.vep.output.factory import AnnotatedVariant
from ensembl.vep.output.factory.write.json import JSONOutputWriter


class TestJSONOutputWriter:
    """JSON-lines output."""

    def test_headers_returns_no_lines(self) -> None:
        """Write no headers for a JSON-lines stream."""
        assert not JSONOutputWriter().headers(("##fileformat=VCFv4.1",))

    def test_line_writes_one_compact_json_record(self) -> None:
        """Write one compact Perl VEP-style JSON object per variant."""
        variant = AnnotatedVariant(
            id="rs1",
            seq_region_name="1",
            start=10,
            end=10,
            strand=1,
            allele_string="A/G",
            input_columns=("1", "10", "rs1", "A", "G"),
        )

        assert JSONOutputWriter().line(variant) == (
            '{"id":"rs1","seq_region_name":"1","start":10,"end":10,"strand":1,'
            '"allele_string":"A/G","assembly_name":null,"most_severe_consequence":"?",'
            '"input":"1\\t10\\trs1\\tA\\tG"}'
        )
