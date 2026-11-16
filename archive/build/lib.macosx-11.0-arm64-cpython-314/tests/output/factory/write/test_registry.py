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
"""Tests for output writer selection."""

from __future__ import annotations

import pytest

from ensembl.vep.output.factory.write import (
    JSONOutputWriter,
    VCFOutputWriter,
    available_output_writers,
    get_output_writer,
    register_output_writer,
)
from ensembl.vep.output.factory.write import registry as registry_module


class TestOutputWriterRegistry:
    """Output-format writer selection."""

    @pytest.mark.parametrize(
        ("output_format", "writer_type"),
        (("json", JSONOutputWriter), ("vcf", VCFOutputWriter)),
    )
    def test_get_output_writer_returns_the_selected_writer(
        self,
        output_format: str,
        writer_type: type[JSONOutputWriter] | type[VCFOutputWriter],
    ) -> None:
        """Select the writer registered for an output format."""
        assert isinstance(get_output_writer(output_format), writer_type)

    def test_get_output_writer_rejects_an_unimplemented_format(self) -> None:
        """Explain which formats have a writer implementation."""
        with pytest.raises(ValueError, match="no writer for 'tab'; supported: json, vcf"):
            get_output_writer("tab")

    def test_register_output_writer_adds_an_extensible_format(self) -> None:
        """Register and construct a writer without changing the central lookup."""
        command_lines: list[str] = []

        @register_output_writer("test.json")
        def build_writer(command_line: str) -> JSONOutputWriter:
            command_lines.append(command_line)
            return JSONOutputWriter()

        try:
            writer = get_output_writer("test.json", command_line="vep run")

            assert isinstance(writer, JSONOutputWriter)
            assert command_lines == ["vep run"]
            assert "test.json" in available_output_writers()
        finally:
            registry_module._OUTPUT_WRITER_REGISTRY.pop("test.json", None)  # pylint: disable=protected-access
