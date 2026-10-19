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
"""Tests for resolving a report's format against the kind of report asked for."""

from __future__ import annotations

import pytest

from ensembl.vep.config.report import GRAPH_FORMATS, REPORT_FORMATS, TABLE_FORMATS, resolve_report_format


class TestVocabulary:
    """What ``--format`` accepts."""

    def test_the_two_kinds_do_not_overlap(self) -> None:
        """A format renders rows or a diagram, never both, which is what makes the check decidable."""
        assert not set(TABLE_FORMATS) & set(GRAPH_FORMATS)

    def test_every_format_belongs_to_one_kind(self) -> None:
        """The option's vocabulary is exactly the two kinds together."""
        assert set(REPORT_FORMATS) == set(TABLE_FORMATS) | set(GRAPH_FORMATS)


class TestDefaults:
    """What is rendered when no format is named."""

    def test_a_table_takes_the_command_default(self) -> None:
        """Each command keeps its own table default."""
        assert resolve_report_format(None, graph=False, table_default="markdown") == "markdown"

    def test_a_diagram_takes_the_diagram_default(self) -> None:
        """--graph on its own draws, rather than demanding a format as well."""
        assert resolve_report_format(None, graph=True, table_default="markdown") == "svg"

    def test_the_diagram_default_can_be_chosen_per_command(self) -> None:
        """A command whose diagrams are better read as source can say so."""
        assert (
            resolve_report_format(None, graph=True, table_default="text", graph_default="mermaid")
            == "mermaid"
        )


class TestCompatibility:
    """Which combinations can be rendered, and what the refusals say."""

    @pytest.mark.parametrize("format_name", TABLE_FORMATS)
    def test_a_table_format_renders_a_table(self, format_name: str) -> None:
        """Any table format is honoured when no diagram was asked for."""
        assert resolve_report_format(format_name, graph=False, table_default="text") == format_name

    @pytest.mark.parametrize("format_name", GRAPH_FORMATS)
    def test_a_diagram_format_renders_a_diagram(self, format_name: str) -> None:
        """Any diagram format is honoured once --graph asks for one."""
        assert resolve_report_format(format_name, graph=True, table_default="text") == format_name

    def test_a_diagram_format_without_the_flag_is_refused(self) -> None:
        """The message names the flag that would make the request renderable."""
        with pytest.raises(ValueError, match="pass --graph"):
            resolve_report_format("svg", graph=False, table_default="text")

    def test_a_table_format_with_the_flag_is_refused(self) -> None:
        """The message lists what a diagram can be drawn as."""
        with pytest.raises(ValueError, match="--graph draws svg, dot, mermaid"):
            resolve_report_format("markdown", graph=True, table_default="text")

    def test_a_report_with_no_diagram_form_refuses_the_flag(self) -> None:
        """A timing table cannot be drawn, so asking says so rather than rendering a table."""
        with pytest.raises(ValueError, match="no diagram form"):
            resolve_report_format(None, graph=True, table_default="text", graph_formats=())

    def test_a_report_with_no_diagram_form_refuses_a_diagram_format(self) -> None:
        """The same request through --format alone is refused for the same reason."""
        with pytest.raises(ValueError, match="pass --graph"):
            resolve_report_format("dot", graph=False, table_default="text", graph_formats=())

    def test_an_unknown_format_lists_the_known_ones(self) -> None:
        """Reached when a command declares choices wider than the vocabulary."""
        with pytest.raises(ValueError, match="unknown format 'pdf'"):
            resolve_report_format("pdf", graph=False, table_default="text")

    def test_a_diagram_format_outside_this_report_is_refused(self) -> None:
        """A command may draw only some of the diagram formats."""
        with pytest.raises(ValueError, match="--graph draws mermaid"):
            resolve_report_format("dot", graph=True, table_default="text", graph_formats=("mermaid",))
