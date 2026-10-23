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
"""Tests for table rendering."""

from __future__ import annotations

import pytest

from ensembl.vep.render.table import EMPTY_CELL, TableFormat, render_table

HEADERS = ("Fact", "Cost tier")
ROWS = (("is_snv", "FREE"), ("overlap_cds", "MAPPING"))


class TestMarkdown:
    """The default format, for pasting into a pull request or Confluence's editor."""

    def test_it_emits_a_header_and_a_rule(self) -> None:
        """A pipe table needs both to render."""
        lines = render_table(HEADERS, ROWS, TableFormat.MARKDOWN).splitlines()
        assert lines[0] == "| Fact | Cost tier |"
        assert lines[1] == "|---|---|"

    def test_every_row_appears(self) -> None:
        """One line per row, after the header and the rule."""
        assert len(render_table(HEADERS, ROWS, TableFormat.MARKDOWN).splitlines()) == 4


class TestConfluence:
    """Classic wiki markup, for pasting into the source view."""

    def test_headers_use_double_pipes(self) -> None:
        """What distinguishes a heading row from a body row."""
        first = render_table(HEADERS, ROWS, TableFormat.CONFLUENCE).splitlines()[0]
        assert first == "||Fact||Cost tier||"

    def test_body_rows_use_single_pipes(self) -> None:
        """The counterpart to the heading row."""
        second = render_table(HEADERS, ROWS, TableFormat.CONFLUENCE).splitlines()[1]
        assert second == "|is_snv|FREE|"


class TestText:
    """Fixed-width output, as the retired scripts printed."""

    def test_columns_are_padded_to_their_widest_value(self) -> None:
        """Otherwise a long value pushes the following column out of alignment."""
        lines = render_table(HEADERS, ROWS, TableFormat.TEXT).splitlines()
        assert lines[0].startswith("Fact       ")

    def test_a_rule_separates_the_header(self) -> None:
        """Matching the section formatter the validator used."""
        assert set(render_table(HEADERS, ROWS, TableFormat.TEXT).splitlines()[1]) <= {"-", " "}

    def test_lines_carry_no_trailing_space(self) -> None:
        """Padding the last column would leave whitespace at the end of every line."""
        for line in render_table(HEADERS, ROWS, TableFormat.TEXT).splitlines():
            assert line == line.rstrip()


class TestHtml:
    """The format a self-contained page embeds."""

    def test_it_emits_a_head_and_a_body(self) -> None:
        """A bare sequence of rows would not be stylable."""
        html = render_table(HEADERS, ROWS, TableFormat.HTML)
        assert "<thead>" in html
        assert "<tbody>" in html

    def test_it_carries_the_class_the_stylesheet_targets(self) -> None:
        """Styling follows ensembl-client's table component rather than an inline style."""
        assert 'class="vep-table"' in render_table(HEADERS, ROWS, TableFormat.HTML)

    def test_cell_text_is_escaped(self) -> None:
        """A term or expression may contain characters that would otherwise close a tag."""
        html = render_table(("Rule",), (("a < b & c",),), TableFormat.HTML)
        assert "a &lt; b &amp; c" in html


class TestSharedBehaviour:
    """What every format does the same way."""

    @pytest.mark.parametrize("table_format", list(TableFormat), ids=lambda f: f.value)
    def test_an_absent_value_is_shown_as_a_placeholder(self, table_format: TableFormat) -> None:
        """An empty cell is indistinguishable from a rendering fault."""
        assert EMPTY_CELL in render_table(("Deps",), ((None,),), table_format)

    @pytest.mark.parametrize("table_format", list(TableFormat), ids=lambda f: f.value)
    def test_a_table_with_no_rows_still_renders_its_headers(self, table_format: TableFormat) -> None:
        """An empty result is a real answer, not an error."""
        assert "Fact" in render_table(HEADERS, (), table_format)

    def test_a_short_row_is_rejected(self) -> None:
        """Silently padding would misalign every column after the gap."""
        with pytest.raises(ValueError, match="Row 0"):
            render_table(HEADERS, (("only",),))

    def test_a_long_row_is_rejected(self) -> None:
        """The counterpart to a short row."""
        with pytest.raises(ValueError, match="Row 0"):
            render_table(HEADERS, (("a", "b", "c"),))
