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
"""Rendering a table in the formats this project publishes to."""

from __future__ import annotations

from collections.abc import Sequence
from enum import Enum
from xml.sax.saxutils import escape

from ensembl.vep.render.theme import Colour, Metric

__all__ = ["TableFormat", "render_table"]

EMPTY_CELL = "--"
"""What an absent value is shown as, matching the retired scripts' convention."""


class TableFormat(Enum):
    """The formats a table can be rendered in."""

    MARKDOWN = "markdown"
    CONFLUENCE = "confluence"
    TEXT = "text"
    HTML = "html"


def _cells(row: Sequence[object]) -> list[str]:
    """Return one row's values as display strings."""
    rendered = []
    for value in row:
        text = "" if value is None else str(value)
        rendered.append(text or EMPTY_CELL)
    return rendered


def _column_widths(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> list[int]:
    """Return the width each column needs to fit its widest value."""
    widths = [len(header) for header in headers]
    for row in rows:
        for index, text in enumerate(_cells(row)):
            widths[index] = max(widths[index], len(text))
    return widths


def _render_markdown(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    """Return a GitHub-flavoured pipe table."""
    lines = ["| " + " | ".join(headers) + " |"]
    lines.append("|" + "|".join("---" for _ in headers) + "|")
    for row in rows:
        lines.append("| " + " | ".join(_cells(row)) + " |")
    return "\n".join(lines)


def _render_confluence(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    """Return classic Confluence wiki markup, for pasting into the source view."""
    lines = ["||" + "||".join(headers) + "||"]
    for row in rows:
        lines.append("|" + "|".join(_cells(row)) + "|")
    return "\n".join(lines)


def _render_text(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    """Return a fixed-width table, as the retired scripts printed."""
    widths = _column_widths(headers, rows)
    header_line = "  ".join(header.ljust(widths[index]) for index, header in enumerate(headers))
    lines = [header_line.rstrip(), "  ".join("-" * width for width in widths)]
    for row in rows:
        cells = _cells(row)
        line = "  ".join(cell.ljust(widths[index]) for index, cell in enumerate(cells))
        lines.append(line.rstrip())
    return "\n".join(lines)


def _render_html(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    """Return an HTML table styled as ``ensembl-client``'s table component.

    Follows the site rather than the old VEP summary page, which striped alternate rows.
    """
    head = "".join(f"<th>{escape(header)}</th>" for header in headers)
    body_rows = []
    for row in rows:
        cells = "".join(f"<td>{escape(cell)}</td>" for cell in _cells(row))
        body_rows.append(f"<tr>{cells}</tr>")
    body = "\n    ".join(body_rows)
    return (
        f'<table class="vep-table">\n  <thead><tr>{head}</tr></thead>\n'
        f"  <tbody>\n    {body}\n  </tbody>\n</table>"
    )


TABLE_STYLES = f"""
.vep-table {{
  border-collapse: collapse;
  background-color: {Colour.WHITE.var};
}}

.vep-table th {{
  font-size: 11px;
  font-weight: {Metric.FONT_WEIGHT_LIGHT.value};
  text-align: left;
  padding: 0 32px 4px 24px;
  box-shadow: inset 0 -1.5px 0 {Colour.BLUE.var};
}}

.vep-table td {{
  padding: 9px 32px 9px 24px;
  vertical-align: middle;
}}

.vep-table tr + tr {{
  border-top: 1px solid {Colour.GREY.var};
}}
"""
"""CSS the HTML format depends on, for a page that embeds a table."""

_RENDERERS = {
    TableFormat.MARKDOWN: _render_markdown,
    TableFormat.CONFLUENCE: _render_confluence,
    TableFormat.TEXT: _render_text,
    TableFormat.HTML: _render_html,
}


def render_table(
    headers: Sequence[str],
    rows: Sequence[Sequence[object]],
    table_format: TableFormat = TableFormat.MARKDOWN,
) -> str:
    """Return a table in one of the published formats.

    Args:
        headers (Sequence[str]): The column headings.
        rows (Sequence[Sequence[object]]): One sequence of values per row, each the same length as
            ``headers``.
        table_format (TableFormat): Which format to render.

    Returns:
        str: The rendered table.

    Raises:
        ValueError: If a row has a different number of values than there are headers.
    """
    for index, row in enumerate(rows):
        if len(row) != len(headers):
            raise ValueError(f"Row {index} has {len(row)} values, expected {len(headers)}")
    return _RENDERERS[table_format](headers, rows)
