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
"""A report: titled sections of tables and charts, rendered to one of three formats.

The one place a run summary, a fact profile, or any later report becomes HTML, text or JSON. A
producer describes *what* it has -- headings, rows, labelled quantities -- and never how it should
look, which is what lets it live outside this package without a renderer of its own.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from enum import Enum

from attrs import field, frozen

from ensembl.vep.render.chart import Slice, render_bars, render_legend, render_pie
from ensembl.vep.render.page import Section, render_page
from ensembl.vep.render.table import TableFormat, render_table

__all__ = ["ChartKind", "ReportFormat", "ReportSection", "render_report"]


class ReportFormat(Enum):
    """The formats a report can be rendered in."""

    HTML = "html"
    TEXT = "text"
    JSON = "json"


class ChartKind(Enum):
    """How a section's quantities are drawn, if at all."""

    NONE = "none"
    PIE = "pie"
    BARS = "bars"


@frozen(kw_only=True, auto_attribs=True)
class ReportSection:
    """One block of a report: a table, a chart, or both.

    Attributes:
        identifier (str): The anchor, and what a page's side menu links to.
        title (str): The heading.
        headers (tuple[str, ...]): Column headings of the table, if any.
        rows (tuple[tuple[object, ...], ...]): The table's rows.
        slices (tuple[Slice, ...]): Labelled quantities, drawn as ``chart`` and, in the text
            formats, listed as a two-column table.
        chart (ChartKind): How the slices are drawn in HTML.
    """

    identifier: str
    title: str
    headers: tuple[str, ...] = field(default=())
    rows: tuple[tuple[object, ...], ...] = field(default=())
    slices: tuple[Slice, ...] = field(default=())
    chart: ChartKind = field(default=ChartKind.NONE)


SLICE_HEADERS = ("Label", "Count")


def _html_blocks(section: ReportSection) -> tuple[str, ...]:
    """Return a section's table and chart as HTML fragments."""
    blocks = []
    if section.rows or section.headers:
        blocks.append(render_table(section.headers, section.rows, TableFormat.HTML))
    if section.slices and section.chart is ChartKind.PIE:
        blocks.append(render_pie(section.slices, section.title))
        blocks.append(render_legend(section.slices))
    elif section.slices:
        blocks.append(render_bars(section.slices, section.title))
    return tuple(blocks)


def _text_blocks(section: ReportSection) -> list[str]:
    """Return a section's heading and its tables as text."""
    blocks = [f"\n## {section.title}\n"]
    if section.rows or section.headers:
        blocks.append(render_table(section.headers, section.rows, TableFormat.TEXT))
    if section.slices:
        rows = [(item.label, item.value) for item in section.slices]
        blocks.append(render_table(SLICE_HEADERS, rows, TableFormat.TEXT))
    return blocks


def _json_section(section: ReportSection) -> dict[str, object]:
    """Return a section as plain JSON-ready data."""
    return {
        "identifier": section.identifier,
        "title": section.title,
        "headers": list(section.headers),
        "rows": [list(row) for row in section.rows],
        "slices": [{"label": item.label, "value": item.value} for item in section.slices],
    }


def render_report(
    title: str, subtitle: str, sections: Sequence[ReportSection], report_format: str | ReportFormat
) -> str:
    """Render a report in one of the published formats.

    Args:
        title (str): The document title.
        subtitle (str): One line under it.
        sections (Sequence[ReportSection]): The blocks, in order.
        report_format (str | ReportFormat): Which format to render.

    Returns:
        str: The rendered report.
    """
    chosen = ReportFormat(report_format)
    if chosen is ReportFormat.HTML:
        page_sections = []
        for section in sections:
            page_sections.append(
                Section(identifier=section.identifier, title=section.title, blocks=_html_blocks(section))
            )
        return render_page(title, page_sections, subtitle)
    if chosen is ReportFormat.JSON:
        payload = {
            "title": title,
            "subtitle": subtitle,
            "sections": [_json_section(section) for section in sections],
        }
        return json.dumps(payload, indent=2, default=str)

    blocks = [f"# {title}", subtitle]
    for section in sections:
        blocks.extend(_text_blocks(section))
    return "\n".join(blocks)
