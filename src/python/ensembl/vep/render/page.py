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
"""A self-contained HTML page.

Mirrors the information architecture of VEP's summary page -- a masthead, a fixed side menu of
anchors, then one section per topic -- with the Ensembl skin and no external requests. Everything
the page needs is inlined, so it opens from a filesystem with no network.
"""

from __future__ import annotations

from collections.abc import Sequence
from xml.sax.saxutils import escape

from attrs import define, field

from ensembl.vep.render.table import TABLE_STYLES
from ensembl.vep.render.theme import Colour, Metric, stylesheet

__all__ = ["Section", "render_page"]

LAYOUT_STYLES = f"""
.vep-masthead {{
  border-bottom: 1px solid {Colour.MEDIUM_LIGHT_GREY.var};
  padding: {Metric.STANDARD_GUTTER.value};
}}

.vep-masthead h1 {{
  margin: 0;
  font-size: 20px;
  font-weight: {Metric.FONT_WEIGHT_LIGHT.value};
}}

.vep-masthead p {{
  margin: 4px 0 0;
  color: {Colour.DARK_GREY.var};
}}

.vep-layout {{
  display: flex;
  align-items: flex-start;
  gap: {Metric.STANDARD_GUTTER.value};
  padding: {Metric.STANDARD_GUTTER.value};
}}

.vep-menu {{
  position: sticky;
  top: {Metric.STANDARD_GUTTER.value};
  min-width: 200px;
  max-height: calc(100vh - 2 * {Metric.STANDARD_GUTTER.value});
  overflow-y: auto;
  border-left: 2px solid {Colour.BLUE.var};
  padding-left: 12px;
}}

.vep-menu ul {{
  list-style: none;
  margin: 0;
  padding: 0;
}}

.vep-menu li {{
  padding: 3px 0;
}}

.vep-menu .vep-menu-group {{
  margin-top: 12px;
  color: {Colour.SOFT_BLACK.var};
  font-weight: {Metric.FONT_WEIGHT_BOLD.value};
}}

.vep-menu .vep-menu-group:first-child {{
  margin-top: 0;
}}

.vep-content {{
  flex: 1;
  min-width: 0;
}}

.vep-section {{
  margin-bottom: {Metric.STANDARD_GUTTER.value};
}}

.vep-section h2 {{
  font-size: 15px;
  font-weight: {Metric.FONT_WEIGHT_NORMAL.value};
  color: {Colour.SOFT_BLACK.var};
  margin: 0 0 12px;
  padding-bottom: 4px;
  box-shadow: inset 0 -1.5px 0 {Colour.BLUE.var};
}}

.vep-figures {{
  display: flex;
  flex-wrap: wrap;
  align-items: flex-start;
  gap: {Metric.STANDARD_GUTTER.value};
}}

.vep-figures > svg {{
  max-width: 100%;
  height: auto;
}}
"""


@define(frozen=True, kw_only=True, auto_attribs=True)
class Section:
    """One anchored block of the page.

    Attributes:
        identifier (str): The anchor, and what the side menu links to.
        title (str): The heading, and the menu entry.
        blocks (tuple[str, ...]): Ready-made HTML fragments -- a table, a chart, a graph -- laid
            out side by side.
        group (str): Optional heading grouping consecutive entries in the side menu.
    """

    identifier: str
    title: str
    blocks: tuple[str, ...] = field(default=())
    group: str = field(default="")


def _menu(sections: Sequence[Section]) -> str:
    """Return the side menu linking to each section."""
    entries = []
    previous_group = ""
    for section in sections:
        if section.group and section.group != previous_group:
            entries.append(f'<li class="vep-menu-group">{escape(section.group)}</li>')
        entries.append(f'<li><a href="#{escape(section.identifier)}">{escape(section.title)}</a></li>')
        previous_group = section.group
    listing = "\n      ".join(entries)
    return f'<nav class="vep-menu">\n    <ul>\n      {listing}\n    </ul>\n  </nav>'


def _section(section: Section) -> str:
    """Return one section as HTML."""
    figures = "\n      ".join(section.blocks)
    return (
        f'<section class="vep-section" id="{escape(section.identifier)}">\n'
        f"    <h2>{escape(section.title)}</h2>\n"
        f'    <div class="vep-figures">\n      {figures}\n    </div>\n'
        f"  </section>"
    )


def render_page(title: str, sections: Sequence[Section], subtitle: str = "") -> str:
    """Return a complete HTML document.

    Args:
        title (str): The document title and masthead heading.
        sections (Sequence[Section]): The blocks to show, in order.
        subtitle (str): Optional line under the heading, e.g. the run that produced the page.

    Returns:
        str: A standalone document that makes no external requests.
    """
    styles = stylesheet() + TABLE_STYLES + LAYOUT_STYLES
    rendered_sections = "\n  ".join(_section(section) for section in sections)
    lead = f"<p>{escape(subtitle)}</p>" if subtitle else ""
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)}</title>
<style>
{styles}
</style>
</head>
<body>
  <header class="vep-masthead">
    <h1>{escape(title)}</h1>
    {lead}
  </header>
  <div class="vep-layout">
  {_menu(sections)}
  <main class="vep-content">
  {rendered_sections}
  </main>
  </div>
</body>
</html>
"""
