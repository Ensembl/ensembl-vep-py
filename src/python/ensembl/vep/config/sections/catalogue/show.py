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
"""Options for listing what the catalogue and fact registry contain."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from attrs import define, field

from ensembl.vep.config.option_spec import option
from ensembl.vep.config.sections.base import ConfigSection

__all__ = ["CatalogueShowConfig"]

SUBJECTS = ("terms", "facts")
"""What can be listed: the catalogue's SO terms, or the fact registry behind them."""

SORT_KEYS = ("name", "impact", "cost_tier")
"""Every order any subject can be listed in.

``name`` applies to both. ``impact`` is a term's property and ``cost_tier`` a fact's, so asking for
one on the other subject falls back to that subject's natural order rather than failing.
"""

TABLE_FORMATS = ("markdown", "confluence", "text", "html")
GRAPH_FORMATS = ("svg", "dot", "mermaid")
"""Kept apart from the table formats: a table format cannot draw a graph, and vice versa."""


@define(frozen=True, kw_only=True, auto_attribs=True)
class CatalogueShowConfig(ConfigSection):
    """What to list, and how to render it.

    Attributes:
        subject (str): Whether to list the catalogue's terms or the registered facts.
        select (str | None): Show only the term or fact with this name, or None for all of them.
            In graph mode it keeps that node and whatever it connects to directly.
        graph (str | None): Draw the subject as a graph in this format, or None to render a table.
            Unselected SVG terms are embedded in a linked HTML catalogue.
        format (str): Table format. Only meaningful when ``graph`` is not given.
        sort (str | None): Row order, or None for the subject's natural order -- catalogue rank for
            terms, cost tier for facts.
        output_file (Path | None): Where to write, or None for stdout.
    """

    section_name: ClassVar[str] = "catalogue_show"

    subject: str = field(default="terms", metadata=option("what to list", choices=SUBJECTS))
    select: str | None = field(default=None, metadata=option("show only this term or fact"))
    graph: str | None = field(
        default=None, metadata=option("draw as a graph in this format", choices=GRAPH_FORMATS)
    )
    format: str = field(default="markdown", metadata=option("table format", choices=TABLE_FORMATS))
    sort: str | None = field(default=None, metadata=option("row order", choices=SORT_KEYS))
    output_file: Path | None = field(default=None, metadata=option("write output here"))
