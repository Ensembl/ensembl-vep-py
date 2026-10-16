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
"""The formats a report can be rendered in, and which of them a request can combine.

One vocabulary for every reporting command, so ``--output_format`` means the same thing to ``explain``,
``catalogue``, ``corpus`` and ``bench``. A report is a table unless ``--graph`` asks for a diagram,
and a format belongs to one kind or the other: resolution refuses the combinations that cannot be
drawn rather than quietly rendering something else.
"""

from __future__ import annotations

__all__ = [
    "GRAPH_FORMATS",
    "REPORT_FORMATS",
    "TABLE_FORMATS",
    "resolve_report_format",
]

TABLE_FORMATS = ("text", "markdown", "confluence", "html")
"""Formats that render rows."""

GRAPH_FORMATS = ("svg", "dot", "mermaid")
"""Formats that render a diagram. ``svg`` is wrapped in a navigable HTML document."""

REPORT_FORMATS = TABLE_FORMATS + GRAPH_FORMATS
"""Every format ``--format`` accepts, whichever kind the command can render."""


def resolve_report_format(
    format_name: str | None,
    *,
    graph: bool,
    table_default: str,
    graph_default: str = "svg",
    graph_formats: tuple[str, ...] = GRAPH_FORMATS,
) -> str:
    """Return the format to render in, or say why the request cannot be met.

    Args:
        format_name (str | None): The requested format, or None to take the default for the kind
            of report asked for.
        graph (bool): Whether a diagram was asked for.
        table_default (str): The format a table takes when none was named.
        graph_default (str): The format a diagram takes when none was named.
        graph_formats (tuple[str, ...]): The diagram formats this report can draw; empty for a
            report that has no diagram form.

    Returns:
        str: One of ``REPORT_FORMATS``.

    Raises:
        ValueError: If a diagram was asked for in a table format or of a report that has none, or
            if a diagram format was named without asking for a diagram.
    """

    if graph and not graph_formats:
        raise ValueError("this report has no diagram form; drop --graph")
    if format_name is None:
        return graph_default if graph else table_default

    if graph:
        if format_name not in graph_formats:
            drawable = ", ".join(graph_formats)
            raise ValueError(f"--graph draws {drawable}, not {format_name!r}")
        return format_name

    if format_name in GRAPH_FORMATS:
        raise ValueError(f"{format_name!r} draws a diagram; pass --graph to ask for one")
    if format_name not in TABLE_FORMATS:
        known = ", ".join(REPORT_FORMATS)
        raise ValueError(f"unknown format {format_name!r}; expected one of {known}")
    return format_name
