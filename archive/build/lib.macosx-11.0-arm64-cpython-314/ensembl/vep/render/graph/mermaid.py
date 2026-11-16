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
"""Rendering a graph as Mermaid flowchart source."""

from __future__ import annotations

import re

from ensembl.vep.render.graph.model import Graph

__all__ = ["render_mermaid"]

_UNSAFE = re.compile(r"[^0-9A-Za-z_]")


def _identifier(raw: str) -> str:
    """Return an identifier Mermaid accepts unquoted."""
    return _UNSAFE.sub("_", raw)


def _label(text: str) -> str:
    """Return label text with the characters Mermaid treats as syntax removed."""
    return text.replace('"', "'").replace("[", "(").replace("]", ")")


def render_mermaid(graph: Graph) -> str:
    """Return a graph as a left-to-right Mermaid flowchart.

    Node state becomes a Mermaid class rather than a colour, so the diagram picks up whatever
    theme the surrounding document defines.
    """
    lines = ["flowchart LR"]

    for node in graph.nodes:
        identifier = _identifier(node.identifier)
        lines.append(f'  {identifier}["{_label(node.label)}"]')
        if node.state:
            lines.append(f"  class {identifier} {_identifier(node.state)};")

    for edge in graph.edges:
        source = _identifier(edge.source)
        target = _identifier(edge.target)
        if edge.label:
            lines.append(f"  {source} -->|{_label(edge.label)}| {target}")
        else:
            lines.append(f"  {source} --> {target}")

    return "\n".join(lines)
