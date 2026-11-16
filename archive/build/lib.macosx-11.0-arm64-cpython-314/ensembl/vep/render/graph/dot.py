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
"""Rendering a graph as Graphviz DOT."""

from __future__ import annotations

from ensembl.vep.render.graph.model import Graph

__all__ = ["render_dot"]


def _quote(text: str) -> str:
    """Return text escaped for a DOT double-quoted string."""
    return text.replace("\\", "\\\\").replace('"', '\\"')


def render_dot(graph: Graph) -> str:
    """Return a graph as a DOT digraph.

    Node semantics are emitted as attributes rather than as an appearance, so a caller styling the
    output with their own Graphviz theme has everything the builder knew.
    """
    lines = [f'digraph "{_quote(graph.title)}" {{', "  rankdir=LR;", "  node [shape=box];"]

    for layer in graph.layers():
        members = []
        for node in graph.nodes_in_layer(layer):
            members.append(f'"{_quote(node.identifier)}"')
        lines.append(f"  {{ rank=same; {' '.join(members)} }}")

    for node in graph.nodes:
        attributes = [f'label="{_quote(node.label)}"']
        if node.kind:
            attributes.append(f'kind="{_quote(node.kind)}"')
        if node.state:
            attributes.append(f'state="{_quote(node.state)}"')
        if node.hits:
            attributes.append(f'hits="{node.hits}"')
        lines.append(f'  "{_quote(node.identifier)}" [{", ".join(attributes)}];')

    for edge in graph.edges:
        attributes = []
        if edge.label:
            attributes.append(f'label="{_quote(edge.label)}"')
        if edge.kind:
            attributes.append(f'kind="{_quote(edge.kind)}"')
        suffix = f' [{", ".join(attributes)}]' if attributes else ""
        lines.append(f'  "{_quote(edge.source)}" -> "{_quote(edge.target)}"{suffix};')

    lines.append("}")
    return "\n".join(lines)
