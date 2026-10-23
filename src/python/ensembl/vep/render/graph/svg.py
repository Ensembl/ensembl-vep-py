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
"""Rendering a laid-out graph as inline SVG.

The one place node semantics become an appearance: a builder says a node was forced false, this
module decides that means a dashed grey box.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from xml.sax.saxutils import escape, quoteattr

from ensembl.vep.render.graph.layout import Layout, layered
from ensembl.vep.render.graph.model import Graph
from ensembl.vep.render.theme import Colour, Font

__all__ = ["STATE_COLOURS", "render_svg"]

STATE_COLOURS = {
    "operator": Colour.LIGHT_GREY,
    "term": Colour.LIGHT_BLUE,
    "subsumed": Colour.DUCKEGG_BLUE,
    "suppressor": Colour.DARK_PINK,
    "computed": Colour.ICE_BLUE,
    "prefilled": Colour.DUCKEGG_BLUE,
    "cached": Colour.NEON_BLUE,
    "forced_false": Colour.LIGHT_GREY,
    "fired": Colour.LIGHT_BLUE,
    "pruned": Colour.MEDIUM_LIGHT_GREY,
    "emitted": Colour.BLUE,
    "true": Colour.LIGHT_BLUE,
    "false": Colour.RED,
    "unreached": Colour.WHITE,
    "impact_high": Colour.DARK_PINK,
    "impact_moderate": Colour.ORANGE,
    "impact_low": Colour.DUCKEGG_BLUE,
    "impact_modifier": Colour.LIGHT_GREY,
}
"""Node state -> fill colour. A state absent from here falls back to the neutral fill.

The one place a builder's semantics become an appearance: ``subsumed`` and ``suppressor`` are the
two directions of suppression, and they are coloured apart so which way it runs is readable without
following the arrowheads.
"""

DEFAULT_FILL = Colour.WHITE
DEFAULT_TEXT = Colour.BLACK
STATE_TEXT_COLOURS = {"false": Colour.WHITE}
STROKE = Colour.DARK_GREY
EDGE_COLOUR = Colour.GREY
LABEL_SIZE = 11
DETAIL_SIZE = 9
_UNSAFE_IDENTIFIER = re.compile(r"[^0-9A-Za-z_-]")


def _fill_for(state: str) -> Colour:
    """Return the fill a node state is drawn with."""
    return STATE_COLOURS.get(state, DEFAULT_FILL)


def _dashed(state: str) -> bool:
    """Return whether a node state is drawn with a dashed outline."""
    return state in {"forced_false", "unreached"}


def _text_for(state: str) -> Colour:
    """Return a readable text colour for a node state."""
    return STATE_TEXT_COLOURS.get(state, DEFAULT_TEXT)


def _node_svg(layout: Layout, links: Mapping[str, str]) -> list[str]:
    """Return one SVG group per placed node."""
    parts = []
    for placement in layout.placements:
        node = placement.node
        text_colour = _text_for(node.state)
        dash = ' stroke-dasharray="4 3"' if _dashed(node.state) else ""
        radius = placement.height // 2 if node.kind == "operator" else 0
        stroke_width = 2 if node.kind == "rule" else 1
        node_parts = [
            f'<rect x="{placement.x}" y="{placement.y}" width="{placement.width}" '
            f'height="{placement.height}" fill="{_fill_for(node.state).value}" '
            f'stroke="{STROKE.value}" stroke-width="{stroke_width}" rx="{radius}"{dash}/>'
        ]
        text_x = placement.x + 8
        node_parts.append(
            f'<text x="{text_x}" y="{placement.y + 15}" font-size="{LABEL_SIZE}" '
            f'fill="{text_colour.value}">{escape(node.label)}</text>'
        )
        if node.detail:
            node_parts.append(
                f'<text x="{text_x}" y="{placement.y + 27}" font-size="{DETAIL_SIZE}" '
                f'fill="{text_colour.value}">{escape(node.detail)}</text>'
            )
        body = "".join(node_parts)
        if node.identifier in links:
            parts.append(f"<a href={quoteattr(links[node.identifier])}>{body}</a>")
        else:
            parts.append(body)
    return parts


def _edge_svg(layout: Layout, marker_identifier: str) -> list[str]:
    """Return one SVG path per edge whose endpoints were both placed."""
    parts = []
    for edge in layout.graph.edges:
        try:
            source = layout.placement(edge.source)
            target = layout.placement(edge.target)
        except KeyError:
            continue
        start_x = source.right
        start_y = source.centre_y
        end_x = target.x
        end_y = target.centre_y
        control = (start_x + end_x) // 2
        parts.append(
            f'<path d="M {start_x} {start_y} C {control} {start_y}, {control} {end_y}, '
            f'{end_x} {end_y}" fill="none" stroke="{EDGE_COLOUR.value}" stroke-width="1" '
            f'marker-end="url(#{marker_identifier})"/>'
        )
        if edge.label:
            label_x = (start_x + end_x) // 2
            label_y = (start_y + end_y) // 2 - 4
            parts.append(
                f'<text x="{label_x}" y="{label_y}" text-anchor="middle" '
                f'font-size="{DETAIL_SIZE}" fill="{Colour.DARK_GREY.value}">{escape(edge.label)}</text>'
            )
    return parts


def _marker_identifier(title: str) -> str:
    """Return an SVG identifier unique to a titled graph on an HTML page."""
    suffix = _UNSAFE_IDENTIFIER.sub("-", title).strip("-") or "graph"
    return f"arrow-{suffix}"


def render_svg(
    graph: Graph,
    layout: Layout | None = None,
    *,
    links: Mapping[str, str] | None = None,
) -> str:
    """Return a graph as a standalone ``<svg>`` element.

    Args:
        graph (Graph): The graph to draw.
        layout (Layout | None): A layout to use, or None to lay the graph out in columns.
        links (Mapping[str, str] | None): Optional node identifiers to link destinations.

    Returns:
        str: The SVG markup, self-contained apart from the font stack it names.
    """
    placed = layout or layered(graph)
    marker_identifier = _marker_identifier(graph.title)
    marker = (
        f'<defs><marker id="{marker_identifier}" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="6" '
        f'markerHeight="6" orient="auto"><path d="M 0 0 L 8 4 L 0 8 z" '
        f'fill="{EDGE_COLOUR.value}"/></marker></defs>'
    )
    body = "\n  ".join(_edge_svg(placed, marker_identifier) + _node_svg(placed, links or {}))
    title = f"<title>{escape(graph.title)}</title>" if graph.title else ""
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {placed.width} {placed.height}" '
        f'width="{placed.width}" height="{placed.height}" '
        f'font-family="{Font.BODY.value}">{title}\n  {marker}\n  {body}\n</svg>'
    )
