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
"""Charts as inline SVG.

Replaces the Google Charts the old VEP summary loaded from a CDN, so a report opens offline. The
legend matches ``ensembl-client``'s legend component: a 10px square marker, 12px from its label.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from xml.sax.saxutils import escape

from attrs import define, field

from ensembl.vep.render.theme import Colour, Font

__all__ = ["Slice", "render_bars", "render_legend", "render_pie", "render_series"]

MARKER_SIZE = 10
"""Legend marker edge length, in pixels, matching the site's legend component."""
MARKER_GAP = 12
LEGEND_ROW_HEIGHT = 20
LABEL_SIZE = 11

PIE_DIAMETER = 200
BAR_HEIGHT = 18
BAR_GAP = 8
BAR_LABEL_WIDTH = 220
BAR_TRACK_WIDTH = 320
SERIES_HEIGHT = 120
MARGIN = 12

DEFAULT_COLOURS = (
    Colour.BLUE,
    Colour.DARK_PINK,
    Colour.DARK_YELLOW,
    Colour.LIME,
    Colour.TEAL,
    Colour.DUCKEGG_BLUE,
    Colour.PURPLE,
    Colour.ORANGE,
    Colour.GREY,
)
"""Fallback categorical palette, for a chart whose slices name no colour of their own."""


@define(frozen=True, kw_only=True, auto_attribs=True)
class Slice:
    """One labelled quantity in a chart.

    Attributes:
        label (str): What the quantity counts.
        value (float): How much.
        colour (Colour | None): The colour to draw it in, or None to take one from the default
            palette by position.
    """

    label: str
    value: float
    colour: Colour | None = field(default=None)


def _colour_at(slices: Sequence[Slice], index: int) -> Colour:
    """Return the colour a slice is drawn in, falling back to the palette by position."""
    declared = slices[index].colour
    if declared is not None:
        return declared
    return DEFAULT_COLOURS[index % len(DEFAULT_COLOURS)]


def _svg(width: int, height: int, body: str, title: str = "") -> str:
    """Wrap chart body markup in a sized ``<svg>`` element."""
    heading = f"<title>{escape(title)}</title>" if title else ""
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="{width}" height="{height}" font-family="{Font.BODY.value}">{heading}\n  '
        f"{body}\n</svg>"
    )


def render_legend(slices: Sequence[Slice]) -> str:
    """Return a legend listing each slice with its marker."""
    rows = []
    for index, item in enumerate(slices):
        y = MARGIN + index * LEGEND_ROW_HEIGHT
        colour = _colour_at(slices, index)
        rows.append(
            f'<rect x="{MARGIN}" y="{y}" width="{MARKER_SIZE}" height="{MARKER_SIZE}" '
            f'fill="{colour.value}"/>'
        )
        rows.append(
            f'<text x="{MARGIN + MARKER_SIZE + MARKER_GAP}" y="{y + MARKER_SIZE}" '
            f'font-size="{LABEL_SIZE}" fill="{Colour.BLACK.value}">'
            f"{escape(item.label)} ({item.value:g})</text>"
        )
    height = MARGIN * 2 + max(len(slices), 1) * LEGEND_ROW_HEIGHT
    return _svg(320, height, "\n  ".join(rows))


def _arc_path(centre: float, radius: float, start: float, end: float) -> str:
    """Return the path of a pie wedge between two angles, in radians."""
    start_x = centre + radius * math.cos(start)
    start_y = centre + radius * math.sin(start)
    end_x = centre + radius * math.cos(end)
    end_y = centre + radius * math.sin(end)
    large_arc = 1 if end - start > math.pi else 0
    return (
        f"M {centre:.2f} {centre:.2f} L {start_x:.2f} {start_y:.2f} "
        f"A {radius:.2f} {radius:.2f} 0 {large_arc} 1 {end_x:.2f} {end_y:.2f} Z"
    )


def render_pie(slices: Sequence[Slice], title: str = "") -> str:
    """Return a pie chart.

    A single slice is drawn as a whole circle, because a wedge spanning the full turn has
    coincident start and end points and would render as nothing.
    """
    total = sum(item.value for item in slices)
    radius = PIE_DIAMETER / 2
    centre = radius + MARGIN
    size = PIE_DIAMETER + MARGIN * 2

    if total <= 0:
        body = f'<circle cx="{centre}" cy="{centre}" r="{radius}" fill="{Colour.LIGHT_GREY.value}"/>'
        return _svg(size, size, body, title)

    drawn = []
    if len(slices) == 1:
        drawn.append(
            f'<circle cx="{centre}" cy="{centre}" r="{radius}" fill="{_colour_at(slices, 0).value}"/>'
        )
    else:
        angle = -math.pi / 2
        for index, item in enumerate(slices):
            sweep = 2 * math.pi * (item.value / total)
            path = _arc_path(centre, radius, angle, angle + sweep)
            drawn.append(f'<path d="{path}" fill="{_colour_at(slices, index).value}"/>')
            angle += sweep
    return _svg(size, size, "\n  ".join(drawn), title)


def render_bars(slices: Sequence[Slice], title: str = "") -> str:
    """Return a horizontal bar chart, one labelled bar per slice."""
    largest = max((item.value for item in slices), default=0)
    rows = []
    for index, item in enumerate(slices):
        y = MARGIN + index * (BAR_HEIGHT + BAR_GAP)
        length = 0 if largest <= 0 else int(BAR_TRACK_WIDTH * item.value / largest)
        rows.append(
            f'<text x="{MARGIN}" y="{y + BAR_HEIGHT - 5}" font-size="{LABEL_SIZE}" '
            f'fill="{Colour.BLACK.value}">{escape(item.label)}</text>'
        )
        rows.append(
            f'<rect x="{MARGIN + BAR_LABEL_WIDTH}" y="{y}" width="{length}" height="{BAR_HEIGHT}" '
            f'fill="{_colour_at(slices, index).value}"/>'
        )
        rows.append(
            f'<text x="{MARGIN + BAR_LABEL_WIDTH + length + 6}" y="{y + BAR_HEIGHT - 5}" '
            f'font-size="{LABEL_SIZE}" fill="{Colour.DARK_GREY.value}">{item.value:g}</text>'
        )
    width = MARGIN * 2 + BAR_LABEL_WIDTH + BAR_TRACK_WIDTH + 60
    height = MARGIN * 2 + max(len(slices), 1) * (BAR_HEIGHT + BAR_GAP)
    return _svg(width, height, "\n  ".join(rows), title)


def render_series(values: Sequence[float], title: str = "", colour: Colour = Colour.BLUE) -> str:
    """Return a filled area chart of one series, for a positional distribution."""
    width = max(len(values), 1) * 3 + MARGIN * 2
    height = SERIES_HEIGHT + MARGIN * 2
    largest = max(values, default=0)

    if not values or largest <= 0:
        baseline = f'<line x1="{MARGIN}" y1="{height - MARGIN}" x2="{width - MARGIN}" '
        baseline += f'y2="{height - MARGIN}" stroke="{Colour.GREY.value}" stroke-width="1"/>'
        return _svg(width, height, baseline, title)

    points = [f"{MARGIN},{height - MARGIN}"]
    for index, value in enumerate(values):
        x = MARGIN + index * 3
        y = height - MARGIN - int(SERIES_HEIGHT * value / largest)
        points.append(f"{x},{y}")
    points.append(f"{MARGIN + (len(values) - 1) * 3},{height - MARGIN}")
    body = f'<polygon points="{" ".join(points)}" fill="{colour.value}" fill-opacity="0.75"/>'
    return _svg(width, height, body, title)
