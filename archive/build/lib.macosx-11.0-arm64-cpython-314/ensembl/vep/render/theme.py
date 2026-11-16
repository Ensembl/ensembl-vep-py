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
"""Ensembl's design tokens.

Transcribed from the frontend source of truth, ``Ensembl/ensembl-client``
(``src/styles/design-tokens.css``), so refreshing the theme is a diff against one upstream file.

The design system is light-only: there is no ``prefers-color-scheme`` block, no ``data-theme`` and
no dark token set anywhere upstream, so emitting one would invent a theme Ensembl does not have.
"""

from __future__ import annotations

from enum import Enum

__all__ = ["Colour", "Font", "Metric", "stylesheet"]


def _custom_property(prefix: str, member_name: str) -> str:
    """Return the CSS custom property a token member is published as."""
    return f"--{prefix}{member_name.lower().replace('_', '-')}"


class Colour(Enum):
    """The Ensembl palette, by its upstream token name."""

    BLACK = "#1b2c39"
    SOFT_BLACK = "#374148"
    WHITE = "#fff"
    OFF_WHITE = "#fefefe"

    BLUE = "#0099ff"
    LIGHT_BLUE = "#33adff"
    ICE_BLUE = "#e5f5ff"
    TEAL = "#327c89"
    DUCKEGG_BLUE = "#96d0c9"
    NEON_BLUE = "#8ef4f8"

    LIGHT_GREY = "#f1f2f4"
    MEDIUM_LIGHT_GREY = "#d4d9de"
    GREY = "#b7c0c8"
    MEDIUM_DARK_GREY = "#9aa7b1"
    DARK_GREY = "#6f8190"

    RED = "#d90000"
    DARK_PINK = "#eb768a"
    PURPLE = "#ca81ff"
    ORANGE = "#ff9900"
    MUSTARD = "#cc9933"
    DARK_YELLOW = "#f8c041"
    GREEN = "#47d147"
    LIME = "#84fa3a"

    @property
    def token(self) -> str:
        """Return the CSS custom property this colour is published as."""
        return _custom_property("color-", self.name)

    @property
    def var(self) -> str:
        """Return a ``var()`` reference to this colour, for use in a stylesheet."""
        return f"var({self.token})"


class Font(Enum):
    """The two font stacks, with their upstream fallbacks."""

    BODY = "Lato, 'Helvetica Neue', Helvetica, Roboto, Arial, sans-serif"
    MONOSPACE = "'IBM Plex Mono', 'Liberation Mono', Courier, monospace"

    @property
    def token(self) -> str:
        """Return the CSS custom property this stack is published as."""
        return _custom_property("font-family-", self.name)

    @property
    def var(self) -> str:
        """Return a ``var()`` reference to this stack, for use in a stylesheet."""
        return f"var({self.token})"


class Metric(Enum):
    """Spacing and type metrics, by their upstream token name."""

    STANDARD_GUTTER = "30px"
    BODY_FONT_SIZE = "13px"
    BODY_LINE_HEIGHT = "1.5"
    FONT_WEIGHT_LIGHT = "300"
    FONT_WEIGHT_NORMAL = "400"
    FONT_WEIGHT_BOLD = "700"

    @property
    def token(self) -> str:
        """Return the CSS custom property this metric is published as."""
        return _custom_property("", self.name)

    @property
    def var(self) -> str:
        """Return a ``var()`` reference to this metric, for use in a stylesheet."""
        return f"var({self.token})"


def stylesheet() -> str:
    """Return the theme as CSS, for inlining into a self-contained page.

    Emits the tokens as custom properties on ``:root`` followed by the handful of element rules
    ``ensembl-client``'s ``base.css`` sets, so a report matches the site without linking to it.
    Fonts are stacks only, with no ``@font-face``: the site self-hosts its woff2 files, and linking
    them would stop the page opening offline.
    """
    declarations = []
    for colour in Colour:
        declarations.append(f"  {colour.token}: {colour.value};")
    for font in Font:
        declarations.append(f"  {font.token}: {font.value};")
    for metric in Metric:
        declarations.append(f"  {metric.token}: {metric.value};")
    tokens = "\n".join(declarations)

    return f""":root {{
{tokens}
}}

*,
*::before,
*::after {{
  box-sizing: border-box;
}}

html {{
  font-size: {Metric.BODY_FONT_SIZE.value};
}}

body {{
  margin: 0;
  padding: 0;
  background: {Colour.OFF_WHITE.var};
  color: {Colour.BLACK.var};
  font-family: {Font.BODY.var};
  line-height: {Metric.BODY_LINE_HEIGHT.value};
  -webkit-font-smoothing: antialiased;
}}

a {{
  color: {Colour.BLUE.var};
  text-decoration: none;
}}

code, kbd, pre {{
  font-family: {Font.MONOSPACE.var};
}}
"""
