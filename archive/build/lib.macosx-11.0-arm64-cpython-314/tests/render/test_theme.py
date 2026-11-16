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
"""Tests for the Ensembl design tokens."""

from __future__ import annotations

import re

import pytest

from ensembl.vep.render.theme import Colour, Font, Metric, stylesheet

HEX_COLOUR = re.compile(r"^#(?:[0-9a-f]{3}|[0-9a-f]{6})$")


class TestColour:
    """The palette, transcribed from ensembl-client."""

    @pytest.mark.parametrize("colour", list(Colour), ids=lambda c: c.name.lower())
    def test_every_colour_is_a_lowercase_hex_value(self, colour: Colour) -> None:
        """A stylesheet is compared against upstream by eye, so casing must not drift."""
        assert HEX_COLOUR.match(colour.value)

    def test_the_token_follows_the_upstream_naming(self) -> None:
        """Refreshing the theme is a diff against one upstream file."""
        assert Colour.DARK_PINK.token == "--color-dark-pink"

    def test_the_var_reference_wraps_the_token(self) -> None:
        """What a rule in the stylesheet actually writes."""
        assert Colour.BLUE.var == "var(--color-blue)"

    def test_colours_are_distinct(self) -> None:
        """Two tokens sharing a value would make a chart ambiguous."""
        values = [colour.value for colour in Colour]
        assert len(values) == len(set(values))


class TestFontAndMetric:
    """The type and spacing tokens."""

    def test_the_body_stack_names_its_fallbacks(self) -> None:
        """The site self-hosts its fonts, so a report must degrade to what is installed."""
        assert Font.BODY.value.startswith("Lato, ")
        assert Font.BODY.value.endswith("sans-serif")

    def test_the_monospace_stack_ends_in_a_generic_family(self) -> None:
        """Same reason as the body stack."""
        assert Font.MONOSPACE.value.endswith("monospace")

    def test_font_tokens_are_namespaced(self) -> None:
        """Matching the upstream custom property."""
        assert Font.MONOSPACE.token == "--font-family-monospace"

    def test_metric_tokens_are_not_namespaced(self) -> None:
        """Upstream publishes them at the top level, e.g. --standard-gutter."""
        assert Metric.STANDARD_GUTTER.token == "--standard-gutter"

    @pytest.mark.parametrize("token", list(Font) + list(Metric), ids=lambda t: t.name.lower())
    def test_every_token_offers_a_var_reference(self, token: Font | Metric) -> None:
        """What a rule in the stylesheet writes, rather than the raw value."""
        assert token.var == f"var({token.token})"


class TestStylesheet:
    """The CSS a self-contained page inlines."""

    def test_every_colour_is_published_as_a_custom_property(self) -> None:
        """A rule referring to a token that was never emitted would silently render unstyled."""
        css = stylesheet()
        for colour in Colour:
            assert f"{colour.token}: {colour.value};" in css

    def test_every_font_and_metric_is_published(self) -> None:
        """Same reason as the colours."""
        css = stylesheet()
        for token in list(Font) + list(Metric):
            assert f"{token.token}: {token.value};" in css

    def test_it_declares_no_font_face(self) -> None:
        """Fonts are stacks only: linking the self-hosted woff2 would stop the page opening offline."""
        assert "@font-face" not in stylesheet()

    def test_it_declares_no_dark_theme(self) -> None:
        """Ensembl's design system is light-only; inventing one would not match the site."""
        css = stylesheet()
        assert "prefers-color-scheme" not in css
        assert "data-theme" not in css

    def test_the_body_rules_follow_the_site(self) -> None:
        """base.css sets the off-white ground, the black text and undecorated blue links."""
        css = stylesheet()
        assert f"background: {Colour.OFF_WHITE.var};" in css
        assert f"color: {Colour.BLACK.var};" in css
        assert "text-decoration: none;" in css
