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
"""Tests for chart rendering."""

from __future__ import annotations

import pytest

from ensembl.vep.render.chart import (
    DEFAULT_COLOURS,
    MARKER_SIZE,
    Slice,
    render_bars,
    render_legend,
    render_pie,
    render_series,
)
from ensembl.vep.render.theme import Colour

SLICES = (
    Slice(label="intron_variant", value=564),
    Slice(label="upstream_gene_variant", value=52),
    Slice(label="feature_elongation", value=8),
)


class TestPie:
    """The chart the summary uses for variant classes and consequences."""

    def test_it_draws_one_wedge_per_slice(self) -> None:
        """A missing wedge would silently under-report a category."""
        assert render_pie(SLICES).count("<path") == len(SLICES)

    def test_a_single_slice_is_drawn_as_a_circle(self) -> None:
        """A wedge spanning a full turn has coincident endpoints and would render as nothing."""
        svg = render_pie((Slice(label="only", value=5),))
        assert "<circle" in svg
        assert "<path" not in svg

    def test_no_slices_yields_an_empty_dial(self) -> None:
        """A run with nothing to report is a real answer."""
        assert "<circle" in render_pie(())

    def test_slices_with_no_value_yield_an_empty_dial(self) -> None:
        """Guards the division by the total."""
        assert "<circle" in render_pie((Slice(label="none", value=0),))

    def test_a_declared_colour_is_used(self) -> None:
        """Consequence charts pass Ensembl's own group colours."""
        svg = render_pie((Slice(label="a", value=1, colour=Colour.TEAL), Slice(label="b", value=1)))
        assert Colour.TEAL.value in svg

    def test_an_undeclared_colour_comes_from_the_palette(self) -> None:
        """So a chart of arbitrary categories is still legible."""
        assert DEFAULT_COLOURS[0].value in render_pie(SLICES)


class TestBars:
    """The chart the summary uses for variants by chromosome."""

    def test_it_draws_one_bar_per_slice(self) -> None:
        """One row each, label and track."""
        assert render_bars(SLICES).count("<rect") == len(SLICES)

    def test_the_largest_value_fills_the_track(self) -> None:
        """Bars are scaled against the largest, not against an absolute width."""
        svg = render_bars((Slice(label="a", value=10), Slice(label="b", value=5)))
        assert 'width="320"' in svg
        assert 'width="160"' in svg

    def test_values_of_zero_render_without_dividing(self) -> None:
        """A chromosome with no variants is a real row."""
        assert 'width="0"' in render_bars((Slice(label="a", value=0),))

    def test_labels_are_escaped(self) -> None:
        """A label may carry characters that would close a tag."""
        assert "&lt;" in render_bars((Slice(label="a<b", value=1),))


class TestSeries:
    """The area chart the summary uses for positional distribution."""

    def test_it_draws_a_filled_polygon(self) -> None:
        """The shape the old page drew with a Google area chart."""
        assert "<polygon" in render_series((11, 0, 7, 3))

    def test_an_empty_series_draws_a_baseline(self) -> None:
        """A chromosome with no variants still gets a row rather than a gap."""
        assert "<line" in render_series(())

    def test_an_all_zero_series_draws_a_baseline(self) -> None:
        """Guards the division by the largest value."""
        assert "<line" in render_series((0, 0, 0))


class TestLegend:
    """The key shared by the charts."""

    def test_it_draws_one_marker_per_slice(self) -> None:
        """A legend shorter than its chart would mislabel a category."""
        assert render_legend(SLICES).count("<rect") == len(SLICES)

    def test_markers_match_the_site_component(self) -> None:
        """ensembl-client's legend uses a 10px square."""
        assert f'width="{MARKER_SIZE}" height="{MARKER_SIZE}"' in render_legend(SLICES)

    def test_each_entry_shows_its_value(self) -> None:
        """The old summary listed counts beside the key."""
        assert "(564)" in render_legend(SLICES)


class TestEveryChart:
    """What all the charts do the same way."""

    @pytest.mark.parametrize(
        "markup",
        [render_pie(SLICES), render_bars(SLICES), render_legend(SLICES), render_series((1, 2))],
    )
    def test_output_is_a_sized_svg_element(self, markup: str) -> None:
        """The page embeds these inline, so each must carry its own dimensions."""
        assert markup.startswith("<svg ")
        assert "viewBox=" in markup
        assert markup.rstrip().endswith("</svg>")

    def test_a_title_becomes_an_svg_title(self) -> None:
        """Gives the figure an accessible name without a separate caption."""
        assert "<title>Consequences</title>" in render_pie(SLICES, "Consequences")
