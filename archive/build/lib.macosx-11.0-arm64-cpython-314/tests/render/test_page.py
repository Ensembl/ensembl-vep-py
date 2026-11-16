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
"""Tests for the self-contained HTML page."""

from __future__ import annotations

import re

import pytest

from ensembl.vep.render.page import Section, render_page

SECTIONS = (
    Section(identifier="classes", title="Variant classes", blocks=("<p>one</p>",)),
    Section(identifier="chromosomes", title="Variants by chromosome", blocks=("<p>two</p>",)),
)


@pytest.fixture(name="page")
def fixture_page() -> str:
    """A page with two sections and a subtitle."""
    return render_page("VEP summary", SECTIONS, subtitle="synthetic corpus")


class TestDocument:
    """The shape of the emitted document."""

    def test_it_is_a_complete_html_document(self, page: str) -> None:
        """A fragment would not open from a filesystem."""
        assert page.startswith("<!DOCTYPE html>")
        assert page.rstrip().endswith("</html>")

    def test_the_title_appears_in_the_head_and_the_masthead(self, page: str) -> None:
        """The browser tab and the page itself should agree."""
        assert "<title>VEP summary</title>" in page
        assert "<h1>VEP summary</h1>" in page

    def test_a_subtitle_is_shown_when_given(self, page: str) -> None:
        """Names the run that produced the page."""
        assert "synthetic corpus" in page

    def test_a_page_without_a_subtitle_omits_the_line(self) -> None:
        """An empty paragraph would leave a gap under the heading."""
        assert "<p></p>" not in render_page("VEP summary", SECTIONS)


class TestSections:
    """How sections and the menu relate."""

    def test_each_section_becomes_an_anchored_block(self, page: str) -> None:
        """The side menu links to these identifiers."""
        for section in SECTIONS:
            assert f'id="{section.identifier}"' in page

    def test_the_menu_links_to_every_section(self, page: str) -> None:
        """A menu entry with no target is a dead link."""
        for section in SECTIONS:
            assert f'href="#{section.identifier}"' in page

    def test_section_blocks_are_embedded_unescaped(self, page: str) -> None:
        """Blocks are ready-made markup -- a chart, a table, a graph."""
        assert "<p>one</p>" in page

    def test_headings_and_titles_are_escaped(self) -> None:
        """A title is text, unlike a block."""
        page = render_page("a & b", (Section(identifier="s", title="x < y"),))
        assert "a &amp; b" in page
        assert "x &lt; y" in page

    def test_a_page_with_no_sections_still_renders(self) -> None:
        """A run that reported nothing is a real answer."""
        assert "<main" in render_page("VEP summary", ())

    def test_consecutive_sections_can_share_a_menu_group(self) -> None:
        """A long report can expose its hierarchy without changing the section headings."""
        sections = (
            Section(identifier="one", title="One", group="Emitted"),
            Section(identifier="two", title="Two", group="Emitted"),
            Section(identifier="three", title="Three", group="Not fired"),
        )
        page = render_page("Grouped", sections)
        assert page.count('class="vep-menu-group"') == 2
        assert page.index("Emitted") < page.index('href="#one"')
        assert page.index("Not fired") < page.index('href="#three"')


class TestSelfContained:
    """The page must open offline."""

    def test_it_makes_no_external_requests(self, page: str) -> None:
        """The old summary loaded Google Charts from a CDN; this one loads nothing."""
        fetched = re.findall(r'(?:src|href)="(https?:)?//[^"]+"', page)
        assert not fetched

    def test_the_only_absolute_url_is_the_svg_namespace(self) -> None:
        """An xmlns is an identifier, not something the browser fetches."""
        page = render_page(
            "t",
            (Section(identifier="s", title="s", blocks=('<svg xmlns="http://www.w3.org/2000/svg"></svg>',)),),
        )
        urls = set(re.findall(r"https?://[^\"' ]+", page))
        assert urls <= {"http://www.w3.org/2000/svg"}

    def test_styles_are_inlined(self, page: str) -> None:
        """A linked stylesheet would be a second file to keep alongside the report."""
        assert "<style>" in page
        assert "stylesheet" not in page

    def test_the_theme_tokens_are_present(self, page: str) -> None:
        """The page carries the palette it refers to."""
        assert "--color-blue: #0099ff;" in page

    def test_the_table_styles_are_present(self, page: str) -> None:
        """A section embedding a table must not depend on a separate stylesheet."""
        assert ".vep-table" in page
