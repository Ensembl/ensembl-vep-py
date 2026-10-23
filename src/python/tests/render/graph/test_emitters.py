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
"""Tests for the graph emitters."""

from __future__ import annotations

from collections.abc import Callable

import pytest

from ensembl.vep.render.graph.dot import render_dot
from ensembl.vep.render.graph.mermaid import render_mermaid
from ensembl.vep.render.graph.model import Edge, Graph, Node
from ensembl.vep.render.graph.svg import STATE_COLOURS, render_svg


@pytest.fixture(name="graph")
def fixture_graph() -> Graph:
    """One computed fact feeding a fired rule, plus a forced-false fact."""
    return Graph(
        title="derivation",
        nodes=(
            Node(identifier="overlap_exon", label="overlap_exon", kind="fact", layer=0, state="computed"),
            Node(
                identifier="overlap_cds",
                label="overlap_cds",
                kind="fact",
                layer=0,
                state="forced_false",
                detail="via overlap_exon",
            ),
            Node(identifier="intron_variant", label="intron_variant", kind="rule", layer=1, state="fired"),
        ),
        edges=(Edge(source="overlap_exon", target="intron_variant", kind="atom", label="atom"),),
    )


class TestSvg:
    """Inline SVG, the format the page embeds."""

    def test_every_node_becomes_a_box(self, graph: Graph) -> None:
        """A missing box would drop a step from the derivation."""
        assert render_svg(graph).count("<rect") == len(graph.nodes)

    def test_every_edge_becomes_a_connector(self, graph: Graph) -> None:
        """Counted by the arrow marker, since the marker definition is itself a path."""
        assert render_svg(graph).count('marker-end="url(#arrow-derivation)"') == len(graph.edges)

    def test_an_edge_label_is_drawn(self, graph: Graph) -> None:
        """The connection meaning should not disappear in SVG output."""
        assert ">atom</text>" in render_svg(graph)

    def test_state_decides_the_fill(self, graph: Graph) -> None:
        """The one place semantics become an appearance."""
        assert STATE_COLOURS["computed"].value in render_svg(graph)

    @pytest.mark.parametrize("state", ["true", "false"])
    def test_boolean_states_have_distinct_fills(self, state: str) -> None:
        """An explanation uses blue and red to distinguish satisfied and blocked conditions."""
        graph = Graph(nodes=(Node(identifier="condition", label="condition", state=state),))
        assert f'fill="{STATE_COLOURS[state].value}"' in render_svg(graph)

    def test_false_nodes_use_readable_text(self) -> None:
        """White text retains contrast against the red false-state fill."""
        graph = Graph(nodes=(Node(identifier="condition", label="condition", state="false"),))
        assert 'fill="#fff">condition</text>' in render_svg(graph)

    def test_a_forced_false_node_is_dashed(self, graph: Graph) -> None:
        """It was never computed, so it should not read as a solid step."""
        assert "stroke-dasharray" in render_svg(graph)

    def test_detail_text_is_drawn(self, graph: Graph) -> None:
        """The reason a fact was forced false is the interesting part."""
        assert "via overlap_exon" in render_svg(graph)

    def test_an_edge_to_an_unplaced_node_is_skipped(self) -> None:
        """A builder may reference a node it chose not to draw."""
        graph = Graph(nodes=(Node(identifier="a", label="a"),), edges=(Edge(source="a", target="absent"),))
        assert "marker-end=" not in render_svg(graph)

    def test_labels_are_escaped(self) -> None:
        """A rule expression contains characters that would close a tag."""
        graph = Graph(nodes=(Node(identifier="r", label="a & b"),))
        assert "a &amp; b" in render_svg(graph)

    def test_a_node_can_link_to_an_html_section(self, graph: Graph) -> None:
        """The catalogue overview links back to each detailed term graph."""
        svg = render_svg(graph, links={"intron_variant": "#intron_variant"})
        assert '<a href="#intron_variant">' in svg

    def test_each_titled_graph_has_its_own_marker(self) -> None:
        """Several inline SVGs must not publish duplicate arrow identifiers."""
        first = render_svg(Graph(title="first"))
        second = render_svg(Graph(title="second"))
        assert 'id="arrow-first"' in first
        assert 'id="arrow-second"' in second


class TestDot:
    """Graphviz source, for a caller with their own styling."""

    def test_it_emits_a_left_to_right_digraph(self, graph: Graph) -> None:
        """Matching the layered layout the other emitters use."""
        dot = render_dot(graph)
        assert dot.startswith('digraph "derivation" {')
        assert "rankdir=LR;" in dot

    def test_layers_become_rank_groups(self, graph: Graph) -> None:
        """So Graphviz reproduces the domain's own columns."""
        assert render_dot(graph).count("rank=same") == len(graph.layers())

    def test_semantics_are_emitted_as_attributes(self, graph: Graph) -> None:
        """A caller styling the output has everything the builder knew."""
        dot = render_dot(graph)
        assert 'kind="fact"' in dot
        assert 'state="forced_false"' in dot

    def test_no_colour_is_emitted(self, graph: Graph) -> None:
        """Styling belongs to the renderer, and DOT output is styled by its consumer."""
        assert "fillcolor" not in render_dot(graph)

    def test_an_edge_carrying_neither_label_nor_kind_gets_no_attributes(self) -> None:
        """An empty attribute list would emit a bare `[]`, which DOT rejects."""
        graph = Graph(
            nodes=(Node(identifier="a", label="a"), Node(identifier="b", label="b")),
            edges=(Edge(source="a", target="b"),),
        )
        assert '"a" -> "b";' in render_dot(graph)

    def test_hit_counts_are_emitted_when_an_overlay_was_applied(self, graph: Graph) -> None:
        """A corpus-reach overlay must survive into the DOT a caller styles."""
        overlaid = graph.with_hits({"overlap_cds": 12})
        assert 'hits="12"' in render_dot(overlaid)

    def test_a_node_never_reached_emits_no_hit_count(self, graph: Graph) -> None:
        """Zero is the default, so emitting it would add noise to every unoverlaid graph."""
        assert "hits=" not in render_dot(graph)

    def test_quotes_in_a_label_are_escaped(self) -> None:
        """A DOT string ends at the first unescaped quote."""
        graph = Graph(nodes=(Node(identifier="a", label='say "hi"'),))
        assert r"\"hi\"" in render_dot(graph)


class TestMermaid:
    """Mermaid source, for embedding in a markdown document."""

    def test_it_emits_a_left_to_right_flowchart(self, graph: Graph) -> None:
        """Matching the other emitters."""
        assert render_mermaid(graph).startswith("flowchart LR")

    def test_every_node_appears(self, graph: Graph) -> None:
        """One declaration each."""
        mermaid = render_mermaid(graph)
        for node in graph.nodes:
            assert f'{node.identifier}["{node.label}"]' in mermaid

    def test_state_becomes_a_class(self, graph: Graph) -> None:
        """So the diagram picks up the surrounding document's theme."""
        assert "class overlap_cds forced_false;" in render_mermaid(graph)

    def test_an_edge_label_is_carried(self, graph: Graph) -> None:
        """Names what the connection means."""
        assert "-->|atom|" in render_mermaid(graph)

    def test_an_unlabelled_edge_is_a_plain_arrow(self) -> None:
        """Most edges carry no label, and Mermaid has no empty-label form."""
        graph = Graph(
            nodes=(Node(identifier="a", label="a"), Node(identifier="b", label="b")),
            edges=(Edge(source="a", target="b"),),
        )
        assert "  a --> b" in render_mermaid(graph)

    def test_identifiers_are_made_safe(self) -> None:
        """A fact name may contain characters Mermaid reads as syntax."""
        graph = Graph(nodes=(Node(identifier="a.b-c", label="a.b-c"),))
        assert "a_b_c[" in render_mermaid(graph)

    def test_brackets_in_a_label_are_replaced(self) -> None:
        """Square brackets would close the node declaration."""
        graph = Graph(nodes=(Node(identifier="a", label="x[0]"),))
        assert "x(0)" in render_mermaid(graph)


class TestEveryEmitter:
    """What all three do the same way."""

    @pytest.mark.parametrize("emitter", [render_svg, render_dot, render_mermaid])
    def test_an_empty_graph_renders(self, emitter: Callable[[Graph], str]) -> None:
        """A graph nothing matched is a real answer, not an error."""
        assert emitter(Graph())
