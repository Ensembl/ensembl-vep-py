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
"""Tests for the graph model."""

from __future__ import annotations

import pytest

from ensembl.vep.render.graph.model import Edge, Graph, Node


@pytest.fixture(name="graph")
def fixture_graph() -> Graph:
    """A three-layer graph with one node per layer and one branch."""
    return Graph(
        title="derivation",
        nodes=(
            Node(identifier="is_snv", label="is_snv", kind="fact", layer=0, state="computed"),
            Node(identifier="overlap_cds", label="overlap_cds", kind="fact", layer=1, state="forced_false"),
            Node(identifier="overlap_exon", label="overlap_exon", kind="fact", layer=1, state="computed"),
            Node(identifier="intron_variant", label="intron_variant", kind="rule", layer=2, state="fired"),
        ),
        edges=(
            Edge(source="is_snv", target="overlap_cds", kind="requires"),
            Edge(source="overlap_exon", target="intron_variant", kind="atom"),
        ),
    )


class TestNodeLookup:
    """Finding a node by its identifier."""

    def test_a_known_identifier_returns_its_node(self, graph: Graph) -> None:
        """Edges refer to nodes by identifier, so the lookup must agree."""
        assert graph.node("overlap_cds").state == "forced_false"

    def test_an_unknown_identifier_is_rejected(self, graph: Graph) -> None:
        """Silently returning nothing would render an edge to nowhere."""
        with pytest.raises(KeyError):
            graph.node("absent")


class TestLayers:
    """Layers come from the domain, not from inference."""

    def test_layers_are_reported_in_ascending_order(self, graph: Graph) -> None:
        """Columns are drawn left to right in this order."""
        assert graph.layers() == (0, 1, 2)

    def test_a_layer_returns_its_nodes_in_declaration_order(self, graph: Graph) -> None:
        """Row order within a column follows the builder, so output is reproducible."""
        labels = [node.label for node in graph.nodes_in_layer(1)]
        assert labels == ["overlap_cds", "overlap_exon"]

    def test_an_empty_layer_returns_nothing(self, graph: Graph) -> None:
        """Asking about a layer no node occupies is not an error."""
        assert not graph.nodes_in_layer(9)

    def test_an_empty_graph_has_no_layers(self) -> None:
        """A graph nothing was added to is a real answer."""
        assert not Graph().layers()


class TestEdges:
    """Reading the edges leaving a node."""

    def test_only_the_edges_leaving_the_node_are_returned(self, graph: Graph) -> None:
        """Direction is what makes the graph a derivation rather than an association."""
        leaving = graph.edges_from("is_snv")
        assert [edge.target for edge in leaving] == ["overlap_cds"]

    def test_a_node_with_no_outgoing_edges_returns_nothing(self, graph: Graph) -> None:
        """A terminal node is normal, not an error."""
        assert not graph.edges_from("intron_variant")


class TestHitsOverlay:
    """Applying corpus-reach counts."""

    def test_named_nodes_carry_their_count(self, graph: Graph) -> None:
        """How a dead rule becomes visible rather than only listed."""
        overlaid = graph.with_hits({"overlap_cds": 12})
        assert overlaid.node("overlap_cds").hits == 12

    def test_unnamed_nodes_report_no_hits(self, graph: Graph) -> None:
        """Never reached is the point of the overlay, so it must be zero rather than absent."""
        overlaid = graph.with_hits({"overlap_cds": 12})
        assert overlaid.node("is_snv").hits == 0

    def test_the_shape_of_the_graph_is_unchanged(self, graph: Graph) -> None:
        """Only what each node reports changes."""
        overlaid = graph.with_hits({"overlap_cds": 12})
        assert len(overlaid.nodes) == len(graph.nodes)
        assert overlaid.edges == graph.edges
        assert overlaid.title == graph.title

    def test_the_original_is_left_alone(self, graph: Graph) -> None:
        """The model is frozen; an overlay returns a copy."""
        graph.with_hits({"is_snv": 5})
        assert graph.node("is_snv").hits == 0
