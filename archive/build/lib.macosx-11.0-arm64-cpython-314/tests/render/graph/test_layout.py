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
"""Tests for laying a graph out in columns."""

from __future__ import annotations

import pytest

from ensembl.vep.render.graph.layout import LAYER_GAP, MARGIN, NODE_GAP, NODE_HEIGHT, layered, tree
from ensembl.vep.render.graph.model import Edge, Graph, Node


@pytest.fixture(name="graph")
def fixture_graph() -> Graph:
    """Two layers, with two nodes in the second."""
    return Graph(
        nodes=(
            Node(identifier="a", label="a", layer=0),
            Node(identifier="b", label="b", layer=1),
            Node(identifier="c", label="c", layer=1),
        )
    )


class TestPlacement:
    """Where each node lands."""

    def test_every_node_is_placed(self, graph: Graph) -> None:
        """An unplaced node would be dropped from the drawing."""
        assert len(layered(graph).placements) == len(graph.nodes)

    def test_layers_become_columns(self, graph: Graph) -> None:
        """The horizontal axis is the domain's own ordering: cost tier, or rank."""
        layout = layered(graph)
        assert layout.placement("a").x == MARGIN
        assert layout.placement("b").x == MARGIN + LAYER_GAP

    def test_nodes_in_a_layer_stack_downwards(self, graph: Graph) -> None:
        """Row order follows declaration order, so output is reproducible."""
        layout = layered(graph)
        assert layout.placement("b").y == MARGIN
        assert layout.placement("c").y == MARGIN + NODE_HEIGHT + NODE_GAP

    def test_an_unknown_identifier_is_rejected(self, graph: Graph) -> None:
        """An edge referring to an unplaced node must be detectable."""
        with pytest.raises(KeyError):
            layered(graph).placement("absent")


class TestAttachmentPoints:
    """Where an edge meets a node."""

    def test_an_edge_leaves_the_right_edge(self, graph: Graph) -> None:
        """Layout is left to right."""
        placement = layered(graph).placement("a")
        assert placement.right == placement.x + placement.width

    def test_an_edge_attaches_at_the_vertical_midpoint(self, graph: Graph) -> None:
        """So a connector meets the box centrally rather than at a corner."""
        placement = layered(graph).placement("a")
        assert placement.centre_y == placement.y + placement.height // 2


class TestCanvas:
    """The size the drawing needs."""

    def test_the_canvas_covers_every_node(self, graph: Graph) -> None:
        """A node outside the viewBox would be clipped."""
        layout = layered(graph)
        for placement in layout.placements:
            assert placement.right <= layout.width
            assert placement.y + placement.height <= layout.height

    def test_an_empty_graph_still_has_a_canvas(self) -> None:
        """A zero-sized viewBox renders as nothing at all."""
        layout = layered(Graph())
        assert layout.width > 0
        assert layout.height > 0


class TestTree:
    """Tree-aware vertical placement."""

    def test_a_parent_is_centred_on_its_children(self) -> None:
        """The shape of a Boolean branch remains visible across columns."""
        graph = Graph(
            nodes=(
                Node(identifier="a", label="a", kind="atom", layer=0),
                Node(identifier="b", label="b", kind="atom", layer=0),
                Node(identifier="all", label="ALL", kind="operator", layer=1),
                Node(identifier="term", label="term", kind="rule", layer=2),
            ),
            edges=(
                Edge(source="a", target="all"),
                Edge(source="b", target="all"),
                Edge(source="all", target="term"),
            ),
        )
        layout = tree(graph)
        child_midpoint = (layout.placement("a").centre_y + layout.placement("b").centre_y) // 2
        assert layout.placement("all").centre_y == child_midpoint
        assert layout.placement("term").centre_y == child_midpoint

    def test_operator_nodes_are_compact(self) -> None:
        """ALL and ANY act as gates rather than full-width fact boxes."""
        graph = Graph(
            nodes=(
                Node(identifier="atom", label="atom", kind="atom", layer=0),
                Node(identifier="operator", label="ALL", kind="operator", layer=1),
            ),
            edges=(Edge(source="atom", target="operator"),),
        )
        layout = tree(graph)
        assert layout.placement("operator").width < layout.placement("atom").width

    def test_more_than_one_root_is_rejected(self) -> None:
        """A disconnected graph cannot be positioned as one expression tree."""
        graph = Graph(nodes=(Node(identifier="a", label="a"), Node(identifier="b", label="b")))
        with pytest.raises(ValueError, match="exactly one root"):
            tree(graph)
