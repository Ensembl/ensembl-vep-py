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
"""Placing a layered graph on a canvas.

These graphs are small -- tens of facts and rules -- and already carry the layer each node belongs
to, so laying them out is arithmetic rather than a general graph-drawing problem. Layers become
columns, and a node's position within its column follows declaration order.
"""

from __future__ import annotations

from attrs import define, field

from ensembl.vep.render.graph.model import Graph, Node

__all__ = ["Layout", "Placement", "layered", "tree"]

LAYER_GAP = 220
"""Horizontal distance between the left edges of two adjacent layers, in pixels."""
NODE_GAP = 18
"""Vertical distance between two nodes in the same layer, in pixels."""
NODE_WIDTH = 170
NODE_HEIGHT = 34
OPERATOR_WIDTH = 64
RULE_WIDTH = 220
MARGIN = 24
"""Space left around the drawing, in pixels."""


@define(frozen=True, kw_only=True, auto_attribs=True)
class Placement:
    """Where one node sits on the canvas.

    Attributes:
        node (Node): The node placed.
        x (int): Left edge, in pixels.
        y (int): Top edge, in pixels.
        width (int): Box width, in pixels.
        height (int): Box height, in pixels.
    """

    node: Node
    x: int
    y: int
    width: int = field(default=NODE_WIDTH)
    height: int = field(default=NODE_HEIGHT)

    @property
    def centre_y(self) -> int:
        """Return the vertical midpoint, where an edge attaches."""
        return self.y + self.height // 2

    @property
    def right(self) -> int:
        """Return the right edge, where an outgoing edge leaves."""
        return self.x + self.width


@define(frozen=True, kw_only=True, auto_attribs=True)
class Layout:
    """A graph with a position for every node.

    Attributes:
        graph (Graph): The graph laid out.
        placements (tuple[Placement, ...]): One entry per node, in the graph's own order.
        width (int): Canvas width, in pixels.
        height (int): Canvas height, in pixels.
    """

    graph: Graph
    placements: tuple[Placement, ...]
    width: int
    height: int

    def placement(self, identifier: str) -> Placement:
        """Return the placement of one node.

        Raises:
            KeyError: If no node has that identifier.
        """
        for placement in self.placements:
            if placement.node.identifier == identifier:
                return placement
        raise KeyError(f"Unknown node '{identifier}'")


def _node_width(node: Node) -> int:
    """Return the width assigned to a node kind."""
    if node.kind == "operator":
        return OPERATOR_WIDTH
    if node.kind == "rule":
        return RULE_WIDTH
    return NODE_WIDTH


def _column_x(column: int, width: int) -> int:
    """Return a node's left edge, centred in its layer column."""
    return MARGIN + column * LAYER_GAP + (NODE_WIDTH - width) // 2


def _canvas_width(placements: list[Placement]) -> int:
    """Return the width covering every placed node."""
    if not placements:
        return MARGIN * 2 + NODE_WIDTH
    return max(placement.right for placement in placements) + MARGIN


def layered(graph: Graph) -> Layout:
    """Lay a graph out in columns, one per layer.

    Layers are taken from the nodes themselves rather than inferred, because the domain already
    knows the answer: a fact's layer is its cost tier, a rule's is its rank.
    """
    layers = graph.layers()
    placements: list[Placement] = []
    tallest_column = 0

    for column, layer in enumerate(layers):
        nodes = graph.nodes_in_layer(layer)
        for row, node in enumerate(nodes):
            width = _node_width(node)
            x = _column_x(column, width)
            y = MARGIN + row * (NODE_HEIGHT + NODE_GAP)
            placements.append(Placement(node=node, x=x, y=y, width=width))
        column_height = len(nodes) * (NODE_HEIGHT + NODE_GAP) - NODE_GAP
        tallest_column = max(tallest_column, column_height)

    width = _canvas_width(placements)
    height = MARGIN * 2 + max(tallest_column, NODE_HEIGHT)
    return Layout(graph=graph, placements=tuple(placements), width=width, height=height)


def tree(graph: Graph) -> Layout:
    """Lay a child-to-parent tree out with every parent centred on its children.

    Edges must point from an operand towards its parent. The one node with no outgoing edge is the
    root. Leaves keep graph declaration order, making the result stable while preserving sibling
    order from the source expression.

    Raises:
        ValueError: If the graph is disconnected, cyclic, or gives a node more than one parent.
    """
    if not graph.nodes:
        return layered(graph)

    identifiers = {node.identifier for node in graph.nodes}
    children: dict[str, list[str]] = {identifier: [] for identifier in identifiers}
    parents: dict[str, str] = {}
    for edge in graph.edges:
        if edge.source not in identifiers or edge.target not in identifiers:
            continue
        if edge.source in parents:
            raise ValueError(f"Node '{edge.source}' has more than one parent")
        parents[edge.source] = edge.target
        children[edge.target].append(edge.source)

    roots = [node.identifier for node in graph.nodes if node.identifier not in parents]
    if len(roots) != 1:
        raise ValueError("A tree must have exactly one root")

    centres: dict[str, int] = {}
    visiting: set[str] = set()
    next_leaf_centre = MARGIN + NODE_HEIGHT // 2

    def place_centre(identifier: str) -> int:
        nonlocal next_leaf_centre
        if identifier in centres:
            return centres[identifier]
        if identifier in visiting:
            raise ValueError("A tree cannot contain a cycle")
        visiting.add(identifier)

        child_identifiers = children[identifier]
        if child_identifiers:
            child_centres = [place_centre(child) for child in child_identifiers]
            centre = (child_centres[0] + child_centres[-1]) // 2
        else:
            centre = next_leaf_centre
            next_leaf_centre += NODE_HEIGHT + NODE_GAP

        visiting.remove(identifier)
        centres[identifier] = centre
        return centre

    place_centre(roots[0])
    if len(centres) != len(graph.nodes):
        raise ValueError("A tree must be connected")

    columns = {layer: column for column, layer in enumerate(graph.layers())}
    placements: list[Placement] = []
    for node in graph.nodes:
        width = _node_width(node)
        x = _column_x(columns[node.layer], width)
        y = centres[node.identifier] - NODE_HEIGHT // 2
        placements.append(Placement(node=node, x=x, y=y, width=width))

    width = _canvas_width(placements)
    height = max(placement.y + placement.height for placement in placements) + MARGIN
    return Layout(graph=graph, placements=tuple(placements), width=width, height=height)
