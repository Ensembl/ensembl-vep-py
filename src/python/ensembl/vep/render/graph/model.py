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
"""A directed graph carrying semantics, not styling.

A builder describes what a node *is* -- its kind, which layer it belongs to, what state it ended
in -- and never what colour it should be. Renderers map those semantics onto an appearance, so
every styling decision stays in one package and a builder needs no knowledge of the theme.
"""

from __future__ import annotations

from attrs import define, field

__all__ = ["Edge", "Graph", "Node"]


@define(frozen=True, kw_only=True, auto_attribs=True)
class Node:
    """One vertex.

    Attributes:
        identifier (str): Unique within a graph, and what edges refer to.
        label (str): The text a renderer draws.
        kind (str): What the node represents, e.g. "fact", "rule", "term".
        layer (int): Which column or row the node belongs to. Layered layouts read this directly;
            it comes from something meaningful in the domain, such as a fact's cost tier or a
            rule's rank.
        state (str): How the node ended up, e.g. "computed", "forced_false", "fired".
        hits (int): How many times the node was reached, for a corpus-reach overlay.
        detail (str): Optional secondary text, e.g. a value or a reason.
    """

    identifier: str
    label: str
    kind: str = field(default="")
    layer: int = field(default=0)
    state: str = field(default="")
    hits: int = field(default=0)
    detail: str = field(default="")


@define(frozen=True, kw_only=True, auto_attribs=True)
class Edge:
    """One directed connection between two nodes.

    Attributes:
        source (str): The identifier the edge leaves.
        target (str): The identifier the edge enters.
        kind (str): What the connection means, e.g. "depends", "requires", "subsumes".
        label (str): Optional text drawn along the edge.
    """

    source: str
    target: str
    kind: str = field(default="")
    label: str = field(default="")


@define(frozen=True, kw_only=True, auto_attribs=True)
class Graph:
    """A set of nodes and the edges between them.

    Attributes:
        nodes (tuple[Node, ...]): The vertices, in the order a builder produced them.
        edges (tuple[Edge, ...]): The connections.
        title (str): What the graph shows.
    """

    nodes: tuple[Node, ...] = field(default=())
    edges: tuple[Edge, ...] = field(default=())
    title: str = field(default="")

    def node(self, identifier: str) -> Node:
        """Return one node by its identifier.

        Raises:
            KeyError: If no node has that identifier.
        """
        for node in self.nodes:
            if node.identifier == identifier:
                return node
        raise KeyError(f"Unknown node '{identifier}'")

    def layers(self) -> tuple[int, ...]:
        """Return the distinct layer numbers present, in ascending order."""
        present = {node.layer for node in self.nodes}
        return tuple(sorted(present))

    def nodes_in_layer(self, layer: int) -> tuple[Node, ...]:
        """Return the nodes of one layer, in declaration order."""
        selected = []
        for node in self.nodes:
            if node.layer == layer:
                selected.append(node)
        return tuple(selected)

    def edges_from(self, identifier: str) -> tuple[Edge, ...]:
        """Return the edges leaving one node."""
        leaving = []
        for edge in self.edges:
            if edge.source == identifier:
                leaving.append(edge)
        return tuple(leaving)

    def with_hits(self, hits_by_identifier: dict[str, int]) -> Graph:
        """Return a copy whose nodes carry the given hit counts.

        How a corpus-reach overlay is applied: the shape of the graph is unchanged, only what each
        node reports about how often it was reached.
        """
        updated = []
        for node in self.nodes:
            hits = hits_by_identifier.get(node.identifier, 0)
            updated.append(
                Node(
                    identifier=node.identifier,
                    label=node.label,
                    kind=node.kind,
                    layer=node.layer,
                    state=node.state,
                    hits=hits,
                    detail=node.detail,
                )
            )
        return Graph(nodes=tuple(updated), edges=self.edges, title=self.title)
