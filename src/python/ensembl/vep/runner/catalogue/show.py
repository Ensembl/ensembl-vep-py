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
"""Consequence catalogue and fact registry rendering."""

from __future__ import annotations

import logging
from collections.abc import Callable
from typing import Any, ClassVar

from ensembl.vep.common import CostTier
from ensembl.vep.annotators.consequence.facts.base import Fact, REGISTRY
from ensembl.vep.annotators.consequence.graphs import (
    catalogue_overview_graph,
    fact_dependency_graph,
    rule_graph,
    rule_graphs,
)
from ensembl.vep.annotators.consequence.rules.catalogue import (
    ConsequenceDefinition,
    load_consequence_definitions,
)
from ensembl.vep.annotators.consequence.rules.expressions import expr_from_json
from ensembl.vep.config import Config
from ensembl.vep.config.report import resolve_report_format
from ensembl.vep.render.graph import Edge, Graph, Node, render_dot, render_mermaid, render_svg, tree
from ensembl.vep.render.page import Section, render_page
from ensembl.vep.render.table import TableFormat, render_table
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.registry import register_runner

__all__ = ["Runner"]

logger = logging.getLogger(__name__)

TERM_HEADERS = ("Rank", "Term", "Accession", "Impact", "When")
FACT_HEADERS = ("Fact", "Cost tier", "Deps", "Requires", "Description")
NO_EXPRESSION = "(no when -- excluded from the annotators)"

GRAPH_EMITTERS = {"dot": render_dot, "mermaid": render_mermaid}
CATALOGUE_OVERVIEW_ID = "catalogue-overview"

SORTABLE_BY_SUBJECT = {"terms": ("name", "impact"), "facts": ("name", "cost_tier")}
"""Which sort keys mean anything for which subject."""


def _description(registered_fact: Fact) -> str:
    """Return the first line of a fact's documentation."""
    for line in (registered_fact.documentation or "").splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return ""


def _term_sort_key(sort: str | None) -> Callable[[ConsequenceDefinition], Any] | None:
    """Return the key a term listing is ordered by, or None for catalogue rank."""
    if sort == "name":
        return lambda definition: definition.term
    if sort == "impact":
        return lambda definition: (definition.impact.value if definition.impact else "", definition.rank)
    return None


def _fact_sort_key(sort: str | None) -> Callable[[tuple[str, Fact]], Any]:
    """Return the key a fact listing is ordered by, or None for cost tier then name."""
    if sort == "name":
        return lambda item: item[0]
    return lambda item: (item[1].cost, item[0])


@register_runner("catalogue.show")
class Runner:
    """Runner for rendering catalogue terms or facts."""

    name: ClassVar[str] = "catalogue.show"
    sections: ClassVar[tuple[str, ...]] = ("logging", "catalogue", "catalogue_show")
    description: ClassVar[str] = "list catalogue terms or registered facts"

    def __init__(self, config: Config) -> None:
        """Hold the configuration."""
        self.config = config

    def run(self) -> ExitCode:
        """Render the requested subject in the requested form."""
        show = self.config.catalogue_show
        self._warn_if_sort_does_not_apply()
        try:
            report_format = resolve_report_format(
                show.output_format, graph=show.graph, table_default="markdown"
            )
            rendered = self._graph(report_format) if show.graph else self._table(report_format)
        except ValueError as unusable:
            logger.error("%s", unusable)
            return ExitCode.USAGE
        except KeyError as unknown:
            logger.error("%s", unknown.args[0])
            return ExitCode.USAGE

        self._emit(rendered)
        return ExitCode.OK

    def _warn_if_sort_does_not_apply(self) -> None:
        """Say so when a sort key means nothing for the chosen subject."""
        show = self.config.catalogue_show
        if show.sort and show.sort not in SORTABLE_BY_SUBJECT[show.subject]:
            logger.warning("--sort %s does not apply to %s; using the default order", show.sort, show.subject)

    def _table(self, table_format: str) -> str:
        """Render the chosen subject as a table."""
        show = self.config.catalogue_show
        if show.subject == "facts":
            headers, rows = FACT_HEADERS, self._fact_rows()
        else:
            headers, rows = TERM_HEADERS, self._term_rows()
        logger.info("%d %s", len(rows), show.subject)
        return render_table(headers, rows, TableFormat(table_format))

    def _graph(self, graph_format: str) -> str:
        """Render the chosen subject as one or more graphs.

        Unselected SVG terms become a navigable HTML document. Text graph formats remain one
        expression after another so callers can process their source directly.
        """
        if graph_format == "svg":
            return self._svg()

        emitter = GRAPH_EMITTERS[graph_format]
        rendered = []
        for graph in self._graphs():
            self._log_graph(graph)
            rendered.append(emitter(graph))
        return "\n\n".join(rendered)

    def _svg(self) -> str:
        """Render a selected graph, or the navigable term catalogue."""
        show = self.config.catalogue_show
        if show.subject == "terms" and show.select is None:
            return self._catalogue_page()

        rendered = []
        for graph in self._graphs():
            self._log_graph(graph)
            layout = tree(graph) if show.subject == "terms" else None
            rendered.append(render_svg(graph, layout))
        return "\n\n".join(rendered)

    def _catalogue_page(self) -> str:
        """Return every term expression and the catalogue overview as an HTML document."""
        definitions = load_consequence_definitions()
        graphs = {graph.title: graph for graph in rule_graphs()}
        sections: list[Section] = []

        for definition in definitions:
            graph = graphs.get(definition.term)
            if graph is None:
                block = "<p>This fallback term has no <code>when</code> expression.</p>"
            else:
                self._log_graph(graph)
                block = render_svg(graph, tree(graph))
            sections.append(Section(identifier=definition.term, title=definition.term, blocks=(block,)))

        overview = catalogue_overview_graph(definitions)
        self._log_graph(overview)
        links = {definition.term: f"#{definition.term}" for definition in definitions}
        sections.append(
            Section(
                identifier=CATALOGUE_OVERVIEW_ID,
                title=overview.title,
                blocks=(render_svg(overview, links=links),),
            )
        )

        subtitle = f"{len(graphs)} rule expressions and {len(definitions) - len(graphs)} fallback terms"
        return render_page("VEP consequence catalogue", sections, subtitle=subtitle)

    @staticmethod
    def _log_graph(graph: Graph) -> None:
        """Log the size of one rendered graph."""
        logger.info("%s: %d nodes, %d edges", graph.title, len(graph.nodes), len(graph.edges))

    def _graphs(self) -> tuple[Graph, ...]:
        """Return the graphs to draw, in catalogue order."""
        show = self.config.catalogue_show
        if show.subject == "facts":
            graph = fact_dependency_graph()
            return (_neighbourhood(graph, show.select) if show.select else graph,)
        if show.select:
            return (rule_graph(show.select),)
        return rule_graphs()

    def _term_rows(self) -> list[tuple[object, ...]]:
        """Return one row per catalogue term, filtered and ordered as configured."""
        show = self.config.catalogue_show
        definitions = list(load_consequence_definitions())
        if show.select:
            definitions = [definition for definition in definitions if definition.term == show.select]
            if not definitions:
                raise KeyError(f"no catalogue term named {show.select!r}")

        sort_key = _term_sort_key(show.sort if show.sort in SORTABLE_BY_SUBJECT["terms"] else None)
        if sort_key is not None:
            definitions.sort(key=sort_key)

        rows: list[tuple[object, ...]] = []
        for definition in definitions:
            rows.append(_term_row(definition))
        return rows

    def _fact_rows(self) -> list[tuple[object, ...]]:
        """Return one row per registered fact, filtered and ordered as configured."""
        show = self.config.catalogue_show
        facts = list(REGISTRY.items())
        if show.select:
            facts = [item for item in facts if item[0] == show.select]
            if not facts:
                raise KeyError(f"no registered fact named {show.select!r}")

        facts.sort(key=_fact_sort_key(show.sort if show.sort in SORTABLE_BY_SUBJECT["facts"] else None))

        rows: list[tuple[object, ...]] = []
        for name, registered_fact in facts:
            rows.append(
                (
                    name,
                    CostTier.label(registered_fact.cost),
                    ", ".join(registered_fact.deps),
                    ", ".join(registered_fact.requires),
                    _description(registered_fact),
                )
            )
        return rows

    def _emit(self, rendered: str) -> None:
        """Write the rendering to the configured destination."""
        destination = self.config.catalogue_show.output_file
        if destination is None:
            print(rendered)
            return
        destination.write_text(rendered, encoding="utf-8")
        logger.info("wrote %s", destination)


def _term_row(definition: ConsequenceDefinition) -> tuple[object, ...]:
    """Return one catalogue term as a table row."""
    formula = NO_EXPRESSION if definition.when is None else repr(expr_from_json(definition.when))
    impact = definition.impact.value if definition.impact else None
    return (definition.rank, definition.term, definition.accession, impact, formula)


def _neighbourhood(graph: Graph, identifier: str) -> Graph:
    """Return the part of a graph around one node: itself and whatever it connects to directly.

    Raises:
        KeyError: If no node has that identifier.
    """
    graph.node(identifier)

    kept = {identifier}
    edges: list[Edge] = []
    for edge in graph.edges:
        if identifier in (edge.source, edge.target):
            edges.append(edge)
            kept.add(edge.source)
            kept.add(edge.target)

    nodes: list[Node] = []
    for node in graph.nodes:
        if node.identifier in kept:
            nodes.append(node)
    return Graph(nodes=tuple(nodes), edges=tuple(edges), title=f"{graph.title}: {identifier}")
