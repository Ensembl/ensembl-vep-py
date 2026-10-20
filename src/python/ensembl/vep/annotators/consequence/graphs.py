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
"""Build graphs from consequence facts and rules.

This module assigns domain meaning to graph nodes and edges. Renderers remain responsible for
visual presentation so ``ensembl.vep.render`` stays independent of consequence annotation.
"""

from __future__ import annotations

from collections.abc import Iterable

from ensembl.vep.common import CostTier
from ensembl.vep.annotators.consequence.facts.base import REGISTRY
from ensembl.vep.annotators.consequence.rules.catalogue import (
    ConsequenceDefinition,
    load_consequence_definitions,
)
from ensembl.vep.annotators.consequence.rules.expressions import And, Atom, Expr, Not, Or
from ensembl.vep.annotators.consequence.rules.rulesets import FactRule, FactRuleSet, build_fact_ruleset
from ensembl.vep.render.graph import Edge, Graph, Node
from ensembl.vep.trace import AnnotationTrace

__all__ = [
    "catalogue_overview_graph",
    "derivation_graphs",
    "fact_dependency_graph",
    "rule_graph",
    "rule_graphs",
]

FACT_KIND = "fact"
RULE_KIND = "rule"
ATOM_KIND = "atom"
OPERATOR_KIND = "operator"
TERM_KIND = "term"

DEPENDS_EDGE = "depends"
REQUIRES_EDGE = "requires"
OPERAND_EDGE = "operand"
WHEN_EDGE = "when"
SUBSUMES_EDGE = "subsumes"


def _tier_layer(cost: int) -> int:
    """Return the column a fact belongs in, by its declared cost tier.

    Tiers are the domain's own left-to-right ordering -- cheapest first -- so the layout needs no
    inference of its own.
    """
    tiers = sorted(CostTier)
    for position, tier in enumerate(tiers):
        if cost <= tier:
            return position
    return len(tiers)


def fact_dependency_graph() -> Graph:
    """Return the fact registry as a graph, columns ordered by cost tier.

    Two edge kinds, and the distinction matters: ``deps`` is what must be computed first, while
    ``requires`` is the implication a staging planner exploits to settle a rule without computing
    the expensive atom at all.
    """
    nodes = []
    edges = []
    for name, registered_fact in sorted(REGISTRY.items()):
        nodes.append(
            Node(
                identifier=name,
                label=name,
                kind=FACT_KIND,
                layer=_tier_layer(registered_fact.cost),
                detail=CostTier.label(registered_fact.cost),
            )
        )
        for dependency in registered_fact.deps:
            edges.append(Edge(source=dependency, target=name, kind=DEPENDS_EDGE))
        for requirement in registered_fact.requires:
            edges.append(Edge(source=requirement, target=name, kind=REQUIRES_EDGE))

    return Graph(nodes=tuple(nodes), edges=tuple(edges), title="Fact dependencies")


def _operator(expression: Expr) -> tuple[str, tuple[Expr, ...]]:
    """Return an operator node's label and its operands, spelled as the catalogue JSON spells them."""
    if isinstance(expression, Not):
        return "NOT", (expression.child,)
    if isinstance(expression, And):
        return "ALL", expression.children
    if isinstance(expression, Or):
        return "ANY", expression.children
    return type(expression).__name__.lower(), ()


def _truth_state(expression: Expr, known: dict[str, bool] | None) -> str:
    """Return the state of an expression whose observed atoms are in ``known``."""
    if known is None:
        return ""

    value = expression.evaluate_partial(known)
    if value is True:
        return "true"
    if value is False:
        return "false"
    return ""


def _expression_tree(
    expression: Expr,
    prefix: str,
    depth: int,
    known: dict[str, bool] | None = None,
) -> tuple[list[Node], list[Edge], str, int]:
    """Build one formula as a tree, and return its nodes, edges, root and deepest level.

    Identifiers embed the path to the node, so an atom appearing in two branches of the same
    formula gets a node in each. That duplication is what keeps a tree readable: a shared atom node
    would pull unrelated branches together and turn the tree back into a tangle.
    """
    if isinstance(expression, Atom):
        leaf = Node(
            identifier=prefix,
            label=expression.name,
            kind=ATOM_KIND,
            layer=depth,
            state=_truth_state(expression, known),
        )
        return [leaf], [], prefix, depth

    if isinstance(expression, Not) and isinstance(expression.child, Atom):
        leaf = Node(
            identifier=prefix,
            label=f"NOT {expression.child.name}",
            kind=ATOM_KIND,
            layer=depth,
            state=_truth_state(expression, known),
            detail="must be false" if known is None else "",
        )
        return [leaf], [], prefix, depth

    label, children = _operator(expression)
    nodes = [
        Node(
            identifier=prefix,
            label=label,
            kind=OPERATOR_KIND,
            layer=depth,
            state=_truth_state(expression, known) or ("operator" if known is None else ""),
        )
    ]
    edges: list[Edge] = []
    deepest = depth

    for index, child in enumerate(children):
        child_nodes, child_edges, child_root, child_deepest = _expression_tree(
            child,
            f"{prefix}.{index}",
            depth + 1,
            known,
        )
        nodes.extend(child_nodes)
        edges.extend(child_edges)
        edges.append(Edge(source=child_root, target=prefix, kind=OPERAND_EDGE))
        deepest = max(deepest, child_deepest)

    return nodes, edges, prefix, deepest


def _rule_detail(rule: FactRule) -> str:
    """Return the metadata shown under a consequence term."""
    parts = [rule.accession, f"rank {rule.rank}", f"tier {rule.tier}"]
    if rule.impact is not None:
        parts.append(rule.impact.value)
    return " · ".join(parts)


def _rule_graph(
    rule: FactRule,
    *,
    known: dict[str, bool] | None = None,
    term_state: str = "term",
    term_detail: str | None = None,
) -> Graph:
    """Return one rule tree, optionally decorated with one overlap's observed truth values."""
    term = rule.term
    nodes, edges, root, deepest = _expression_tree(rule.when, f"{term}:when", 0, known)

    # The tree is built root-first; flip it so the deepest atoms sit in the leftmost column.
    flipped = []
    for node in nodes:
        flipped.append(
            Node(
                identifier=node.identifier,
                label=node.label,
                kind=node.kind,
                layer=deepest - node.layer,
                state=node.state,
                detail=node.detail,
            )
        )

    flipped.append(
        Node(
            identifier=term,
            label=term,
            kind=RULE_KIND,
            layer=deepest + 1,
            state=term_state,
            detail=term_detail if term_detail is not None else _rule_detail(rule),
        )
    )
    edges.append(Edge(source=root, target=term, kind=WHEN_EDGE, label="when"))

    return Graph(nodes=tuple(flipped), edges=tuple(edges), title=term)


def rule_graph(term: str, ruleset: FactRuleSet | None = None) -> Graph:
    """Return one rule as a Boolean expression tree.

    The formula reads left to right: atoms feed ``ALL`` and ``ANY`` operators, which feed the term.
    Negation of an atom stays on that leaf; a separate ``NOT`` operator is reserved for a compound
    expression. Subsumption belongs to ``catalogue_overview_graph`` rather than this tree.

    Args:
        term (str): The SO term to draw.
        ruleset (FactRuleSet | None): The rules to read, or None for the packaged catalogue's.

    Returns:
        Graph: The term's tree.

    Raises:
        KeyError: If no rule is named ``term``.
    """
    rules = ruleset or build_fact_ruleset()
    if term not in rules.by_term:
        raise KeyError(f"no catalogue rule named {term!r}")
    return _rule_graph(rules.by_term[term])


def _subsumption_layers(definitions: tuple[ConsequenceDefinition, ...]) -> dict[str, int]:
    """Return left-to-right layers following the direction of suppression."""
    known_terms = {definition.term for definition in definitions}
    outgoing: dict[str, list[str]] = {definition.term: [] for definition in definitions}
    incoming_count = {definition.term: 0 for definition in definitions}

    for definition in definitions:
        for suppressed in definition.subsumes:
            if suppressed not in known_terms:
                continue
            outgoing[definition.term].append(suppressed)
            incoming_count[suppressed] += 1

    layers = {definition.term: 0 for definition in definitions}
    pending = [definition.term for definition in definitions if incoming_count[definition.term] == 0]
    visited = 0
    for term in pending:
        visited += 1
        for suppressed in outgoing[term]:
            layers[suppressed] = max(layers[suppressed], layers[term] + 1)
            incoming_count[suppressed] -= 1
            if incoming_count[suppressed] == 0:
                pending.append(suppressed)

    if visited != len(definitions):
        raise ValueError("Cannot draw a catalogue whose subsumption graph contains a cycle")
    return layers


def catalogue_overview_graph(
    definitions: Iterable[ConsequenceDefinition] | None = None,
) -> Graph:
    """Return every catalogue term and the direction of subsumption between them."""
    catalogue = tuple(definitions or load_consequence_definitions())
    layers = _subsumption_layers(catalogue)
    known_terms = {definition.term for definition in catalogue}
    nodes: list[Node] = []
    edges: list[Edge] = []

    for definition in catalogue:
        impact = definition.impact.value if definition.impact is not None else ""
        state = f"impact_{impact.lower()}" if impact else ""
        detail = f"rank {definition.rank} · tier {definition.tier}"
        if impact:
            detail = f"{detail} · {impact}"
        nodes.append(
            Node(
                identifier=definition.term,
                label=definition.term,
                kind=TERM_KIND,
                layer=layers[definition.term],
                state=state,
                detail=detail,
            )
        )
        for suppressed in definition.subsumes:
            if suppressed in known_terms:
                edges.append(Edge(source=definition.term, target=suppressed, kind=SUBSUMES_EDGE))

    return Graph(nodes=tuple(nodes), edges=tuple(edges), title="Catalogue overview")


def rule_graphs(ruleset: FactRuleSet | None = None) -> tuple[Graph, ...]:
    """Return one tree per rule, in catalogue rank order.

    One graph each rather than a single merged diagram: the rules share atoms but not structure, so
    merging them produces a tangle that shows neither formula clearly.
    """
    rules = ruleset or build_fact_ruleset()
    graphs = []
    for rule in rules.rules:
        graphs.append(rule_graph(rule.term, rules))
    return tuple(graphs)


def derivation_graphs(trace: AnnotationTrace, ruleset: FactRuleSet | None = None) -> tuple[Graph, ...]:
    """Return every rule decorated with its result for one overlap.

    Each term keeps its own expression tree. Observed True and False expression nodes carry those
    states for the renderer; nodes whose atoms were never observed remain neutral. The graphs are
    grouped by outcome -- emitted first, then firing terms removed by subsumption, then rules that
    did not fire -- and retain rank order within each group.

    Args:
        trace (AnnotationTrace): The derivation to draw.
        ruleset (FactRuleSet | None): The rules the trace was produced against, or None for the
            packaged catalogue's.

    Returns:
        tuple[Graph, ...]: One decorated expression tree per rule.
    """
    rules = ruleset or build_fact_ruleset()
    emitted = set(trace.emitted_terms)
    events = {event.term: event for event in trace.rule_events}
    known = {event.name: bool(event.value) for event in trace.fact_events}
    pruned_by: dict[str, list[str]] = {}
    for prune_event in trace.prune_events:
        pruned_by.setdefault(prune_event.dropped_term, []).append(prune_event.subsumed_by)

    def outcome(rule: FactRule) -> int:
        event = events.get(rule.term)
        if rule.term in emitted:
            return 0
        if event is not None and event.fired:
            return 1
        return 2

    graphs = []
    for rule in sorted(rules.rules, key=lambda rule: (outcome(rule), rule.rank)):
        event = events.get(rule.term)
        if rule.term in emitted:
            term_state = "emitted"
            term_detail = "emitted"
        elif event is not None and event.fired:
            term_state = "pruned"
            suppressors = ", ".join(pruned_by.get(rule.term, ()))
            term_detail = f"pruned by {suppressors}" if suppressors else "pruned"
        else:
            term_state = ""
            term_detail = "not fired" if event is not None else "not evaluated"

        graphs.append(
            _rule_graph(
                rule,
                known=known,
                term_state=term_state,
                term_detail=term_detail,
            )
        )
    return tuple(graphs)
