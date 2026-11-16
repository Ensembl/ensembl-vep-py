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
"""Tests for the fact and rule graph builders."""

from __future__ import annotations

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the facts and rules these graphs are built from.
from ensembl.vep.annotators.consequence.facts.base import REGISTRY
from ensembl.vep.annotators.consequence.graphs import (
    DEPENDS_EDGE,
    REQUIRES_EDGE,
    SUBSUMES_EDGE,
    catalogue_overview_graph,
    derivation_graphs,
    fact_dependency_graph,
    rule_graph,
    rule_graphs,
)
from ensembl.vep.annotators.consequence.rules.catalogue import load_consequence_definitions
from ensembl.vep.annotators.consequence.rules.expressions import all_of, no
from ensembl.vep.annotators.consequence.rules.rulesets import FactRule, FactRuleSet, build_fact_ruleset
from ensembl.vep.render.graph import Graph
from ensembl.vep.trace import AnnotationTrace, FactEvent, FactOrigin, PruneEvent, RuleEvent


@pytest.fixture(name="facts")
def fixture_facts() -> Graph:
    """The fact dependency graph."""
    return fact_dependency_graph()


class TestFactGraph:
    """The fact registry as a graph."""

    def test_every_registered_fact_is_a_node(self, facts: Graph) -> None:
        """A missing node would hide a fact from the diagram entirely."""
        assert {node.identifier for node in facts.nodes} == set(REGISTRY)

    def test_facts_are_layered_by_cost_tier(self, facts: Graph) -> None:
        """The domain's own left-to-right ordering: cheapest first.

        Six layers, one per ``CostTier`` member actually used by a registered fact --
        ``_full_cds_peptide_alleles`` is the sole ``TRANSCRIPT_SEQUENCE`` occupant of the last one.
        """
        assert facts.layers() == (0, 1, 2, 3, 4, 5)

    def test_a_free_fact_sits_in_the_first_layer(self, facts: Graph) -> None:
        """FREE is the cheapest tier, so it is the leftmost column."""
        assert facts.node("is_snv").layer == 0

    def test_the_tier_name_is_carried_as_detail(self, facts: Graph) -> None:
        """So a renderer can label the column without knowing about CostTier."""
        assert facts.node("is_snv").detail == "FREE"

    def test_both_edge_kinds_are_distinguished(self, facts: Graph) -> None:
        """``deps`` is what must be computed first; ``requires`` is what can settle a rule early."""
        kinds = {edge.kind for edge in facts.edges}
        assert kinds == {DEPENDS_EDGE, REQUIRES_EDGE}

    def test_every_edge_endpoint_is_a_node(self, facts: Graph) -> None:
        """An edge to an unknown identifier would be silently dropped by the renderer."""
        identifiers = {node.identifier for node in facts.nodes}
        for edge in facts.edges:
            assert edge.source in identifiers
            assert edge.target in identifiers

    def test_nodes_carry_no_styling(self, facts: Graph) -> None:
        """Appearance is the renderer's decision; a builder emits only semantics."""
        for node in facts.nodes:
            assert node.kind == "fact"
            assert not node.state


class TestRuleGraph:
    """One rule as a Boolean expression tree."""

    def test_the_formula_becomes_a_tree(self) -> None:
        """Operators are nodes, so the formula's shape survives into the diagram."""
        graph = rule_graph("splice_donor_5th_base_variant")
        kinds = {node.kind for node in graph.nodes}
        assert "operator" in kinds
        assert "atom" in kinds

    def test_atoms_sit_left_of_the_operators_they_feed(self) -> None:
        """The tree reads left to right, deepest atoms first."""
        graph = rule_graph("splice_donor_5th_base_variant")
        atom_layers = {node.layer for node in graph.nodes if node.kind == "atom"}
        operator_layers = {node.layer for node in graph.nodes if node.kind == "operator"}
        assert max(atom_layers) < min(operator_layers)

    def test_the_term_is_the_rightmost_of_the_formula(self) -> None:
        """Everything in the formula feeds it."""
        graph = rule_graph("splice_donor_5th_base_variant")
        term = graph.node("splice_donor_5th_base_variant")
        formula_layers = {node.layer for node in graph.nodes if node.kind in {"atom", "operator"}}
        assert term.layer > max(formula_layers)

    def test_nesting_is_preserved(self) -> None:
        """A rule with an inner disjunction is deeper than a flat conjunction."""
        flat = rule_graph("splice_donor_5th_base_variant")
        nested = rule_graph("stop_gained")
        assert len(nested.layers()) > len(flat.layers())

    def test_an_atom_used_twice_gets_a_node_each_time(self) -> None:
        """A shared node would pull unrelated branches together and untree the tree."""
        graph = rule_graph("stop_gained")
        labels = [node.label for node in graph.nodes]
        assert sum(label.endswith("peptide_has_X") for label in labels) == 2

    def test_negation_of_an_atom_stays_on_the_leaf(self) -> None:
        """A separate NOT gate would add a column without clarifying the expression."""
        graph = rule_graph("stop_gained")
        labels = {node.label for node in graph.nodes}
        assert "NOT is_structural" in labels
        assert "NOT" not in labels

    def test_the_term_carries_its_catalogue_metadata(self) -> None:
        """The selected result identifies the exact catalogue entry."""
        detail = rule_graph("stop_gained").node("stop_gained").detail
        assert detail == "SO:0001587 · rank 4 · tier 3 · HIGH"

    def test_node_identifiers_are_unique(self) -> None:
        """Duplicated labels are intended; duplicated identifiers would merge them again."""
        graph = rule_graph("stop_gained")
        identifiers = [node.identifier for node in graph.nodes]
        assert len(identifiers) == len(set(identifiers))

    def test_the_formula_root_does_not_collide_with_the_term(self) -> None:
        """They shared an identifier once, which turned the ``when`` edge into a self-edge."""
        graph = rule_graph("splice_donor_5th_base_variant")
        when_edges = [edge for edge in graph.edges if edge.kind == "when"]
        assert len(when_edges) == 1
        assert when_edges[0].source != when_edges[0].target

    def test_an_unknown_term_is_rejected(self) -> None:
        """Named clearly, since the runner turns this into a usage error."""
        with pytest.raises(KeyError, match="no catalogue rule named"):
            rule_graph("nope")

    def test_a_supplied_ruleset_is_used(self) -> None:
        """So a candidate change can be drawn before it is written."""
        ruleset = build_fact_ruleset()
        assert rule_graph("intron_variant", ruleset).nodes

    def test_subsumption_is_not_mixed_into_the_expression(self) -> None:
        """Relations between terms belong to the catalogue overview."""
        graph = rule_graph("splice_region_variant")
        assert SUBSUMES_EDGE not in {edge.kind for edge in graph.edges}
        assert {node.kind for node in graph.nodes} <= {"atom", "operator", "rule"}


class TestCatalogueOverview:
    """The complete catalogue subsumption graph."""

    def test_every_catalogue_term_is_shown(self) -> None:
        """The overview includes the two fallback terms without expressions."""
        definitions = load_consequence_definitions()
        graph = catalogue_overview_graph(definitions)
        assert {node.identifier for node in graph.nodes} == {definition.term for definition in definitions}

    def test_every_declared_subsumption_is_shown(self) -> None:
        """The overview is a rendering of the catalogue relationship, not a sample."""
        definitions = load_consequence_definitions()
        declared = {
            (definition.term, suppressed) for definition in definitions for suppressed in definition.subsumes
        }
        graph = catalogue_overview_graph(definitions)
        drawn = {(edge.source, edge.target) for edge in graph.edges}
        assert drawn == declared

    def test_suppression_always_moves_to_the_right(self) -> None:
        """The overview needs no backward arrows crossing the graph."""
        graph = catalogue_overview_graph()
        for edge in graph.edges:
            assert graph.node(edge.source).layer < graph.node(edge.target).layer

    def test_impact_is_available_to_the_renderer(self) -> None:
        """Overview colours can distinguish consequence impact."""
        assert catalogue_overview_graph().node("stop_gained").state == "impact_high"


class TestRuleGraphs:
    """One diagram per rule, rather than one merged diagram."""

    def test_there_is_one_graph_per_compiled_rule(self) -> None:
        """Merging them would show neither formula clearly."""
        assert len(rule_graphs()) == len(build_fact_ruleset().rules)

    def test_they_come_in_catalogue_rank_order(self) -> None:
        """The order the rules are evaluated and pruned in."""
        titles = [graph.title for graph in rule_graphs()]
        assert titles == [rule.term for rule in build_fact_ruleset().rules]

    def test_each_is_titled_with_its_term(self) -> None:
        """So a reader can tell the diagrams apart when they are emitted back to back."""
        for graph in rule_graphs():
            assert graph.title


@pytest.fixture(name="derivation")
def fixture_derivation() -> tuple[AnnotationTrace, FactRuleSet]:
    """Three rule outcomes over a mix of true, false and unobserved atoms."""
    rules = FactRuleSet(
        rules=(
            FactRule(
                term="blocked_term",
                accession="SO:0000001",
                rank=1,
                impact=None,
                when=all_of("is_transcript", "unobserved_fact"),
            ),
            FactRule(
                term="pruned_term",
                accession="SO:0000002",
                rank=2,
                impact=None,
                when=all_of("overlap_cds"),
            ),
            FactRule(
                term="emitted_term",
                accession="SO:0000003",
                rank=3,
                impact=None,
                when=all_of(no("is_structural"), "overlap_cds"),
            ),
        )
    )
    trace = AnnotationTrace(
        fact_events=(
            FactEvent(
                name="is_transcript",
                value=False,
                cost=0,
                origin=FactOrigin.COMPUTED,
            ),
            FactEvent(
                name="is_structural",
                value=False,
                cost=0,
                origin=FactOrigin.COMPUTED,
            ),
            FactEvent(
                name="overlap_cds",
                value=True,
                cost=0,
                origin=FactOrigin.COMPUTED,
            ),
        ),
        rule_events=(
            RuleEvent(
                term="blocked_term",
                rank=1,
                fired=False,
                settled_early=True,
                reason="settled early",
            ),
            RuleEvent(
                term="pruned_term",
                rank=2,
                fired=True,
                settled_early=False,
                reason="evaluated",
            ),
            RuleEvent(
                term="emitted_term",
                rank=3,
                fired=True,
                settled_early=False,
                reason="evaluated",
            ),
        ),
        prune_events=(PruneEvent(dropped_term="pruned_term", subsumed_by="emitted_term"),),
        emitted_terms=("emitted_term",),
    )
    return trace, rules


class TestDerivationGraphs:
    """Every rule expression decorated with one overlap's result."""

    def test_terms_are_grouped_by_outcome(self, derivation: tuple[AnnotationTrace, FactRuleSet]) -> None:
        """Outcome is more useful than rank when reviewing one completed annotation."""
        trace, rules = derivation
        graphs = derivation_graphs(trace, rules)
        assert [graph.title for graph in graphs] == ["emitted_term", "pruned_term", "blocked_term"]

    def test_conditions_use_the_rule_syntax(self, derivation: tuple[AnnotationTrace, FactRuleSet]) -> None:
        """Truth belongs to colour, leaving labels free of values and status symbols."""
        trace, rules = derivation
        graph = derivation_graphs(trace, rules)[0]
        labels = {node.label for node in graph.nodes}
        assert "NOT is_structural" in labels
        assert all("=" not in label and "✓" not in label for label in labels)

    def test_expression_nodes_carry_their_observed_truth(
        self, derivation: tuple[AnnotationTrace, FactRuleSet]
    ) -> None:
        """Renderers can consistently map true to blue and false to red."""
        trace, rules = derivation
        graphs = {graph.title: graph for graph in derivation_graphs(trace, rules)}
        emitted_states = {node.label: node.state for node in graphs["emitted_term"].nodes}
        blocked_states = {node.label: node.state for node in graphs["blocked_term"].nodes}
        assert emitted_states["NOT is_structural"] == "true"
        assert emitted_states["overlap_cds"] == "true"
        assert emitted_states["ALL"] == "true"
        assert blocked_states["is_transcript"] == "false"
        assert blocked_states["ALL"] == "false"

    def test_an_unobserved_condition_stays_neutral(
        self, derivation: tuple[AnnotationTrace, FactRuleSet]
    ) -> None:
        """A value the trace never saw must not be invented merely to complete the tree."""
        trace, rules = derivation
        blocked = derivation_graphs(trace, rules)[2]
        assert next(node for node in blocked.nodes if node.label == "unobserved_fact").state == ""

    def test_term_nodes_distinguish_each_outcome(
        self, derivation: tuple[AnnotationTrace, FactRuleSet]
    ) -> None:
        """A non-firing term remains neutral while emitted and pruned terms retain their status."""
        trace, rules = derivation
        graphs = {graph.title: graph for graph in derivation_graphs(trace, rules)}
        assert graphs["emitted_term"].node("emitted_term").state == "emitted"
        assert graphs["pruned_term"].node("pruned_term").state == "pruned"
        assert graphs["pruned_term"].node("pruned_term").detail == "pruned by emitted_term"
        assert graphs["blocked_term"].node("blocked_term").state == ""
