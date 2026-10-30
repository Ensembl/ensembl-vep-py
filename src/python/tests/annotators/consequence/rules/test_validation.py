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

# pylint: disable=redefined-outer-name
"""Tests for consequence rule validation."""

from __future__ import annotations

import sys

import attrs
import pytest

from ensembl.vep.annotators.consequence.facts.base import verify_requires
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.rules import validation as validation_module
from ensembl.vep.annotators.consequence.rules.catalogue import EXCLUDED_TERMS, ConsequenceDefinition
from ensembl.vep.annotators.consequence.rules.expressions import And, Atom, Expr, Not
from ensembl.vep.annotators.consequence.rules.rulesets import FactRule, FactRuleSet
from ensembl.vep.annotators.consequence.rules.validation import (
    FactAssignmentSpace,
    RulePairFinding,
    RuleReachability,
    SolverChecks,
    StructuralChecks,
    ValidationReport,
)
from ensembl.vep.corpus.synthetic import SyntheticCorpus


class _UnsupportedExpr(Expr):
    """An expression node the z3 translator does not handle."""

    def atoms(self) -> set[str]:
        return set()

    def evaluate(self, fact_context: FactContext) -> bool:
        return False

    def evaluate_partial(self, known: dict[str, bool]) -> bool | None:
        return None


class TestTier1StructuralChecks:
    """Structural checks that do not evaluate formulas."""

    def test_real_ruleset_has_no_structural_errors(self, transcript_ruleset: FactRuleSet) -> None:
        """The real ruleset has no structural errors."""
        errors, _ = StructuralChecks.ruleset_problems(transcript_ruleset)
        assert not errors

    def test_real_ruleset_has_no_subsumption_cycles(self, transcript_ruleset: FactRuleSet) -> None:
        """The real ruleset has no subsumption cycles."""
        assert not StructuralChecks.subsumption_cycles(transcript_ruleset)

    def test_two_cycle_is_detected(self) -> None:
        """A direct two-rule cycle is detected."""
        rule_a = FactRule(
            term="a", accession="SO:1", rank=1, impact=None, when=Atom("is_snv"), subsumes=("b",)
        )
        rule_b = FactRule(
            term="b", accession="SO:2", rank=2, impact=None, when=Atom("is_snv"), subsumes=("a",)
        )
        cycles = StructuralChecks.subsumption_cycles(FactRuleSet(rules=(rule_a, rule_b)))
        assert len(cycles) == 1
        assert "a" in cycles[0] and "b" in cycles[0]

    def test_three_cycle_is_detected(self) -> None:
        """A longer cycle is still caught by the DFS cycle detector."""
        rule_a = FactRule(
            term="a", accession="SO:1", rank=1, impact=None, when=Atom("is_snv"), subsumes=("b",)
        )
        rule_b = FactRule(
            term="b", accession="SO:2", rank=2, impact=None, when=Atom("is_snv"), subsumes=("c",)
        )
        rule_c = FactRule(
            term="c", accession="SO:3", rank=3, impact=None, when=Atom("is_snv"), subsumes=("a",)
        )
        ruleset = FactRuleSet(rules=(rule_a, rule_b, rule_c))
        errors, _ = StructuralChecks.ruleset_problems(ruleset)
        assert not any("mutual suppression" in error for error in errors)
        cycles = StructuralChecks.subsumption_cycles(ruleset)
        assert len(cycles) == 1
        assert all(term in cycles[0] for term in ("a", "b", "c"))

    def test_no_cycles_when_subsumes_graph_is_acyclic(self) -> None:
        """An acyclic subsumes graph has no cycles."""
        rule_a = FactRule(
            term="a", accession="SO:1", rank=1, impact=None, when=Atom("is_snv"), subsumes=("b",)
        )
        rule_b = FactRule(term="b", accession="SO:2", rank=2, impact=None, when=Atom("is_snv"))
        assert not StructuralChecks.subsumption_cycles(FactRuleSet(rules=(rule_a, rule_b)))

    def test_catalogue_when_presence_matches_excluded_terms_exactly(self) -> None:
        """``missing_when_terms`` matches excluded terms exactly on current catalogue."""
        assert not StructuralChecks.missing_when_terms()

    def test_so_parent_drift_report_runs_and_returns_strings(self, transcript_ruleset: FactRuleSet) -> None:
        """SO-parent drift reporting runs and returns strings."""
        drift = StructuralChecks.subsumes_so_parent_drift(transcript_ruleset)
        assert isinstance(drift, list)
        assert all(isinstance(entry, str) for entry in drift)

    def test_a_subsumes_target_outside_the_ruleset_is_skipped_by_the_cycle_walk(self) -> None:
        """An unknown target is validate_ruleset's error to report."""
        rule = FactRule(
            term="a", accession="SO:1", rank=1, impact=None, when=Atom("is_snv"), subsumes=("nowhere",)
        )
        assert not StructuralChecks.subsumption_cycles(FactRuleSet(rules=(rule,)))

    def test_one_cycle_reachable_from_two_roots_is_reported_once(self) -> None:
        """A cycle reachable from two roots is reported once."""
        entry_point = FactRule(
            term="entry", accession="SO:0", rank=1, impact=None, when=Atom("is_snv"), subsumes=("a",)
        )
        rule_a = FactRule(
            term="a", accession="SO:1", rank=2, impact=None, when=Atom("is_snv"), subsumes=("b",)
        )
        rule_b = FactRule(
            term="b", accession="SO:2", rank=3, impact=None, when=Atom("is_snv"), subsumes=("a",)
        )
        cycles = StructuralChecks.subsumption_cycles(FactRuleSet(rules=(entry_point, rule_a, rule_b)))
        assert len(cycles) == 1

    def test_missing_when_report_flags_an_excluded_term_that_grew_a_when(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Missing-when reporting flags drift in either direction."""
        definitions = [
            ConsequenceDefinition(
                term=min(EXCLUDED_TERMS),
                accession="SO:1",
                rank=1,
                tier=1,
                impact=None,
                when={"all": ["is_snv"]},
            ),
            ConsequenceDefinition(
                term="not_excluded_and_has_no_when",
                accession="SO:2",
                rank=2,
                tier=1,
                impact=None,
                when=None,
            ),
        ]
        monkeypatch.setattr(validation_module, "load_consequence_definitions", lambda: definitions)

        problems = StructuralChecks.missing_when_terms()

        assert any("listed in EXCLUDED_TERMS" in problem for problem in problems)
        assert any("not listed in EXCLUDED_TERMS" in problem for problem in problems)

    def test_so_parent_drift_skips_a_rule_the_catalogue_does_not_define(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A rule with no catalogue entry has nothing to compare."""
        monkeypatch.setattr(validation_module, "load_consequence_definitions", list)
        rule = FactRule(
            term="not_in_the_catalogue", accession="SO:1", rank=1, impact=None, when=Atom("is_snv")
        )
        assert not StructuralChecks.subsumes_so_parent_drift(FactRuleSet(rules=(rule,)))


class TestTier1ValidateRuleset:
    """Ruleset validation errors and warnings."""

    def test_ruleset_problems_reports_the_expected_rank_inversion_warnings(
        self,
        transcript_ruleset: FactRuleSet,
    ) -> None:
        """The current catalogue's subsumes data yields the expected warnings."""
        _, warnings = StructuralChecks.ruleset_problems(transcript_ruleset)
        expected_pairs = {
            ("stop_lost", "stop_gained"),
            ("splice_donor_region_variant", "splice_region_variant"),
            ("incomplete_terminal_codon_variant", "frameshift_variant"),
            ("incomplete_terminal_codon_variant", "inframe_insertion"),
            ("incomplete_terminal_codon_variant", "inframe_deletion"),
            ("stop_retained_variant", "frameshift_variant"),
            ("stop_retained_variant", "missense_variant"),
            ("stop_retained_variant", "stop_gained"),
        }
        assert len(warnings) == len(expected_pairs)
        for suppressor, suppressed in expected_pairs:
            assert any(suppressor in warning and suppressed in warning for warning in warnings), (
                suppressor,
                suppressed,
            )

    def test_ruleset_problems_catches_an_unknown_fact_reference(self) -> None:
        """An unknown fact reference is rejected."""
        broken_rule = FactRule(
            term="broken_term",
            accession="SO:9999999",
            rank=1,
            impact=None,
            when=Atom("this_fact_does_not_exist"),
        )
        broken_ruleset = FactRuleSet(rules=(broken_rule,))
        errors, _warnings = StructuralChecks.ruleset_problems(broken_ruleset)
        assert any("unknown fact" in error for error in errors)
        assert {field.name for field in attrs.fields(FactRule)} >= {"term", "when", "subsumes"}

    def test_ruleset_problems_catches_mutual_suppression(self) -> None:
        """Mutual suppression is rejected."""
        rule_a = FactRule(
            term="a", accession="SO:1", rank=1, impact=None, when=Atom("is_snv"), subsumes=("b",)
        )
        rule_b = FactRule(
            term="b", accession="SO:2", rank=2, impact=None, when=Atom("is_snv"), subsumes=("a",)
        )
        ruleset = FactRuleSet(rules=(rule_a, rule_b))
        errors, _warnings = StructuralChecks.ruleset_problems(ruleset)
        assert any("mutual suppression" in error for error in errors)

    def test_ruleset_problems_catches_duplicate_rank(self) -> None:
        """Duplicate ranks are rejected."""
        rule_a = FactRule(term="a", accession="SO:1", rank=1, impact=None, when=Atom("is_snv"))
        rule_b = FactRule(term="b", accession="SO:2", rank=1, impact=None, when=Atom("is_snv"))
        ruleset = FactRuleSet(rules=(rule_a, rule_b))
        errors, _warnings = StructuralChecks.ruleset_problems(ruleset)
        assert any("duplicate rank" in error for error in errors)

    def test_ruleset_problems_catches_a_reference_to_an_internal_fact(self) -> None:
        """A rule may only reference public atoms."""
        rule = FactRule(term="a", accession="SO:1", rank=1, impact=None, when=Atom("_codon_alleles"))
        errors, _warnings = StructuralChecks.ruleset_problems(FactRuleSet(rules=(rule,)))
        assert any("internal fact" in error for error in errors)

    def test_ruleset_problems_catches_a_subsumes_target_outside_the_ruleset(self) -> None:
        """A subsumes target outside the ruleset is rejected."""
        rule = FactRule(
            term="a", accession="SO:1", rank=1, impact=None, when=Atom("is_snv"), subsumes=("nowhere",)
        )
        errors, _warnings = StructuralChecks.ruleset_problems(FactRuleSet(rules=(rule,)))
        assert any("subsumes unknown term" in error for error in errors)


class TestTier2AtomSupport:
    """Atom support underpins Tier 2 enumeration."""

    def test_support_is_a_superset_of_the_rules_own_atoms(self, transcript_ruleset: FactRuleSet) -> None:
        """Support is a superset of the rule's own atoms."""
        for rule in transcript_ruleset.rules:
            assert set(rule.atoms()) <= set(FactAssignmentSpace.for_rule(rule).support)

    def test_requires_closure_pulls_in_transitive_dependencies(self, transcript_ruleset: FactRuleSet) -> None:
        """Transitive dependencies appear in the support."""
        rule = transcript_ruleset.by_term["stop_gained"]
        support = set(FactAssignmentSpace.for_rule(rule).support)
        assert "overlap_exon" in support
        assert "overlap_feature" in support

    def test_axiom_partner_atoms_are_pulled_in(self, transcript_ruleset: FactRuleSet) -> None:
        """Axiom partner atoms are included in the support."""
        rule = transcript_ruleset.by_term["inframe_insertion"]
        support = set(FactAssignmentSpace.for_rule(rule).support)
        assert "is_insertion" in support
        assert "is_deletion" in support

    def test_an_atom_outside_the_registry_contributes_no_closure(self) -> None:
        """An unknown atom contributes no additional closure."""
        rule = FactRule(term="a", accession="SO:1", rank=1, impact=None, when=Atom("not_a_registered_fact"))
        assert FactAssignmentSpace.for_rule(rule).support == ("not_a_registered_fact",)
        # It also declares no `requires`, so neither assignment over it can be inconsistent.
        assert RuleReachability.for_rule(rule).valid_assignment_count == 2

    def test_an_equivalence_axiom_pair_constrains_the_enumeration(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Declared-equivalent atoms may not disagree."""
        rule = FactRule(
            term="a",
            accession="SO:1",
            rank=1,
            impact=None,
            when=And(Atom("is_snv"), Not(child=Atom("is_insertion"))),
        )
        assert RuleReachability.for_rule(rule).is_reachable

        monkeypatch.setattr(validation_module, "_EQUIVALENCE_AXIOM_PAIRS", (("is_snv", "is_insertion"),))

        assert not RuleReachability.for_rule(rule).is_reachable


class TestTier2Reachability:
    """Per-rule brute-force reachability enumeration."""

    def test_evaluate_partial_never_returns_none_over_a_rules_own_support(
        self, transcript_ruleset: FactRuleSet
    ) -> None:
        """The support is constructed to fully settle every rule's formula (no Kleene unknowns)."""
        for rule in transcript_ruleset.rules:
            result = RuleReachability.for_rule(rule)
            assert result.valid_assignment_count > 0

    def test_missense_variant_is_reachable(self, transcript_ruleset: FactRuleSet) -> None:
        """``missense_variant`` is reachable."""
        result = RuleReachability.for_rule(transcript_ruleset.by_term["missense_variant"])
        assert result.is_reachable

    def test_a_formula_that_is_its_own_negation_is_unreachable(self) -> None:
        """A trivially self-contradictory rule is unreachable."""
        rule = FactRule(
            term="impossible",
            accession="SO:0",
            rank=1,
            impact=None,
            when=Atom("is_snv") & ~Atom("is_snv"),
        )
        result = RuleReachability.for_rule(rule)
        assert not result.is_reachable
        assert result.satisfying_assignment_count == 0

    def test_unreachable_in_returns_reachability_objects(self, transcript_ruleset: FactRuleSet) -> None:
        """Every unreachable entry is well-typed."""
        for result in RuleReachability.unreachable_in(transcript_ruleset):
            assert isinstance(result, RuleReachability)
            assert not result.is_reachable


class TestTier2PairwiseFindings:
    """Pairwise contradictory and redundant findings."""

    def test_contradictory_pairs_never_include_a_subsumption_edge(
        self, transcript_ruleset: FactRuleSet
    ) -> None:
        """Contradictory pairs never include a subsumption edge."""
        for finding in RulePairFinding.contradictory_in(transcript_ruleset):
            rule_a = transcript_ruleset.by_term[finding.term_a]
            rule_b = transcript_ruleset.by_term[finding.term_b]
            assert finding.term_b not in rule_a.subsumes
            assert finding.term_a not in rule_b.subsumes

    def test_coding_floor_term_subsumption_is_confirmed_redundant(
        self, transcript_ruleset: FactRuleSet
    ) -> None:
        """The stop_gained to coding_sequence_variant edge is provably entailed."""
        findings = RulePairFinding.redundant_in(transcript_ruleset)
        redundant_edges = {(finding.term_a, finding.term_b) for finding in findings}
        assert ("stop_gained", "coding_sequence_variant") in redundant_edges

    def test_redundant_pairs_are_always_declared_subsumption_edges(
        self, transcript_ruleset: FactRuleSet
    ) -> None:
        """Redundant pairs are declared subsumption edges."""
        for finding in RulePairFinding.redundant_in(transcript_ruleset):
            rule_a = transcript_ruleset.by_term[finding.term_a]
            assert finding.term_b in rule_a.subsumes

    def test_transitive_requirements_are_propagated_when_finding_contradictions(self) -> None:
        """Transitive requirements can make two rule expressions contradictory."""
        utr_rule = FactRule(
            term="utr",
            accession="SO:1",
            rank=1,
            impact=None,
            when=Atom("overlap_five_prime_utr"),
        )
        upstream_rule = FactRule(
            term="upstream",
            accession="SO:2",
            rank=2,
            impact=None,
            when=Atom("is_upstream"),
        )

        findings = RulePairFinding.contradictory_in(FactRuleSet(rules=(utr_rule, upstream_rule)))

        assert [(finding.term_a, finding.term_b) for finding in findings] == [("utr", "upstream")]

    def test_transitive_requirements_are_propagated_when_finding_redundancy(self) -> None:
        """Transitive requirements can prove a declared subsumption edge."""
        utr_rule = FactRule(
            term="utr",
            accession="SO:1",
            rank=1,
            impact=None,
            when=Atom("overlap_five_prime_utr"),
            subsumes=("feature",),
        )
        feature_rule = FactRule(
            term="feature",
            accession="SO:2",
            rank=2,
            impact=None,
            when=Atom("overlap_feature"),
        )

        findings = RulePairFinding.redundant_in(FactRuleSet(rules=(utr_rule, feature_rule)))

        assert [(finding.term_a, finding.term_b) for finding in findings] == [("utr", "feature")]

    def test_a_subsumes_target_outside_the_ruleset_yields_no_finding(self) -> None:
        """A subsumes target outside the ruleset yields no finding."""
        rule = FactRule(
            term="a", accession="SO:1", rank=1, impact=None, when=Atom("is_snv"), subsumes=("nowhere",)
        )
        assert not RulePairFinding.redundant_in(FactRuleSet(rules=(rule,)))

    def test_an_edge_whose_union_support_is_too_large_is_skipped(
        self, transcript_ruleset: FactRuleSet
    ) -> None:
        """A zero support cap skips every edge."""
        assert not RulePairFinding.redundant_in(transcript_ruleset, max_union_support=0)


class TestValidationReport:
    """The combined Tier 1 and Tier 2 report."""

    def test_report_runs_over_the_real_ruleset_and_is_well_typed(
        self, transcript_ruleset: FactRuleSet
    ) -> None:
        """The report runs over the real ruleset."""
        report = ValidationReport.for_ruleset(transcript_ruleset)
        assert not report.errors
        assert not report.missing_when_problems
        # has_errors reflects errors/missing_when_problems/unreachable_rules only -- contradictory
        # and redundant findings must never affect it, since they are informational by design.
        assert report.has_errors == bool(report.unreachable_rules)


class TestTier3Z3Validation:
    """Whole-ruleset SAT queries, posed to z3 rather than enumerated.

    Tier 3 answers the same two questions as Tier 2 -- is a rule ever satisfiable, and does a
    declared ``subsumes`` edge hold as entailment -- without being bounded by one rule's own atom
    support, so the two tiers agreeing on the real catalogue is the check worth making.
    """

    def test_the_real_ruleset_reports_no_unreachable_rules(self, transcript_ruleset: FactRuleSet) -> None:
        """The shipped catalogue has no unreachable rules."""
        findings = SolverChecks.findings(transcript_ruleset)
        assert not [finding for finding in findings if finding.startswith("UNREACHABLE")]

    def test_z3_and_enumeration_agree_on_which_rules_are_unreachable(
        self, transcript_ruleset: FactRuleSet
    ) -> None:
        """The two tiers agree on unreachable rules."""
        z3_unreachable = {
            finding.split()[2]
            for finding in SolverChecks.findings(transcript_ruleset)
            if finding.startswith("UNREACHABLE")
        }
        enumerated_unreachable = {
            result.term for result in RuleReachability.unreachable_in(transcript_ruleset)
        }
        assert z3_unreachable == enumerated_unreachable

    def test_not_entailed_findings_are_informational_and_name_a_declared_edge(
        self, transcript_ruleset: FactRuleSet
    ) -> None:
        """Every NOT ENTAILED finding names a declared edge."""
        declared_edges = {
            (rule.term, subsumed_term)
            for rule in transcript_ruleset.rules
            for subsumed_term in rule.subsumes
            if subsumed_term in transcript_ruleset.by_term
        }
        for finding in SolverChecks.findings(transcript_ruleset):
            if not finding.startswith("NOT ENTAILED"):
                continue
            words = finding.split()
            assert (words[3], words[5]) in declared_edges, finding

    def test_an_unsatisfiable_rule_is_reported_unreachable(self) -> None:
        """An unsatisfiable rule is reported unreachable."""
        rule = FactRule(
            term="contradiction",
            accession="SO:1",
            rank=1,
            impact=None,
            when=And(Atom("is_snv"), Not(child=Atom("is_snv"))),
        )
        findings = SolverChecks.findings(FactRuleSet(rules=(rule,)))
        assert findings == ["UNREACHABLE (z3): contradiction is unsatisfiable under the declared axioms"]

    def test_an_unentailed_subsumption_edge_is_reported(self) -> None:
        """An unentailed subsumption edge is reported."""
        suppressor = FactRule(
            term="a", accession="SO:1", rank=1, impact=None, when=Atom("is_snv"), subsumes=("b",)
        )
        suppressed = FactRule(term="b", accession="SO:2", rank=2, impact=None, when=Atom("is_insertion"))
        findings = SolverChecks.findings(FactRuleSet(rules=(suppressor, suppressed)))
        assert any(finding.startswith("NOT ENTAILED (z3): a subsumes b") for finding in findings)

    def test_an_entailed_subsumption_edge_is_not_reported(self) -> None:
        """An entailed subsumption edge is not reported."""
        suppressor = FactRule(
            term="a", accession="SO:1", rank=1, impact=None, when=Atom("is_snv"), subsumes=("b",)
        )
        suppressed = FactRule(term="b", accession="SO:2", rank=2, impact=None, when=Atom("is_snv"))
        assert not SolverChecks.findings(FactRuleSet(rules=(suppressor, suppressed)))

    def test_a_subsumes_target_outside_the_ruleset_is_skipped(self) -> None:
        """A subsumes target outside the ruleset is skipped."""
        rule = FactRule(
            term="a", accession="SO:1", rank=1, impact=None, when=Atom("is_snv"), subsumes=("nowhere",)
        )
        assert not SolverChecks.findings(FactRuleSet(rules=(rule,)))

    def test_an_unsupported_expression_node_is_rejected(self) -> None:
        """Unsupported expression nodes are rejected."""
        rule = FactRule(term="a", accession="SO:1", rank=1, impact=None, when=_UnsupportedExpr())
        with pytest.raises(TypeError, match="unsupported expression node type"):
            SolverChecks.findings(FactRuleSet(rules=(rule,)))

    def test_an_atom_outside_the_registry_carries_no_requires_axiom(self) -> None:
        """An unknown atom carries no requires axiom."""
        rule = FactRule(term="a", accession="SO:1", rank=1, impact=None, when=Atom("not_a_registered_fact"))
        assert not SolverChecks.findings(FactRuleSet(rules=(rule,)))

    def test_an_equivalence_axiom_pair_is_asserted_to_the_solver(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """An equivalence axiom pair is asserted to the solver."""
        suppressor = FactRule(
            term="a", accession="SO:1", rank=1, impact=None, when=Atom("is_snv"), subsumes=("b",)
        )
        suppressed = FactRule(term="b", accession="SO:2", rank=2, impact=None, when=Atom("is_insertion"))
        ruleset = FactRuleSet(rules=(suppressor, suppressed))
        assert SolverChecks.findings(ruleset)

        monkeypatch.setattr(validation_module, "_EQUIVALENCE_AXIOM_PAIRS", (("is_snv", "is_insertion"),))

        assert not SolverChecks.findings(ruleset)

    def test_a_missing_solver_raises_and_names_the_extra(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Missing z3 raises and names the extra."""
        monkeypatch.setitem(sys.modules, "z3", None)
        with pytest.raises(RuntimeError, match="validation"):
            SolverChecks.findings(FactRuleSet(rules=()))


class TestCatalogueEntryChecks:
    """Checks that run over raw catalogue entries."""

    def test_duplicate_terms_and_ranks_are_both_reported(self) -> None:
        """Duplicate terms and ranks are both reported."""
        entries = [
            {"term": "a", "accession": "SO:1", "rank": 1, "tier": 1},
            {"term": "a", "accession": "SO:2", "rank": 1, "tier": 1},
        ]
        errors = StructuralChecks.term_and_rank_uniqueness(entries)
        assert any("duplicate catalogue term: a" in error for error in errors)
        assert any("duplicate catalogue rank 1" in error for error in errors)

    def test_unique_terms_and_ranks_report_nothing(self) -> None:
        """Unique terms and ranks report nothing."""
        entries = [
            {"term": "a", "accession": "SO:1", "rank": 1, "tier": 1},
            {"term": "b", "accession": "SO:2", "rank": 2, "tier": 1},
        ]
        assert not StructuralChecks.term_and_rank_uniqueness(entries)

    def test_a_pre_existing_unreachable_rule_is_downgraded_to_a_warning(self) -> None:
        """A pre-existing unreachable rule is downgraded to a warning."""
        dead_entry = {
            "term": "contradiction",
            "accession": "SO:1",
            "rank": 1,
            "tier": 1,
            "when": {"all": ["is_snv", {"not": "is_snv"}]},
        }
        report = ValidationReport.for_catalogue_entries([dead_entry], [dead_entry])

        assert not report.unreachable_rules
        assert any("pre-existing unreachable rule" in warning for warning in report.warnings)

    def test_a_newly_unreachable_rule_stays_build_breaking(self) -> None:
        """A newly unreachable rule stays build breaking."""
        live_entry = {
            "term": "contradiction",
            "accession": "SO:1",
            "rank": 1,
            "tier": 1,
            "when": {"all": ["is_snv"]},
        }
        dead_entry = dict(live_entry, when={"all": ["is_snv", {"not": "is_snv"}]})
        report = ValidationReport.for_catalogue_entries([dead_entry], [live_entry])

        assert [result.term for result in report.unreachable_rules] == ["contradiction"]
        assert report.has_errors


class TestVerifyRequiresOverTheSyntheticCorpus:
    """``verify_requires`` over real contexts."""

    def test_every_requires_claim_holds_over_the_boundary_and_fuzz_corpus(
        self, small_default_corpus: SyntheticCorpus
    ) -> None:
        """Every ``requires`` claim holds over the boundary and fuzz corpus."""
        fact_contexts = [
            FactContext(overlap=overlap) for _name, overlap in small_default_corpus.named_overlaps
        ]
        assert not verify_requires(fact_contexts)
