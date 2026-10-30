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
"""Compile and evaluate rank-ordered consequence rules.

Rules do not reference other rules. Subsumption expresses suppression between firing terms, while
``FactRule.cheap_guard`` derives inexpensive preconditions directly from each formula. Keeping both
relationships declarative prevents hand-maintained predicates or guard tables from drifting.
"""

from __future__ import annotations

from collections.abc import Iterable
from functools import cache

from attrs import field, frozen

from ensembl.vep.annotators.consequence.rules.catalogue import (
    DEFAULT_TERM,
    ConsequenceDefinition,
    load_consequence_definitions,
)
from ensembl.vep.annotators.consequence.rules.expressions import And, Atom, Expr, expr_from_json
from ensembl.vep.annotators.consequence.facts.base import closure_cost
from ensembl.vep.common import ConsequenceImpact, CostTier
from ensembl.vep.models import Feature
from ensembl.vep.models.overlap import OverlapConsequence
from ensembl.vep.models.variation import VariationAllele

_CODING_ATOM = "overlap_cds"
"""The fact a rule must read to count as coding."""

__all__ = [
    "DEFAULT_TERM",
    "FactRule",
    "FactRuleSet",
    "build_fact_ruleset",
    "build_fact_ruleset_from_definitions",
    "coding_terms",
    "prune",
    "subsumes_topological_order",
]


@frozen(kw_only=True, auto_attribs=True)
class FactRule:
    """Represent one SO term as a Boolean formula and catalogue metadata.

    Attributes:
        term (str): The SO term, e.g. "missense_variant".
        accession (str): The SO accession, e.g. "SO:0001583".
        rank (int): Severity rank; 1 is most severe.
        tier (int): Tier grouping, carried through from the catalogue.
        impact (ConsequenceImpact | None): The IMPACT rating.
        when (Expr): The boolean formula over fact atoms that must hold for this term to apply.
        subsumes (tuple[str, ...]): Term names this rule suppresses when evaluated as true.
        feature_type (type[Feature]): The overlapped feature type this rule applies to.
        overlapping_feature_type (type[Feature]): The overlapping feature type (variant allele).
        notes (str | None): Free-text notes, if any.
    """

    term: str
    accession: str
    rank: int
    tier: int = field(default=3)
    impact: ConsequenceImpact | None
    when: Expr
    subsumes: tuple[str, ...] = field(default=())
    feature_type: type[Feature] = field(default=Feature)
    overlapping_feature_type: type[Feature] = field(default=VariationAllele)
    notes: str | None = field(default=None)

    def atoms(self) -> set[str]:
        """Return every fact atom referenced by this rule.

        Returns:
            set[str]: The referenced fact atom names.
        """
        return self.when.atoms()

    def cheap_guard(self) -> frozenset[str]:
        """Return inexpensive top-level positive conjuncts.

        This is the automatically-derived equivalent of VEP's hand-maintained ``include`` hash.

        Returns:
            frozenset[str]: Names of top-level atoms in a top-level conjunction whose
                closure cost is at most BOUNDS.
        """
        conjuncts: tuple[Expr, ...] = self.when.children if isinstance(self.when, And) else (self.when,)
        cheap_atoms = set()
        for conjunct in conjuncts:
            if isinstance(conjunct, Atom) and closure_cost(conjunct.name) <= CostTier.BOUNDS:
                cheap_atoms.add(conjunct.name)
        return frozenset(cheap_atoms)

    def to_overlap_consequence(self) -> OverlapConsequence:
        """Convert this rule into an OverlapConsequence output object.

        Returns:
            OverlapConsequence: The corresponding output consequence, with ``is_default=False``.
        """
        return OverlapConsequence(
            term=self.term,
            accession=self.accession,
            is_default=False,
            rank=self.rank,
            tier=self.tier,
            impact=self.impact,
            feature_type=self.feature_type,
            overlapping_feature_type=self.overlapping_feature_type,
        )


@frozen(kw_only=True, auto_attribs=True)
class FactRuleSet:
    """Store rank-ordered fact rules and term lookup.

    Covers **every** overlapped-feature class, and there is no per-class subset to select: a
    rule restricted to one class says so in its own formula, via the FREE class facts in
    consequence/facts/kind.py. See build_fact_ruleset.

    Attributes:
        rules (tuple[FactRule, ...]): Every rule, sorted by ascending rank.
        by_term (dict[str, FactRule]): Rules keyed by term name.
    """

    rules: tuple[FactRule, ...]
    by_term: dict[str, FactRule] = field(init=False)

    def __attrs_post_init__(self) -> None:
        object.__setattr__(self, "by_term", {rule.term: rule for rule in self.rules})


def build_fact_ruleset_from_definitions(
    definitions: Iterable[ConsequenceDefinition],
) -> FactRuleSet:
    """Compile a rule set from catalogue definitions.

    The single compilation path: ``build_fact_ruleset`` reads the packaged catalogue and delegates
    here, and the catalogue editor compiles its candidate entries the same way, so a change cannot
    validate against a rule set built differently from the one the annotators will use.

    Entries with no ``when`` (see EXCLUDED_TERMS in consequence/rules/catalogue.py) are skipped,
    since there is no fact-atom formula to compile for them.

    Args:
        definitions (Iterable[ConsequenceDefinition]): The entries to compile.

    Returns:
        FactRuleSet: The compiled, rank-ordered rule set.
    """
    compilable = [definition for definition in definitions if definition.when is not None]
    rules = tuple(
        FactRule(
            term=definition.term,
            accession=definition.accession,
            rank=definition.rank,
            impact=definition.impact,
            when=expr_from_json(definition.when),
            subsumes=definition.subsumes,
            tier=definition.tier,
            notes=definition.notes,
        )
        for definition in sorted(compilable, key=lambda definition: definition.rank)
    )
    return FactRuleSet(rules=rules)


def build_fact_ruleset() -> FactRuleSet:
    """Compile a FactRuleSet from the packaged consequence catalogue.

    Returns:
        FactRuleSet: The compiled, rank-ordered rule set, covering every feature class.
    """
    return build_fact_ruleset_from_definitions(load_consequence_definitions())


@cache
def coding_terms() -> frozenset[str]:
    """Return the catalogue terms whose rule reads ``overlap_cds``.

    Derived rather than tabulated: a rule that reads the CDS is a coding rule, and that is the
    whole definition.
    """
    return frozenset(rule.term for rule in build_fact_ruleset().rules if _CODING_ATOM in rule.when.atoms())


def subsumes_topological_order(ruleset: FactRuleSet) -> tuple[FactRule, ...]:
    """Order the rules so a rule comes after a rule that can suppress it.

    After every rule has been evaluated, the subsumed terms are dropped. A rule whose term is
    certain to be suppressed by a higher-ranked term therefore evaluated for nothing. Evaluating
    in this order instead lets a scan skip such a rule outright.

    The ``subsumes`` graph is acyclic. The ``RulePairFinding.contradictory_in`` and
    ``find_subsumption_cycles`` in rules/validation.py both reject cycles. So, a depth-first
    post-order is a valid topological sort. A cycle, were one ever introduced, it would not
    degrade to an infinite loop, and the conformance suite would catch the resulting output difference.

    Args:
        ruleset (FactRuleSet): The rule set to order.

    Returns:
        tuple[FactRule, ...]: Every rule in ``ruleset``, ordered so that a rule appears after
            every rule that (transitively) subsumes it.
    """
    rules_by_rank = sorted(ruleset.rules, key=lambda rule: rule.rank)
    suppressors: dict[str, list[str]] = {rule.term: [] for rule in rules_by_rank}
    for rule in rules_by_rank:
        for subsumed_term in rule.subsumes:
            if subsumed_term in suppressors:
                suppressors[subsumed_term].append(rule.term)

    ordered: list[FactRule] = []
    visit_state: dict[str, int] = {}  # 0/absent = unvisited, 1 = in progress, 2 = emitted

    def visit(term: str) -> None:
        if visit_state.get(term, 0) != 0:
            return
        visit_state[term] = 1
        for suppressor_term in suppressors[term]:
            visit(suppressor_term)
        visit_state[term] = 2
        ordered.append(ruleset.by_term[term])

    for rule in rules_by_rank:
        visit(rule.term)
    return tuple(ordered)


def prune(ruleset: FactRuleSet, hit_terms: list[str], *, apply_default: bool = True) -> list[str]:
    """Drop terms superseded by another term that fired.

    Args:
        ruleset (FactRuleSet): The rule set ``hit_terms`` were matched against.
        hit_terms (list[str]): Every term whose rule matched, in any order.
        apply_default (bool): If True, substitute ``[DEFAULT_TERM]`` when nothing survives pruning.
            If False, return an empty list instead and let the caller decide whether/how to apply a
            default term.

    Returns:
        list[str]: ``hit_terms`` with every subsumed term removed. Empty (if ``apply_default`` is
            False) or ``[DEFAULT_TERM]`` (otherwise) when nothing survives.
    """
    dropped: set[str] = set()
    for term in hit_terms:
        dropped.update(ruleset.by_term[term].subsumes)
    kept = [term for term in hit_terms if term not in dropped]
    if kept or not apply_default:
        return kept
    return [DEFAULT_TERM]
