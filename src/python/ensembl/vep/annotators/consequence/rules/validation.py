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
"""Validate consequence rule structure and satisfiability.

Validation has three tiers:

* ``StructuralChecks`` validates references, ranks, subsumption, and catalogue metadata.
* ``RuleReachability`` and ``RulePairFinding`` search bounded ``FactAssignmentSpace`` instances.
* Optional ``SolverChecks`` asks equivalent whole-rule-set questions through an SAT solver.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import Any

from attrs import field, frozen

from ensembl.vep.annotators.consequence.rules.catalogue import (
    EXCLUDED_TERMS,
    ConsequenceDefinition,
    load_consequence_definitions,
)
from ensembl.vep.annotators.consequence.facts.base import REGISTRY
from ensembl.vep.annotators.consequence.rules.expressions import And, Atom, Expr, Not, Or
from ensembl.vep.annotators.consequence.rules.rulesets import (
    FactRule,
    FactRuleSet,
    build_fact_ruleset_from_definitions,
)
from ensembl.vep.common import ConsequenceImpact

__all__ = [
    "FactAssignmentSpace",
    "RulePairFinding",
    "RuleReachability",
    "SolverChecks",
    "StructuralChecks",
    "ValidationReport",
]

# --------------------------------------------------------------------------
# Tier 1: structural checks (no formula evaluation).
# --------------------------------------------------------------------------


class StructuralChecks:
    """Run validation checks that require no formula evaluation.

    Grouped rather than left as loose functions because they share one subject - the declared
    shape of a rule set.
    """

    @classmethod
    def ruleset_problems(cls, ruleset: FactRuleSet) -> tuple[list[str], list[str]]:
        """Run structural checks over a rule set.

        Rank inversions are warnings rather than errors on purpose: VEP has real ones
        (``splice_donor_region_variant`` is ranked below ``splice_region_variant`` yet suppresses
        it), so ranking order and suppression order are genuinely two different orders. Making that
        visible is more useful than pretending one implies the other.

        Args:
            ruleset (FactRuleSet): The rule set to validate.

        Returns:
            tuple[list[str], list[str]]: ``(errors, warnings)``. A non-empty ``errors`` list means
                the rule set is broken (unknown/internal fact reference, duplicate rank, mutual/cyclic
                suppression); a non-empty ``warnings`` list flags rank inversions, which are expected
                and not a failure.
        """
        errors: list[str] = []
        warnings: list[str] = []
        seen_ranks: dict[int, str] = {}

        for rule in ruleset.rules:
            if rule.rank in seen_ranks:
                errors.append(f"duplicate rank {rule.rank}: {rule.term} and {seen_ranks[rule.rank]}")
            seen_ranks[rule.rank] = rule.term

            for atom_name in rule.atoms():
                if atom_name not in REGISTRY:
                    errors.append(f"{rule.term} references unknown fact {atom_name!r}")
                elif atom_name.startswith("_"):
                    errors.append(f"{rule.term} references internal fact {atom_name!r}")

            for subsumed_term in rule.subsumes:
                if subsumed_term not in ruleset.by_term:
                    errors.append(f"{rule.term} subsumes unknown term {subsumed_term!r}")
                elif ruleset.by_term[subsumed_term].rank < rule.rank:
                    warnings.append(
                        f"rank inversion: {rule.term} (rank {rule.rank}) suppresses higher-ranked "
                        f"{subsumed_term} (rank {ruleset.by_term[subsumed_term].rank})"
                    )

        # A term that suppresses something which suppresses it back could never settle. This is the
        # 2-cycle special case of subsumption_cycles below; kept as its own check (rather than
        # folded away) so the "mutual suppression" wording -- and the tests that assert on it -- keep
        # working, while subsumption_cycles additionally catches longer cycles.
        for rule in ruleset.rules:
            for subsumed_term in rule.subsumes:
                subsumed_rule = ruleset.by_term.get(subsumed_term)
                if subsumed_rule is not None and rule.term in subsumed_rule.subsumes:
                    errors.append(f"mutual suppression between {rule.term} and {subsumed_term}")

        errors.extend(cls.subsumption_cycles(ruleset))
        return errors, warnings

    @classmethod
    def subsumption_cycles(cls, ruleset: FactRuleSet) -> list[str]:
        """Detect every cycle in the ``subsumes`` graph via DFS.

        Stronger than a pairwise mutual-suppression check (which only catches 2-cycles: A subsumes B
        and B subsumes A): a longer cycle (A subsumes B subsumes C subsumes A) would let ``prune``
        drop every term in the cycle whenever any one of them fires, which is never the intended
        suppression shape and is worth surfacing on its own even though the pairwise check in
        ruleset_problems also happens to catch the 2-cycle case.

        Args:
            ruleset (FactRuleSet): The rule set to inspect.

        Returns:
            list[str]: One description per distinct cycle found (deduplicated by member-set, so a
                cycle reached via two different DFS roots is reported once), or an empty list if the
                ``subsumes`` graph is acyclic.
        """
        reported_cycles: set[frozenset[str]] = set()
        cycle_descriptions: list[str] = []
        visit_state: dict[str, int] = {}  # 0 = unvisited (default), 1 = in progress, 2 = done
        path: list[str] = []

        def visit(term: str) -> None:
            visit_state[term] = 1
            path.append(term)
            rule = ruleset.by_term.get(term)
            if rule is not None:
                for subsumed_term in rule.subsumes:
                    if subsumed_term not in ruleset.by_term:
                        continue  # unknown target: already reported by validate_ruleset
                    subsumed_state = visit_state.get(subsumed_term, 0)
                    if subsumed_state == 1:
                        cycle_start_index = path.index(subsumed_term)
                        cycle_terms = tuple(path[cycle_start_index:])
                        dedup_key = frozenset(cycle_terms)
                        if dedup_key not in reported_cycles:
                            reported_cycles.add(dedup_key)
                            cycle_descriptions.append(
                                "subsumption cycle: " + " -> ".join((*cycle_terms, subsumed_term))
                            )
                    elif subsumed_state == 0:
                        visit(subsumed_term)
            path.pop()
            visit_state[term] = 2

        for term in ruleset.by_term:
            if visit_state.get(term, 0) == 0:
                visit(term)
        return cycle_descriptions

    @classmethod
    def missing_when_terms(cls) -> list[str]:
        """Cross-check EXCLUDED_TERMS against the catalogue.

        ``build_fact_ruleset`` silently skips any catalogue entry without a ``when``, on the
        documented assumption that every such entry is listed in ``EXCLUDED_TERMS``. This check makes
        that assumption checked rather than trusted: a catalogue edit that adds a ``when``-less term
        without updating ``EXCLUDED_TERMS`` (or vice versa) would otherwise silently under- or
        over-cover the ruleset.

        Returns:
            list[str]: One description per catalogue entry whose ``when``-presence disagrees with its
                ``EXCLUDED_TERMS`` membership, or an empty list if they agree everywhere.
        """
        problems: list[str] = []
        for definition in load_consequence_definitions():
            has_when = definition.when is not None
            is_documented_exclusion = definition.term in EXCLUDED_TERMS
            if has_when and is_documented_exclusion:
                problems.append(f"{definition.term} has a 'when' expression but is listed in EXCLUDED_TERMS")
            if not has_when and not is_documented_exclusion:
                problems.append(
                    f"{definition.term} has no 'when' expression but is not listed in EXCLUDED_TERMS"
                )
        return problems

    @classmethod
    def subsumes_so_parent_drift(cls, ruleset: FactRuleSet) -> list[str]:
        """Report where the declared ``subsumes`` graph disagrees with the SO ``so_supported_parents``.

        These are two genuinely different relations that happen to overlap a lot: ``subsumes`` is
        VEP's empirical suppression behaviour (term A, if it fires, hides term B), while
        ``so_supported_parents`` is the Sequence Ontology's ``is_a`` parent graph (term A *is a kind
        of* term B), restricted to parents also represented in this catalogue. A rule can subsume a
        term it is not an SO child of (e.g. suppressing a sibling), and a rule can have an SO parent it
        does not need to suppress (the parent's own ``when`` may already be false whenever the child's
        is true, with no pruning required). So this report is informational context for a catalogue
        reviewer, not a pass/fail check -- drift here is common and often fine; it is a prompt to look,
        not a bug report.

        Args:
            ruleset (FactRuleSet): The rule set to report drift for.

        Returns:
            list[str]: One description per rule with at least one asymmetry between its ``subsumes``
                set and its SO-parent-derived term set, or an empty list if every rule's two sets
                match exactly.
        """
        all_definitions = load_consequence_definitions()
        definitions_by_term = {definition.term: definition for definition in all_definitions}
        term_by_accession = {definition.accession: definition.term for definition in all_definitions}

        drift: list[str] = []
        for rule in ruleset.rules:
            definition = definitions_by_term.get(rule.term)
            if definition is None:
                continue
            so_parent_terms = {
                term_by_accession[accession]
                for accession in definition.so_supported_parents
                if accession in term_by_accession
            }
            subsumed_terms = set(rule.subsumes)
            subsumes_without_so_parent = sorted(subsumed_terms - so_parent_terms)
            so_parent_without_subsumes = sorted(so_parent_terms - subsumed_terms)
            if subsumes_without_so_parent or so_parent_without_subsumes:
                drift.append(
                    f"{rule.term}: subsumes-only={subsumes_without_so_parent}, "
                    f"so-parent-only={so_parent_without_subsumes}"
                )
        return drift

    @classmethod
    def term_and_rank_uniqueness(cls, entries: list[dict[str, Any]]) -> list[str]:
        """Check term/rank uniqueness across *every* entry, including ``when``-less ones.

        ruleset_problems only sees entries that made it into a FactRuleSet (i.e. have a ``when``), so
        a duplicate rank between two ``when``-less entries -- or between a ``when``-less entry and a
        compiled rule -- would otherwise go unnoticed.

        Args:
            entries (list[dict[str, Any]]): The candidate catalogue's ``consequences`` list.

        Returns:
            list[str]: One description per duplicate term or rank found, or an empty list if every
                term and rank is unique.
        """
        errors: list[str] = []
        seen_terms: set[str] = set()
        for entry in entries:
            term = entry["term"]
            if term in seen_terms:
                errors.append(f"duplicate catalogue term: {term}")
            seen_terms.add(term)

        seen_ranks: dict[int, str] = {}
        for entry in entries:
            rank = int(entry.get("rank", 0))
            if rank in seen_ranks:
                errors.append(f"duplicate catalogue rank {rank}: {entry['term']} and {seen_ranks[rank]}")
            seen_ranks[rank] = entry["term"]

        return errors


# --------------------------------------------------------------------------
# Tier 2: per-rule enumeration and per-pair propositional search.
#
# A rule's atom *support* is its formula's own atoms, extended by the requires-closure (every
# atom that some support atom's declared `requires` would force, transitively) and by a small
# hand-declared exclusivity/equivalence axiom set covering relationships the fact layer's
# `requires` mechanism does not itself express (e.g. "is_insertion and is_deletion cannot both
# hold" is true by construction of the two fact functions, but neither fact declares the other in
# its `requires`, since `requires` only expresses "if X is False, so am I", not general mutual
# exclusion). Enumerating every assignment over that support and discarding the ones the axioms
# rule out gives every *possible* combination of fact values a real Context could produce, which
# is exactly what "is this rule ever true" needs to be sound. See the module docstring for why
# this -- not equality saturation -- is the right tool for that question.
# --------------------------------------------------------------------------

_EXCLUSIVITY_AXIOM_PAIRS: tuple[tuple[str, str], ...] = (
    # alternate_allele_length > reference_allele_length and reference_allele_length >
    # alternate_allele_length cannot both hold.
    ("is_insertion", "is_deletion"),
    # is_upstream/is_downstream both explicitly return False whenever the variant overlaps the
    # transcript at all (facts/position.py).
    ("overlap_feature", "is_upstream"),
    ("overlap_feature", "is_downstream"),
    # overlap_exon requires overlap_feature, and is_upstream/is_downstream require the opposite,
    # so this pair is implied by the two axioms above.
    ("overlap_exon", "is_upstream"),
    ("overlap_exon", "is_downstream"),
    # A fixed transcript has one 5' side and one 3' side; a variant that misses the transcript
    # cannot be simultaneously "before the 5' end" and "after the 3' end" of the same feature.
    ("is_upstream", "is_downstream"),
)
"""Atom pairs that can never both be True, true by construction of the fact functions involved
and not an approximation"""

_EQUIVALENCE_AXIOM_PAIRS: tuple[tuple[str, str], ...] = ()
"""Atom pairs that always carry the same value (duplicate fact formulas under different names)."""


@frozen(kw_only=True, auto_attribs=True)
class FactAssignmentSpace:
    """Define consistent Boolean fact assignments over an atom support.

    Centralises ``requires`` closure, declared axioms, enumeration, and satisfiability searches
    used by single-rule and pairwise validation.

    Attributes:
        support (tuple[str, ...]): The atom names this space is defined over, sorted.
    """

    support: tuple[str, ...]

    @classmethod
    def for_rule(cls, rule: FactRule) -> FactAssignmentSpace:
        """Return the space a rule's formula must be evaluated over to be judged soundly.

        Args:
            rule (FactRule): The rule to compute a support for.

        Returns:
            FactAssignmentSpace: The space over ``rule.atoms()``, extended by the ``requires`` closure
                and by any declared axiom atoms that touch it, sorted for a deterministic
                enumeration order.
        """
        support = cls._extend_with_relevant_axiom_atoms(cls._requires_closure(rule.atoms()))
        return cls(support=tuple(sorted(support)))

    def union(self, other: FactAssignmentSpace) -> FactAssignmentSpace:
        """Return the space covering both supports, for a question about a pair of rules."""
        return type(self)(support=tuple(sorted(set(self.support) | set(other.support))))

    def allows(self, assignment: dict[str, bool]) -> bool:
        """Return whether an assignment over this support violates no ``requires``/axiom.

        Args:
            assignment (dict[str, bool]): A full assignment over some atom support.

        Returns:
            bool: False if the assignment sets some atom True while a declared ``requires`` of it is
                False in the same assignment, or violates an exclusivity/equivalence axiom; True
                otherwise (i.e. this is an assignment a real Context could actually produce, as far as
                the declared constraints can tell).
        """
        for atom_name, value in assignment.items():
            if not value:
                continue
            registered_fact = REGISTRY.get(atom_name)
            if registered_fact is None:
                continue
            for required_atom in registered_fact.requires:
                if assignment.get(required_atom) is False:
                    return False
        for atom_a, atom_b in _EXCLUSIVITY_AXIOM_PAIRS:
            if assignment.get(atom_a) and assignment.get(atom_b):
                return False
        for atom_a, atom_b in _EQUIVALENCE_AXIOM_PAIRS:
            value_a, value_b = assignment.get(atom_a), assignment.get(atom_b)
            if value_a is not None and value_b is not None and value_a != value_b:
                return False
        return True

    def consistent_assignments(self) -> Iterator[dict[str, bool]]:
        """Yield every assignment over this support consistent with ``requires``/the declared axioms.

        Yields:
            dict[str, bool]: One consistent full assignment per iteration.
        """
        for bitmask in range(1 << len(self.support)):
            assignment = {
                atom_name: bool(bitmask & (1 << index)) for index, atom_name in enumerate(self.support)
            }
            if self.allows(assignment):
                yield assignment

    def satisfies(self, expression: Expr) -> bool:
        """Return whether an expression has a consistent assignment over ``support``.

        This searches different combinations of values of atoms making up an expression to check
        whether given the logical constraints of the expression, it is possible for the expression
        to be satisfied. It is just checking if it is possible for the expression to be satisfied,
        so it stops as soon as it finds a valid set of assignments.

        Args:
            expression (Expr): The expression whose satisfiability is being checked.

        Returns:
            bool: True if at least one consistent assignment satisfies ``expression``.
        """
        support_set = frozenset(self.support)
        forced_values: dict[tuple[str, bool], list[tuple[str, bool]]] = {}

        def add_forced_value(source: str, source_value: bool, target: str, target_value: bool) -> None:
            if source in support_set and target in support_set:
                forced_values.setdefault((source, source_value), []).append((target, target_value))

        for atom_name in self.support:
            registered_fact = REGISTRY.get(atom_name)
            if registered_fact is None:
                continue
            for required_atom in registered_fact.requires:
                add_forced_value(atom_name, True, required_atom, True)
                add_forced_value(required_atom, False, atom_name, False)

        for atom_a, atom_b in _EXCLUSIVITY_AXIOM_PAIRS:
            add_forced_value(atom_a, True, atom_b, False)
            add_forced_value(atom_b, True, atom_a, False)

        for atom_a, atom_b in _EQUIVALENCE_AXIOM_PAIRS:
            for value in (False, True):
                add_forced_value(atom_a, value, atom_b, value)
                add_forced_value(atom_b, value, atom_a, value)

        def extend(assignment: dict[str, bool], atom_name: str, value: bool) -> dict[str, bool] | None:
            extended = assignment.copy()
            pending = [(atom_name, value)]

            while pending:
                pending_atom, pending_value = pending.pop()
                if pending_atom in extended:
                    if extended[pending_atom] != pending_value:
                        return None
                    continue

                extended[pending_atom] = pending_value
                pending.extend(forced_values.get((pending_atom, pending_value), ()))

            return extended

        expression_atoms = sorted(
            expression.atoms(),
            key=lambda atom_name: len(forced_values.get((atom_name, True), ())),
            reverse=True,
        )

        def search(assignment: dict[str, bool]) -> bool:
            value = expression.evaluate_partial(assignment)
            if value is not None:
                return value

            for atom_name in expression_atoms:
                if atom_name in assignment:
                    continue
                for atom_value in (True, False):
                    extended = extend(assignment, atom_name, atom_value)
                    if extended is not None and search(extended):
                        return True
                return False

            # A custom Expr may remain unknown after all of its declared atoms are assigned. Preserve
            # the exhaustive enumerator's behaviour for that uncommon extension point.
            return any(expression.evaluate_partial(candidate) for candidate in self.consistent_assignments())

        return search({})

    @staticmethod
    def _requires_closure(atom_names: Iterable[str]) -> frozenset[str]:
        """Extend a set of atom names with everything their declared ``requires`` would pull in.

        Args:
            atom_names (Iterable[str]): The starting atom names.

        Returns:
            frozenset[str]: ``atom_names`` plus every atom reachable by following ``requires`` edges
                to a fixed point.
        """
        support = set(atom_names)
        changed = True
        while changed:
            changed = False
            for atom_name in tuple(support):
                registered_fact = REGISTRY.get(atom_name)
                if registered_fact is None:
                    continue
                for required_atom in registered_fact.requires:
                    if required_atom not in support:
                        support.add(required_atom)
                        changed = True
        return frozenset(support)

    @staticmethod
    def _extend_with_relevant_axiom_atoms(support: frozenset[str]) -> frozenset[str]:
        """Pull in the other half of any declared axiom pair that already touches ``support``.

        Args:
            support (frozenset[str]): The atom support so far.

        Returns:
            frozenset[str]: ``support`` extended with the partner atom of every
                _EXCLUSIVITY_AXIOM_PAIRS/_EQUIVALENCE_AXIOM_PAIRS entry that already
                has at least one member in ``support``.
        """
        extended = set(support)
        for atom_a, atom_b in (*_EXCLUSIVITY_AXIOM_PAIRS, *_EQUIVALENCE_AXIOM_PAIRS):
            if atom_a in extended or atom_b in extended:
                extended.add(atom_a)
                extended.add(atom_b)
        return frozenset(extended)


_MAX_UNION_SUPPORT_FOR_PAIR_ENUMERATION = 18
"""Maximum combined atom support for a CONTRADICTORY or REDUNDANT pair search.

The search propagates constraints and stops as soon as the expression is settled, but remains
exponential in the worst case. The cap lets a future, larger catalogue skip oversized pairs
instead of allowing an unbounded search.
"""


@frozen(kw_only=True, auto_attribs=True)
class RuleReachability:
    """Represent the Tier 2 reachability result for one rule.

    Attributes:
        term (str): The rule's SO term.
        support (tuple[str, ...]): The atom support (see atom_support_for_rule) the rule
            was enumerated over.
        valid_assignment_count (int): Number of assignments over ``support`` consistent with the
            declared ``requires``/axioms.
        satisfying_assignment_count (int): Of those, the number for which the rule's ``when``
            evaluates True.
    """

    term: str
    support: tuple[str, ...]
    valid_assignment_count: int
    satisfying_assignment_count: int

    @property
    def is_reachable(self) -> bool:
        """Return whether a valid assignment satisfies the rule."""
        return self.satisfying_assignment_count > 0

    @classmethod
    def for_rule(cls, rule: FactRule) -> RuleReachability:
        """Determine whether a rule is reachable within its fact-assignment space.

        Args:
            rule (FactRule): The rule to check.

        Returns:
            RuleReachability: The enumeration result. ``evaluate_partial`` never returns None here,
                because ``support`` is constructed to be a superset of ``rule.atoms()``.
        """
        space = FactAssignmentSpace.for_rule(rule)
        valid_assignment_count = 0
        satisfying_assignment_count = 0
        for assignment in space.consistent_assignments():
            valid_assignment_count += 1
            if rule.when.evaluate_partial(assignment):
                satisfying_assignment_count += 1
        return cls(
            term=rule.term,
            support=space.support,
            valid_assignment_count=valid_assignment_count,
            satisfying_assignment_count=satisfying_assignment_count,
        )

    @classmethod
    def unreachable_in(cls, ruleset: FactRuleSet) -> list[RuleReachability]:
        """Return every rule whose ``when`` formula is unsatisfiable under the Tier 2 axioms.

        An unreachable rule is dead code: it can never contribute a term to any output, which -- unlike
        the CONTRADICTORY/REDUNDANT findings below -- is a genuine defect rather than informational
        context, so callers (the validation script, the test suite) should treat a non-empty result as
        a build-breaking error.

        Args:
            ruleset (FactRuleSet): The rule set to check.

        Returns:
            list[RuleReachability]: The unreachable rules' results (``is_reachable`` False for all of
                them), or an empty list if every rule is reachable.
        """
        return [result for rule in ruleset.rules if not (result := cls.for_rule(rule)).is_reachable]


@frozen(kw_only=True, auto_attribs=True)
class RulePairFinding:
    """Record one pairwise rule-validation finding.

    Attributes:
        kind (str): Either ``"CONTRADICTORY"`` (the two rules' formulas can never both hold) or
            ``"REDUNDANT"`` (rule ``term_a`` subsumes ``term_b`` *and* ``term_a``'s formula
            logically entails ``term_b``'s, confirming the subsumption is logically grounded).
            Both kinds are informational governance signal, not build-breaking errors -- see
            RulePairFinding.contradictory_in/RulePairFinding.redundant_in.
        term_a (str): The first rule's term.
        term_b (str): The second rule's term.
        message (str): Human-readable description.
    """

    kind: str
    term_a: str
    term_b: str
    message: str

    @classmethod
    def contradictory_in(
        cls, ruleset: FactRuleSet, *, max_union_support: int = _MAX_UNION_SUPPORT_FOR_PAIR_ENUMERATION
    ) -> list[RulePairFinding]:
        """Find rule pairs whose formulas can never both hold, and that do not subsume each other.

        This is informational rather than an error: two rules never co-firing is often simply how the
        ruleset is meant to work (e.g. ``stop_gained`` and ``synonymous_variant`` are naturally
        mutually exclusive), not a defect. It is useful context for a catalogue reviewer -- a pair
        that *should* be able to co-fire but this check finds contradictory would be worth a second
        look -- rather than something a script should fail the build over.

        Args:
            ruleset (FactRuleSet): The rule set to check every pair of.
            max_union_support (int): Skip a pair whose combined atom support exceeds this size, so a
                future larger catalogue degrades to "fewer pairs checked" rather than a slow/unbounded
                enumeration.

        Returns:
            list[RulePairFinding]: One ``"CONTRADICTORY"`` finding per such pair.
        """
        findings: list[RulePairFinding] = []
        rules = ruleset.rules
        spaces = {rule.term: FactAssignmentSpace.for_rule(rule) for rule in rules}

        for index_a, rule_a in enumerate(rules):
            for index_b in range(index_a + 1, len(rules)):
                rule_b = rules[index_b]
                if rule_b.term in rule_a.subsumes or rule_a.term in rule_b.subsumes:
                    continue
                union = spaces[rule_a.term].union(spaces[rule_b.term])
                if len(union.support) > max_union_support:
                    continue
                can_both_hold = union.satisfies(And(rule_a.when, rule_b.when))
                if not can_both_hold:
                    findings.append(
                        RulePairFinding(
                            kind="CONTRADICTORY",
                            term_a=rule_a.term,
                            term_b=rule_b.term,
                            message=f"{rule_a.term} and {rule_b.term} never co-fire under the Tier 2 axioms",
                        )
                    )
        return findings

    @classmethod
    def redundant_in(
        cls, ruleset: FactRuleSet, *, max_union_support: int = _MAX_UNION_SUPPORT_FOR_PAIR_ENUMERATION
    ) -> list[RulePairFinding]:
        """Find declared ``subsumes`` edges whose entailment is confirmed by the formulas themselves.

        For each declared ``A subsumes B``, check whether ``A.when`` logically entails ``B.when`` over
        every valid assignment of their combined support (i.e. whenever A's formula could possibly be
        true, B's formula is necessarily true too). Where it holds, the subsumption is *logically
        redundant* in the sense that A firing already guarantees B's condition holds -- the ``subsumes``
        declaration is provably sound, not merely asserted. Where it does not hold under this bounded
        axiom set, that is *not* reported as an error: most of this catalogue's subsumption edges rely
        on domain knowledge broader than the small declared axiom set here (e.g. "a deletion spanning
        the whole transcript necessarily also overlaps every intron and UTR" is true, but nothing in
        _EXCLUSIVITY_AXIOM_PAIRS encodes it), so a pair not appearing in this report is
        inconclusive, not suspect.

        Args:
            ruleset (FactRuleSet): The rule set to check every declared ``subsumes`` edge of.
            max_union_support (int): Skip an edge whose combined atom support exceeds this size.

        Returns:
            list[RulePairFinding]: One ``"REDUNDANT"`` finding per confirmed-entailed edge.
        """
        findings: list[RulePairFinding] = []
        spaces = {rule.term: FactAssignmentSpace.for_rule(rule) for rule in ruleset.rules}

        for rule_a in ruleset.rules:
            for subsumed_term in rule_a.subsumes:
                rule_b = ruleset.by_term.get(subsumed_term)
                if rule_b is None:
                    continue
                union = spaces[rule_a.term].union(spaces[rule_b.term])
                if len(union.support) > max_union_support:
                    continue
                counterexample = And(rule_a.when, Not(rule_b.when))
                entailed = not union.satisfies(counterexample)
                if entailed:
                    findings.append(
                        RulePairFinding(
                            kind="REDUNDANT",
                            term_a=rule_a.term,
                            term_b=subsumed_term,
                            message=(
                                f"{rule_a.term}'s when provably entails {subsumed_term}'s when -- the "
                                "declared subsumption is logically sound, confirmed under the Tier 2 axioms"
                            ),
                        )
                    )
        return findings


@frozen(kw_only=True, auto_attribs=True)
class ValidationReport:
    """Combine structural and bounded-search validation results.

    Attributes:
        errors (list[str]): From validate_ruleset; non-empty means the ruleset is broken.
        warnings (list[str]): From validate_ruleset; expected rank inversions.
        missing_when_problems (list[str]): From report_missing_when_terms; non-empty means
            the catalogue and ``EXCLUDED_TERMS`` have drifted apart -- treat as an error.
        so_parent_drift (list[str]): From report_subsumes_so_parent_drift; informational.
        unreachable_rules (list[RuleReachability]): From RuleReachability.unreachable_in; non-empty
            means a rule can never fire -- treat as an error.
        contradictory_pairs (list[RulePairFinding]): From RulePairFinding.contradictory_in;
            informational.
        redundant_pairs (list[RulePairFinding]): From RulePairFinding.redundant_in;
            informational.
    """

    errors: list[str]
    warnings: list[str]
    missing_when_problems: list[str]
    so_parent_drift: list[str]
    unreachable_rules: list[RuleReachability]
    contradictory_pairs: list[RulePairFinding]
    redundant_pairs: list[RulePairFinding] = field(factory=list)

    @property
    def has_errors(self) -> bool:
        """Return whether the report contains a build-breaking finding."""
        return bool(self.errors or self.missing_when_problems or self.unreachable_rules)

    @classmethod
    def for_ruleset(cls, ruleset: FactRuleSet) -> ValidationReport:
        """Run every Tier 1 + Tier 2 check over a rule set. Pure Python, always on; no z3 required.

        Args:
            ruleset (FactRuleSet): The rule set to validate.

        Returns:
            ValidationReport: The combined result. See ValidationReport.has_errors for the
                build-breaking subset.
        """
        errors, warnings = StructuralChecks.ruleset_problems(ruleset)
        return cls(
            errors=errors,
            warnings=warnings,
            missing_when_problems=StructuralChecks.missing_when_terms(),
            so_parent_drift=StructuralChecks.subsumes_so_parent_drift(ruleset),
            unreachable_rules=RuleReachability.unreachable_in(ruleset),
            contradictory_pairs=RulePairFinding.contradictory_in(ruleset),
            redundant_pairs=RulePairFinding.redundant_in(ruleset),
        )

    @classmethod
    def for_catalogue_entries(
        cls, entries: list[dict[str, Any]], baseline_entries: list[dict[str, Any]]
    ) -> ValidationReport:
        """Validate a candidate ``consequences`` list before it is ever written to disk.

        Reuses ensembl.vep.annotators.consequence.rules.validation's structural (Tier 1) and enumeration
        (Tier 2) checks over a FactRuleSet built from the *candidate* entries (see
        _ruleset_from_entries) -- not the packaged catalogue, which is what makes this safe to
        call before a real term has been added or removed on disk. The two checks that hardcode a call
        to the packaged load_consequence_definitions() (report_missing_when_terms,
        report_subsumes_so_parent_drift) are intentionally not reused here, since they are defined in
        terms of the checked-in EXCLUDED_TERMS/catalogue rather than a candidate in-memory one; a
        needs-review, ``when``-less addition is instead surfaced separately, at add-time, by
        resolve_when's NEEDS REVIEW banner.

        Unreachable-rule findings are diffed against ``baseline_entries`` (the catalogue before this
        add/remove was applied) rather than treated as build-breaking outright: the checked-in
        catalogue already carries one known, pre-existing dead rule (``feature_truncation``). Treating
        every unreachable rule as build-breaking here would permanently refuse every write against the
        real catalogue, for a reason unrelated to whatever term is being added or removed. A rule
        unreachable in *both* the baseline and the candidate is downgraded to an informational warning;
        one that only became unreachable in the candidate is a genuine regression introduced by this
        operation and stays build-breaking.

        Args:
            entries (list[dict[str, Any]]): The candidate catalogue's ``consequences`` list.
            baseline_entries (list[dict[str, Any]]): The catalogue's ``consequences`` list *before*
                this operation, used only to distinguish pre-existing unreachable rules from new ones.

        Returns:
            ValidationReport: errors/unreachable_rules non-empty means refuse to write; warnings are
                informational only (see ValidationReport.has_errors).
        """
        ruleset = cls._ruleset_from_entries(entries)
        structural_errors, warnings = StructuralChecks.ruleset_problems(ruleset)
        errors = StructuralChecks.term_and_rank_uniqueness(entries) + structural_errors

        baseline_ruleset = cls._ruleset_from_entries(baseline_entries)
        baseline_unreachable_terms = frozenset(
            result.term for result in RuleReachability.unreachable_in(baseline_ruleset)
        )
        all_unreachable: list[RuleReachability] = RuleReachability.unreachable_in(ruleset)
        new_unreachable = [
            result for result in all_unreachable if result.term not in baseline_unreachable_terms
        ]
        warnings = list(warnings) + [
            f"pre-existing unreachable rule (not introduced by this change, not build-breaking here): "
            f"{result.term}"
            for result in all_unreachable
            if result.term in baseline_unreachable_terms
        ]

        return cls(
            errors=errors,
            warnings=warnings,
            missing_when_problems=[],
            so_parent_drift=[],
            unreachable_rules=new_unreachable,
            contradictory_pairs=[],
            redundant_pairs=[],
        )

    @classmethod
    def _ruleset_from_entries(cls, entries: list[dict[str, Any]]) -> FactRuleSet:
        """Compile a rule set from raw catalogue entry dicts.

        Parses each entry into a ConsequenceDefinition and hands them to the one compilation path, so
        a candidate cannot be validated against a rule set built differently from the one the
        annotators will evaluate.
        """
        definitions = []
        for entry in entries:
            impact = entry.get("impact")
            definitions.append(
                ConsequenceDefinition(
                    term=entry["term"],
                    accession=entry["accession"],
                    rank=int(entry.get("rank", 0)),
                    tier=int(entry.get("tier", 3)),
                    impact=ConsequenceImpact(impact) if impact else None,
                    when=entry.get("when"),
                    subsumes=tuple(entry.get("subsumes", ())),
                    so_supported_parents=tuple(entry.get("so_supported_parents", ())),
                    notes=entry.get("notes"),
                )
            )
        return build_fact_ruleset_from_definitions(definitions)


# --------------------------------------------------------------------------
# Tier 3: whole-ruleset SAT queries via z3 (optional dev dependency).
#
# Same two query shapes as Tier 2 (reachability, subsumption-edge soundness), same axiom set, but
# posed to a solver instead of enumerated -- so a query is not limited to one rule/pair's bounded
# atom support the way Tier 2's brute force is. z3 is imported lazily, only when this function is
# actually called (--solver=z3), so it stays an optional dependency
# ([project.optional-dependencies] validation = ["z3-solver"] in pyproject.toml) rather than a
# hard one.
# --------------------------------------------------------------------------


class SolverChecks:
    """Run whole-rule-set validation through z3."""

    @classmethod
    def findings(cls, ruleset: FactRuleSet) -> list[str]:
        """Run the whole-ruleset reachability and subsumption-soundness queries via z3.

        Args:
            ruleset (FactRuleSet): The rule set to check.

        Returns:
            list[str]: Findings, prefixed ``"UNREACHABLE (z3): ..."`` (a rule whose formula is
                unsatisfiable under the axioms -- a real defect) or ``"NOT ENTAILED (z3): ..."`` (a
                declared ``subsumes`` edge not backed by formula entailment -- informational, see
                RulePairFinding.redundant_in's docstring for why that is not automatically a bug).

        Raises:
            RuntimeError: If ``z3`` is not installed. Install the optional ``validation`` extra
                (``pip install '.[validation]'`` / ``uv sync --extra validation``) to use this.
        """
        try:
            import z3  # pylint: disable=import-outside-toplevel
        except ImportError as import_error:
            raise RuntimeError(
                "z3-solver is not installed. Install the optional 'validation' extra "
                "(e.g. `uv sync --extra validation`) to use --solver=z3."
            ) from import_error

        atom_names_in_ruleset = sorted({atom_name for rule in ruleset.rules for atom_name in rule.atoms()})
        z3_atoms = {atom_name: z3.Bool(atom_name) for atom_name in atom_names_in_ruleset}

        axioms: list[Any] = []
        for atom_name in atom_names_in_ruleset:
            registered_fact = REGISTRY.get(atom_name)
            if registered_fact is None:
                continue
            for required_atom in registered_fact.requires:
                if required_atom in z3_atoms:
                    axioms.append(z3.Implies(z3_atoms[atom_name], z3_atoms[required_atom]))
        for atom_a, atom_b in _EXCLUSIVITY_AXIOM_PAIRS:
            if atom_a in z3_atoms and atom_b in z3_atoms:
                axioms.append(z3.Not(z3.And(z3_atoms[atom_a], z3_atoms[atom_b])))
        for atom_a, atom_b in _EQUIVALENCE_AXIOM_PAIRS:
            if atom_a in z3_atoms and atom_b in z3_atoms:
                axioms.append(z3_atoms[atom_a] == z3_atoms[atom_b])

        findings: list[str] = []

        for rule in ruleset.rules:
            solver = z3.Solver()
            solver.add(*axioms)
            solver.add(cls._expression_to_z3(rule.when, z3_atoms, z3))
            if solver.check() == z3.unsat:
                findings.append(f"UNREACHABLE (z3): {rule.term} is unsatisfiable under the declared axioms")

        for rule in ruleset.rules:
            for subsumed_term in rule.subsumes:
                subsumed_rule = ruleset.by_term.get(subsumed_term)
                if subsumed_rule is None:
                    continue
                solver = z3.Solver()
                solver.add(*axioms)
                solver.add(cls._expression_to_z3(rule.when, z3_atoms, z3))
                solver.add(z3.Not(cls._expression_to_z3(subsumed_rule.when, z3_atoms, z3)))
                if solver.check() == z3.sat:
                    findings.append(
                        f"NOT ENTAILED (z3): {rule.term} subsumes {subsumed_term} but {rule.term}'s when "
                        f"does not entail {subsumed_term}'s when (informational)"
                    )
        return findings

    @classmethod
    def _expression_to_z3(cls, expression: Expr, atoms: dict[str, Any], z3_module: Any) -> Any:
        """Translate an Expr tree into a z3 boolean expression.

        Typed with ``Any`` throughout rather than a z3 stub: z3 is an optional dependency (see the
        module docstring), so this module must import cleanly whether or not it -- and therefore any
        type stub for it -- is installed.

        Args:
            expression (Expr): The formula to translate.
            atoms (dict[str, Any]): Fact atom name -> ``z3.BoolRef`` lookup table.
            z3_module (Any): The imported ``z3`` module (passed through rather than imported here
                again, since the caller already did the lazy import).

        Returns:
            Any: The equivalent ``z3.BoolRef``.

        Raises:
            TypeError: If ``expression`` is a node type this translator does not support.
        """
        if isinstance(expression, Atom):
            return atoms[expression.name]
        if isinstance(expression, Not):
            return z3_module.Not(cls._expression_to_z3(expression.child, atoms, z3_module))
        if isinstance(expression, And):
            return z3_module.And(
                *(cls._expression_to_z3(child, atoms, z3_module) for child in expression.children)
            )
        if isinstance(expression, Or):
            return z3_module.Or(
                *(cls._expression_to_z3(child, atoms, z3_module) for child in expression.children)
            )
        raise TypeError(f"unsupported expression node type: {type(expression).__name__}")
