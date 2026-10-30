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
"""Reconstruct annotation traces from completed fact evaluation.

The fact memo, rule set, and fired terms provide most of a derivation after annotation. Live
recording retains only read counts, cache origins, and the partial assignment used by staged
annotation because those values cannot be recovered by replay.

Fact origins are reconstructed as follows:

* **computed** -- in the memo, and no cache claimed it.
* **prefilled / cached** -- in the memo, and the recorder says which cache put it there.
* **forced False** -- False in the staged assignment and absent from the memo, so it was proved
  False by ``requires`` propagation without being computed.

Call ``memo_fact_events`` before reconstructing rules. Rule reconstruction reads facts again, so
capturing read counts afterwards would include replay activity.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.facts.base import REGISTRY
from ensembl.vep.annotators.consequence.rules.expressions import And, Atom, Expr, Not, Or
from ensembl.vep.annotators.consequence.rules.rulesets import FactRule
from ensembl.vep.trace.events import FactEvent, FactOrigin, PruneEvent, RuleEvent

__all__ = [
    "explain_partial_settlement",
    "forced_false_events",
    "guarded_rule_events",
    "memo_fact_events",
    "partial_rule_events",
    "plain_rule_events",
    "subsumption_events",
]


def memo_fact_events(fact_context: FactContext) -> tuple[FactEvent, ...]:
    """Return memoised facts in materialisation order.

    Dictionary insertion order records when each fact first entered the memo.

    Args:
        fact_context (FactContext): Evaluated context containing the memo and optional recorder.

    Returns:
        tuple[FactEvent, ...]: Fact events in materialisation order.
    """
    recorder = fact_context.recorder
    read_counts = recorder.read_counts if recorder is not None else {}
    served = recorder.served if recorder is not None else {}

    events = []
    for name, value in fact_context.memo.items():
        origin = served.get(name, FactOrigin.COMPUTED)
        registered_fact = REGISTRY.get(name)
        cost = registered_fact.cost if registered_fact is not None else 0
        events.append(
            FactEvent(
                name=name,
                value=value,
                cost=cost if origin is FactOrigin.COMPUTED else 0,
                origin=origin,
                reads=read_counts.get(name, 0),
            )
        )
    return tuple(events)


def forced_false_events(known: Mapping[str, bool], memo: Mapping[str, object]) -> tuple[FactEvent, ...]:
    """Return facts proved false through ``requires`` propagation.

    A false value in the closed partial assignment that is absent from the memo was settled
    without computation.

    Args:
        known (Mapping[str, bool]): Partial assignment used for rule evaluation.
        memo (Mapping[str, object]): Completed fact memo.

    Returns:
        tuple[FactEvent, ...]: Forced-false events with avoided cost and forcing atoms.
    """
    events = []
    for name, value in known.items():
        if value is not False or name in memo:
            continue
        registered_fact = REGISTRY.get(name)
        if registered_fact is None:
            continue
        forced_by = tuple(required for required in registered_fact.requires if known.get(required) is False)
        events.append(
            FactEvent(
                name=name,
                value=False,
                cost=registered_fact.cost,
                origin=FactOrigin.FORCED_FALSE,
                forced_by=forced_by,
            )
        )
    return tuple(events)


# pylint: disable-next=too-many-return-statements
def explain_partial_settlement(expr: Expr, known: dict[str, bool]) -> list[tuple[str, bool]]:
    """Return atom-value pairs that settled a partial expression.

    This follows ``Expr.evaluate_partial`` short-circuiting until it reaches the false conjunct or
    true disjunct that determined the result.

    Args:
        expr (Expr): Expression already settled by partial evaluation.
        known (dict[str, bool]): Partial assignment that settled the expression.

    Returns:
        list[tuple[str, bool]]: Settling atoms in evaluation order, or an empty list.
    """
    if isinstance(expr, Atom):
        value = known.get(expr.name)
        return [(expr.name, value)] if value is not None else []
    if isinstance(expr, Not):
        return explain_partial_settlement(expr.child, known)
    if isinstance(expr, And):
        for child in expr.children:
            if child.evaluate_partial(known) is False:
                return explain_partial_settlement(child, known)
        return []
    if isinstance(expr, Or):
        for child in expr.children:
            if child.evaluate_partial(known) is True:
                return explain_partial_settlement(child, known)
        return []
    return []


def plain_rule_events(rules: Sequence[FactRule], hit_terms: Sequence[str]) -> tuple[RuleEvent, ...]:
    """Return rule events for fully evaluated formulas.

    Use this for strategies without staging or a cheap-guard prefilter.

    Args:
        rules (Sequence[FactRule]): Rules considered, in rank order.
        hit_terms (Sequence[str]): Terms whose rules evaluated true.

    Returns:
        tuple[RuleEvent, ...]: Rule events in rank order.
    """
    hit_terms_set = set(hit_terms)
    events = []
    for rule in rules:
        events.append(
            RuleEvent(
                term=rule.term,
                rank=rule.rank,
                fired=rule.term in hit_terms_set,
                settled_early=False,
                reason=repr(rule.when),
            )
        )
    return tuple(events)


def guarded_rule_events(
    rules: Sequence[FactRule],
    candidate_terms: set[str],
    hit_terms: Sequence[str],
    guard_flags: Mapping[str, bool],
) -> tuple[RuleEvent, ...]:
    """Return rule events that distinguish guard rejection from full evaluation.

    A rule outside ``candidate_terms`` was rejected by ``FactRule.cheap_guard`` and therefore
    reports failed guard atoms instead of its unevaluated formula.

    Args:
        rules (Sequence[FactRule]): Rules in rank order.
        candidate_terms (set[str]): Terms whose rules passed the guard.
        hit_terms (Sequence[str]): Candidate terms whose rules evaluated true.
        guard_flags (Mapping[str, bool]): Guard values for the overlap.

    Returns:
        tuple[RuleEvent, ...]: Rule events in rank order.
    """
    hit_terms_set = set(hit_terms)
    events = []
    for rule in rules:
        if rule.term in candidate_terms:
            events.append(
                RuleEvent(
                    term=rule.term,
                    rank=rule.rank,
                    fired=rule.term in hit_terms_set,
                    settled_early=False,
                    reason=repr(rule.when),
                )
            )
            continue
        failed_atoms = sorted(atom for atom in rule.cheap_guard() if not guard_flags.get(atom, False))
        reason = (
            "guard failed: " + ", ".join(f"{atom}=False" for atom in failed_atoms)
            if failed_atoms
            else "guard failed"
        )
        events.append(
            RuleEvent(term=rule.term, rank=rule.rank, fired=False, settled_early=True, reason=reason)
        )
    return tuple(events)


def partial_rule_events(
    rules: Sequence[FactRule], known: dict[str, bool], hit_terms: Sequence[str]
) -> tuple[RuleEvent, ...]:
    """Return rule events from partial and complete verdicts.

    Rules settled by ``known`` report the determining atoms. Other rules report that they required
    full evaluation.

    Args:
        rules (Sequence[FactRule]): Rules in rank order.
        known (dict[str, bool]): Partial assignment used by staged evaluation.
        hit_terms (Sequence[str]): Terms whose complete verdict was true.

    Returns:
        tuple[RuleEvent, ...]: Rule events in rank order.
    """
    hit_terms_set = set(hit_terms)
    events = []
    for rule in rules:
        verdict = rule.when.evaluate_partial(known)
        fired = rule.term in hit_terms_set
        if verdict is None:
            events.append(
                RuleEvent(
                    term=rule.term,
                    rank=rule.rank,
                    fired=fired,
                    settled_early=False,
                    reason="fell through to full evaluation (needed expensive tier)",
                )
            )
            continue
        settling_atoms = explain_partial_settlement(rule.when, known)
        reason = (
            "settled early: " + ", ".join(f"{name}={value}" for name, value in settling_atoms)
            if settling_atoms
            else f"settled early: {rule.when!r}"
        )
        events.append(
            RuleEvent(term=rule.term, rank=rule.rank, fired=fired, settled_early=True, reason=reason)
        )
    return tuple(events)


def subsumption_events(by_term: Mapping[str, FactRule], hit_terms: Sequence[str]) -> tuple[PruneEvent, ...]:
    """Return the pruning decisions made by ``rules.rulesets.prune``.

    This read-only replay avoids adding trace state to the pruning path.

    Args:
        by_term (Mapping[str, FactRule]): Rules indexed by term.
        hit_terms (Sequence[str]): Terms whose rules fired before pruning.

    Returns:
        tuple[PruneEvent, ...]: Pruning events for dropped terms.
    """
    hit_terms_set = set(hit_terms)
    events = []
    for term in hit_terms:
        for subsumed_term in by_term[term].subsumes:
            if subsumed_term in hit_terms_set:
                events.append(PruneEvent(dropped_term=subsumed_term, subsumed_by=term))
    return tuple(events)
