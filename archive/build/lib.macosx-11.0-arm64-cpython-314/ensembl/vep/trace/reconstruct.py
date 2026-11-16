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
"""Rebuilding a derivation after the fact, from what the run left behind.

The annotation path records almost nothing. What it leaves behind -- the ``FactContext`` memo
table, the rule set, the terms that fired -- is enough to say how every fact became known and why
every rule reached its verdict, so the explanation is derived here instead of being accumulated
event by event on the hottest function in the system.

Three things are *not* derivable, and only those are gathered live (``trace/events.py``'s
``TraceRecorder``): how many times each fact was read, which facts a cache served rather than the
fact layer computing, and -- for a staged annotator -- the partial assignment it settled rules
from. Everything else follows:

* **computed** -- in the memo, and no cache claimed it.
* **prefilled / cached** -- in the memo, and the recorder says which cache put it there.
* **forced False** -- False in the staged assignment and absent from the memo, so it was proved
  False by ``requires`` propagation without ever being computed. Recovering it this way is why
  nothing has to duplicate the propagation logic to explain it.

**Call ``memo_fact_events`` before any other reconstruction step.** The columnar annotators replay
their plan to recover per-rule verdicts, and the linear one recomputes its guard flags; both go
back through ``get_fact``, so read counts taken afterwards would double-count the real path's.
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
    """Return one event per memoised fact, in first-materialisation order.

    The memo is an ordinary dict, so its insertion order already *is* the order the derivation
    materialised its facts -- a prefilled fact appears where the cache put it, a computed one
    where the fact layer reached it.

    Args:
        fact_context (FactContext): The already-evaluated context. Its recorder supplies the read
            counts and the cache origins; with no recorder every fact reads as computed and unread.

    Returns:
        tuple[FactEvent, ...]: The memoised facts as events.
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
    """Return an event per fact ``requires`` propagation proved False without computing it.

    A staged annotator hands its rule loop a partial assignment already closed under the declared
    ``requires`` implications. Every name in it that is False and absent from the memo was settled
    by that closure alone, which is exactly the set to report -- and the reason no propagation
    logic has to be duplicated to report it.

    Args:
        known (Mapping[str, bool]): The partial assignment the annotator evaluated rules against.
        memo (Mapping[str, object]): The context's memo table.

    Returns:
        tuple[FactEvent, ...]: The forced-False facts, carrying the cost they avoided and the
            required atoms that forced them.
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
    """Return the (atom, value) pairs from ``known`` that settled ``expr``.

    Mirrors ``Expr.evaluate_partial``'s Kleene short-circuiting (``rules/expressions.py``) but
    returns *which* atoms caused the settlement -- the first False conjunct of an And, the first
    True disjunct of an Or -- walked through Not/And/Or until a settling Atom is reached.

    Takes a plain dict rather than a Mapping to match ``Expr.evaluate_partial``, which it
    mirrors and calls.

    Args:
        expr (Expr): A formula ``evaluate_partial(known)`` already settled.
        known (dict[str, bool]): The partial assignment that settled it.

    Returns:
        list[tuple[str, bool]]: The settling atoms and their values, in evaluation order. Empty if
            ``known`` did not settle ``expr``, or for a node type this does not recognise.
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
    """Return a verdict per rule for an annotator that evaluated every formula in full.

    For the strategies with no staging and no cheap-guard prefilter (``consequence.columnar``):
    every rule really was evaluated, so the reason is simply the formula.

    Args:
        rules (Sequence[FactRule]): Every rule considered, in rank order.
        hit_terms (Sequence[str]): The terms whose rule evaluated True.

    Returns:
        tuple[RuleEvent, ...]: One event per rule, in rank order.
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
    """Return a verdict per rule, telling a cheap-guard rejection apart from a real one.

    For the digest prefilter in ``consequence.linear`` and anything given a zone map: a rule
    outside ``candidate_terms`` never had its formula evaluated -- ``FactRule.cheap_guard()``
    rejected it from cheap facts alone -- so its reason names the guard atoms that failed rather
    than a formula evaluation that never happened.

    Args:
        rules (Sequence[FactRule]): Every rule in the rule set, in rank order.
        candidate_terms (set[str]): Terms whose rule passed the guard and was evaluated.
        hit_terms (Sequence[str]): The terms, among the candidates, whose rule evaluated True.
        guard_flags (Mapping[str, bool]): Every guard atom's value for this overlap.

    Returns:
        tuple[RuleEvent, ...]: One event per rule, in rank order.
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
    """Return a verdict per rule from a partial (Kleene) assignment plus the real verdicts.

    For the staged annotators: a rule settled by ``known`` alone never reached the expensive tier,
    so its reason names the settling atoms; one ``evaluate_partial`` could not settle fell through
    to full evaluation, and says so rather than repeating a formula it did not decide on.

    Args:
        rules (Sequence[FactRule]): Every rule in the rule set, in rank order.
        known (dict[str, bool]): The partial assignment staging settled from cheap and
            forced-False facts, exactly as passed to ``rule.when.evaluate_partial``.
        hit_terms (Sequence[str]): The terms whose real, fully-settled verdict was True.

    Returns:
        tuple[RuleEvent, ...]: One event per rule, in rank order.
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
    """Return the pruning decisions ``rules.rulesets.prune`` makes internally.

    A read-only replay of the same logic, so ``prune``'s signature and behaviour are untouched.

    Args:
        by_term (Mapping[str, FactRule]): ``FactRuleSet.by_term`` for the rule set the terms were
            matched against.
        hit_terms (Sequence[str]): Every term whose rule fired, pre-pruning.

    Returns:
        tuple[PruneEvent, ...]: One event per dropped term.
    """
    hit_terms_set = set(hit_terms)
    events = []
    for term in hit_terms:
        for subsumed_term in by_term[term].subsumes:
            if subsumed_term in hit_terms_set:
                events.append(PruneEvent(dropped_term=subsumed_term, subsumed_by=term))
    return tuple(events)
