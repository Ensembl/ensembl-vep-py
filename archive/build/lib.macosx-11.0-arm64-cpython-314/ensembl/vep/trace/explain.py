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
"""Comparing two derivations of the same overlap.

The question a developer has the day a conformance test fails is not *that* two annotators
disagree but *where*: which fact, which rule verdict, which subsumption edge they took
differently. Two traces of the same overlap answer it directly, because both record the whole
chain and not only its result.

``diff_terms`` is the same idea against a term list somebody else assigned. It is half of the
oracle diff the plan reserves for already-annotated input; the other half is ``output/parser/``
reading the CSQ terms out of the file, which does not exist yet. It works today against any term
list, which is what the ``--compare`` path uses when one side cannot trace.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from ensembl.vep.trace.annotation import AnnotationTrace
from ensembl.vep.trace.events import RuleEvent

__all__ = ["diff_terms", "diff_traces"]


def diff_traces(trace_a: AnnotationTrace, trace_b: AnnotationTrace) -> str:
    """Compare two derivations of one overlap and report where they first part company.

    Args:
        trace_a (AnnotationTrace): The first annotator's derivation.
        trace_b (AnnotationTrace): The second annotator's derivation.

    Returns:
        str: Whether the emitted terms agree, any fact both sides valued differently, the first
            rule (in rank order) whose verdict disagrees, and the rules only one side recorded.
    """
    name_a = trace_a.annotator_name or "A"
    name_b = trace_b.annotator_name or "B"
    lines = [f"comparing {name_a}  vs  {name_b}"]

    if trace_a.emitted_terms == trace_b.emitted_terms:
        lines.append(f"emitted terms AGREE: {list(trace_a.emitted_terms)}")
    else:
        lines.append(
            "emitted terms DIVERGE:\n"
            f"  {name_a}: {list(trace_a.emitted_terms)}\n"
            f"  {name_b}: {list(trace_b.emitted_terms)}"
        )

    lines.extend(_fact_value_divergences(trace_a, trace_b, name_a, name_b))
    lines.extend(_rule_divergences(trace_a, trace_b, name_a, name_b))
    return "\n".join(lines)


def _fact_value_divergences(
    trace_a: AnnotationTrace, trace_b: AnnotationTrace, name_a: str, name_b: str
) -> list[str]:
    """Report facts both derivations knew and disagreed on the value of.

    Only shared facts are compared. Two strategies reaching the same answer by different routes
    legitimately touch different *sets* of facts -- that is what staging is for -- so a fact only
    one side reached is not a divergence.
    """
    values_a = {event.name: event.value for event in trace_a.fact_events}
    values_b = {event.name: event.value for event in trace_b.fact_events}
    disagreements = []
    for name in sorted(set(values_a) & set(values_b)):
        if values_a[name] != values_b[name]:
            disagreements.append(f"  {name}: {name_a}={values_a[name]!r}  {name_b}={values_b[name]!r}")
    if not disagreements:
        return []
    return ["\nFACT VALUE DIVERGENCE:", *disagreements]


def _rule_divergences(
    trace_a: AnnotationTrace, trace_b: AnnotationTrace, name_a: str, name_b: str
) -> list[str]:
    """Report the first rule whose verdict differs, and any rule only one side recorded."""
    events_a = {event.term: event for event in trace_a.rule_events}
    events_b = {event.term: event for event in trace_b.rule_events}
    shared_terms = sorted(set(events_a) & set(events_b), key=lambda term: events_a[term].rank)

    lines = []
    first_divergent_term = next(
        (term for term in shared_terms if events_a[term].fired != events_b[term].fired), None
    )
    if first_divergent_term is None:
        lines.append("no rule-verdict divergence among the rules both traces recorded")
    else:
        event_a, event_b = events_a[first_divergent_term], events_b[first_divergent_term]
        lines.append(
            f"\nFIRST RULE-VERDICT DIVERGENCE: [{event_a.rank}] {first_divergent_term}\n"
            f"  {name_a}: {'FIRED' if event_a.fired else 'FALSE'} -- {event_a.reason}\n"
            f"  {name_b}: {'FIRED' if event_b.fired else 'FALSE'} -- {event_b.reason}"
        )

    only_in_a = sorted(set(events_a) - set(events_b))
    only_in_b = sorted(set(events_b) - set(events_a))
    if only_in_a:
        lines.append(f"rules traced only by {name_a}: {only_in_a}")
    if only_in_b:
        lines.append(f"rules traced only by {name_b}: {only_in_b}")
    return lines


def diff_terms(trace: AnnotationTrace, expected_terms: Sequence[str], source: str) -> str:
    """Explain how a derivation differs from the terms something else assigned to the same pair.

    The rule ledger is what makes this more than a set difference: for a term the other side
    emitted and we did not, this reports whether our rule was FALSE (and on what), or fired and
    was then subsumed. A term the catalogue does not model at all is reported as unexplainable
    rather than as a rule failure -- the catalogue is a documented subset of VEP, so part of any
    such diff is expected and correct, and this reports rather than judges.

    Args:
        trace (AnnotationTrace): Our derivation of the pair.
        expected_terms (Sequence[str]): The terms the other side assigned.
        source (str): What produced them, for the report, e.g. "consequence.linear" or "VEP CSQ".

    Returns:
        str: One line per disagreeing term, each carrying the reason from our own ledger.
    """
    ours = tuple(trace.emitted_terms)
    theirs = tuple(expected_terms)
    if set(ours) == set(theirs):
        return f"terms AGREE with {source}: {list(ours)}"

    events = {event.term: event for event in trace.rule_events}
    pruned_by = {event.dropped_term: event.subsumed_by for event in trace.prune_events}

    lines = [f"terms DIVERGE from {source}:"]
    for term in sorted(set(theirs) - set(ours)):
        lines.append(f"  only in {source}: {term} -- {_why_absent(term, events, pruned_by)}")
    for term in sorted(set(ours) - set(theirs)):
        event = events.get(term)
        reason = event.reason if event is not None else "no rule event recorded"
        lines.append(f"  only ours: {term} -- fired on {reason}")
    return "\n".join(lines)


def _why_absent(term: str, events: Mapping[str, RuleEvent], pruned_by: Mapping[str, str]) -> str:
    """Return why a term the other side emitted is missing from ours."""
    if term in pruned_by:
        return f"our rule fired, then was subsumed by {pruned_by[term]}"
    event = events.get(term)
    if event is None:
        return "not modelled by this catalogue (unexplainable, not a rule failure)"
    return f"our rule was FALSE -- {event.reason}"
