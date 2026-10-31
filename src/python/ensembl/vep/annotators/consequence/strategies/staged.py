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
"""Evaluate cheap facts before requesting expensive facts.

The strategy computes atoms below ``DEFAULT_CHEAP_ATOM_MAX_CLOSURE_COST``, extends the assignment
through ``implied_false``, and partially evaluates each rule using three-valued logic.

A variant in an intron has ``overlap_exon = False``, which makes ``overlap_cds`` false, which settles
``stop_gained``, ``frameshift_variant``, ``missense_variant`` and the rest of the coding tier as
false without translating sequence.
"""

from __future__ import annotations

from collections.abc import Sequence
import logging
from typing import TYPE_CHECKING

from attrs import field, frozen

from ensembl.vep.annotators.consequence.strategies.base import ConsequenceAnnotator
from ensembl.vep.annotators.consequence.facts.base import closure_cost, get_fact, implied_false
from ensembl.vep.annotators.consequence.context.annotator import AnnotatorContext
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.registry import register_annotator
from ensembl.vep.annotators.consequence.rules.rulesets import (
    FactRule,
    FactRuleSet,
    build_fact_ruleset,
    prune,
    subsumes_topological_order,
)
from ensembl.vep.trace import AnnotationTrace, TraceRecorder, describe_overlap
from ensembl.vep.trace import reconstruct
from ensembl.vep.models.overlap import FeatureOverlap, OverlapConsequence, TranscriptVariationAllele

if TYPE_CHECKING:
    from ensembl.vep.caches.codon_cache import CodonCache
    from ensembl.vep.caches.zone_map import ZoneMapCache

__all__ = ["CheapAtomAlphabet", "DEFAULT_CHEAP_ATOM_MAX_CLOSURE_COST", "StagedFactsAnnotator"]

logger = logging.getLogger(__name__)

DEFAULT_CHEAP_ATOM_MAX_CLOSURE_COST = 6


@frozen(kw_only=True, auto_attribs=True)
class CheapAtomAlphabet:
    """Store fact atoms that are inexpensive enough to evaluate eagerly.

    Attributes:
        atom_names (tuple[str, ...]): Every atom referenced by any rule in the rule set this was
            built from whose closure_cost is at most the threshold it was built with, sorted in
            ascending closure cost order so the cheapest facts are computed first.
    """

    atom_names: tuple[str, ...] = field(default=())

    @classmethod
    def build(
        cls,
        rules: tuple[FactRule, ...],
        max_closure_cost: int = DEFAULT_CHEAP_ATOM_MAX_CLOSURE_COST,
    ) -> CheapAtomAlphabet:
        """Collect the cheap atoms referenced anywhere in ``rules``, cheapest first.

        Args:
            rules (tuple[FactRule, ...]): The rules to collect atoms from.
            max_closure_cost (int): Only atoms whose closure_cost is at most this count as cheap.

        Returns:
            CheapAtomAlphabet: The sorted alphabet of cheap atoms.
        """
        referenced_atom_names = {atom_name for rule in rules for atom_name in rule.atoms()}
        cheap_atom_names = (
            atom_name for atom_name in referenced_atom_names if closure_cost(atom_name) <= max_closure_cost
        )
        return cls(atom_names=tuple(sorted(cheap_atom_names, key=closure_cost)))


@register_annotator("consequence.staged")
class StagedFactsAnnotator(ConsequenceAnnotator[FactRuleSet]):
    """Evaluate all applicable terms while deferring expensive facts.

    Attributes:
        name (str): Registry name, "consequence.staged".
        ruleset (FactRuleSet): The compiled rule set this annotator evaluates.
        settled_early (int): Number of (rule, overlap) decisions settled from cheap facts alone,
            across every call so far.
        needed_expensive (int): Number of (rule, overlap) decisions that stayed undecided after
            staging and had to fall through to evaluate, across every call so far.
        cheap_atom_max_closure_cost (int): Closure-cost ceiling for the atoms computed up front.
        early_subsumption_skip (bool): Whether the subsumes-topological early-skip scan is
            enabled. Off by default.
        skipped_by_subsumption (int): Number of (rule, overlap) evaluations elided by the
            early-skip scan. Stays 0 when ``early_subsumption_skip`` is False.
    """

    name = "consequence.staged"

    def __init__(
        self,
        ruleset: FactRuleSet | None = None,
        annotator_context: AnnotatorContext | None = None,
        early_subsumption_skip: bool = False,
        *,
        cheap_atom_max_closure_cost: int = DEFAULT_CHEAP_ATOM_MAX_CLOSURE_COST,
        zone_map: ZoneMapCache | None = None,
        codon_cache: CodonCache | None = None,
    ) -> None:
        """Initialize the annotator.

        Args:
            ruleset (FactRuleSet | None): The rule set to evaluate, or None to build the default
                rule set from the checked-in consequence catalogue.
            annotator_context (AnnotatorContext | None): Run-scoped caches to evaluate against,
                or None (the default) for a fresh, private one.
            early_subsumption_skip (bool): Evaluate rules in subsumes-topological order and skip any rule
                already suppressed by one that evaluated True.
            cheap_atom_max_closure_cost (int): Closure-cost ceiling for the atoms computed up
                front, defaulting to DEFAULT_CHEAP_ATOM_MAX_CLOSURE_COST. Raising it settles more rules
                from the cheap pass but computes facts some overlaps never need, lowering it does the
                reverse. Output is unaffected either way
            zone_map (ZoneMapCache | None): Optional zone-map cache.
            codon_cache (CodonCache | None): Optional shared-codon cache.
        """
        super().__init__(
            ruleset if ruleset is not None else build_fact_ruleset(),
            annotator_context,
            zone_map=zone_map,
            codon_cache=codon_cache,
        )
        self.cheap_atom_max_closure_cost = cheap_atom_max_closure_cost
        self._cheap_atom_alphabet = CheapAtomAlphabet.build(self.ruleset.rules, cheap_atom_max_closure_cost)
        self.early_subsumption_skip = early_subsumption_skip
        self._subsumes_order = subsumes_topological_order(self.ruleset)
        self.settled_early = 0
        self.needed_expensive = 0
        self.skipped_by_subsumption = 0

    def _known_from_cheap_facts(self, fact_context: FactContext) -> dict[str, bool]:
        """Compute every cheap atom, then propagate the declared ``requires`` implications.

        Args:
            fact_context (FactContext): The overlap being evaluated.

        Returns:
            dict[str, bool]: The cheap atoms' values, extended with everything
                implied_false proves False from them..
        """
        cheap_known = {
            atom_name: bool(get_fact(fact_context, atom_name))
            for atom_name in self._cheap_atom_alphabet.atom_names
        }
        return implied_false(cheap_known)

    def _evaluate_with_known(
        self, fact_context: FactContext, known: dict[str, bool]
    ) -> list[OverlapConsequence]:
        """Evaluate every rule against a partial assignment, falling through when undecided.

        Also where ``last_trace`` is set, because ``known`` is what the forced-False facts are
        recovered from and this is the only place that holds it.

        Args:
            fact_context (FactContext): The overlap being evaluated.
            known (dict[str, bool]): A partial assignment, typically from _known_from_cheap_facts
                and optionally extended by a subclass.

        Returns:
            list[OverlapConsequence]: The consequences that apply, in rank order.
        """
        applicable_rules = self.ruleset.rules
        hit_terms: list[str] = []
        if self.early_subsumption_skip:
            suppressed: set[str] = set()
            fired: set[str] = set()
            applicable_terms = {rule.term for rule in applicable_rules}
            for rule in self._subsumes_order:
                if rule.term not in applicable_terms:
                    continue
                # Being suppressed is not on its own a licence to skip. prune() drops a term
                # whenever any *hit* term subsumes it, and a hit term that would get dropped
                # here still should contribute its own subsumes to that set.
                if rule.term in suppressed and suppressed.issuperset(rule.subsumes):
                    self.skipped_by_subsumption += 1
                    continue
                if self._rule_verdict(rule, fact_context, known):
                    fired.add(rule.term)
                    suppressed.update(rule.subsumes)
            # Restore rank order: the scan above ran in topological order, but the output
            # contract is rank-ordered terms.
            hit_terms = [rule.term for rule in applicable_rules if rule.term in fired]
        else:
            for rule in applicable_rules:
                if self._rule_verdict(rule, fact_context, known):
                    hit_terms.append(rule.term)
        kept_terms = prune(self.ruleset, hit_terms, apply_default=False)
        consequences = [self.ruleset.by_term[term].to_overlap_consequence() for term in kept_terms]

        self.last_trace = None
        if fact_context.recorder is not None:
            self.last_trace = AnnotationTrace(
                fact_events=reconstruct.memo_fact_events(fact_context)
                + reconstruct.forced_false_events(known, fact_context.memo),
                rule_events=reconstruct.partial_rule_events(applicable_rules, known, hit_terms),
                prune_events=reconstruct.subsumption_events(self.ruleset.by_term, hit_terms),
                emitted_terms=tuple(kept_terms),
                annotator_name=self.name,
                overlap_label=describe_overlap(fact_context.overlap),
            )
            logger.debug(self.last_trace.format_explanation())

        return consequences

    def _rule_verdict(self, rule: FactRule, fact_context: FactContext, known: dict[str, bool]) -> bool:
        """Resolve one rule's verdict, reaching the expensive fact tier only when it has to.

        Args:
            rule (FactRule): The rule to decide.
            fact_context (FactContext): The overlap being evaluated.
            known (dict[str, bool]): The partial assignment to try first.

        Returns:
            bool: Whether the rule's formula holds for this overlap.
        """
        verdict = rule.when.evaluate_partial(known)
        if verdict is None:
            self.needed_expensive += 1
            return bool(rule.when.evaluate(fact_context))
        self.settled_early += 1
        return bool(verdict)

    def annotate_feature_overlap(self, overlap: FeatureOverlap) -> list[OverlapConsequence]:
        """Calculate consequences for an overlap.

        Args:
            overlap (FeatureOverlap): The overlap to annotate. Must be a
                TranscriptVariationAllele.

        Returns:
            list[OverlapConsequence]: The consequences that apply, in rank order. Empty if
                nothing matched -- the caller applies a default/fallback term.

        Raises:
            TypeError: If ``overlap`` is not a TranscriptVariationAllele.
        """
        if not isinstance(overlap, TranscriptVariationAllele):
            raise TypeError(
                f"{type(self).__name__} only annotates TranscriptVariationAllele overlaps, got "
                f"{type(overlap).__name__}"
            )

        recorder = TraceRecorder() if logger.isEnabledFor(logging.DEBUG) else None
        fact_context = FactContext(
            overlap=overlap, recorder=recorder, annotator_context=self.annotator_context
        )
        self.annotator_context.prepare(fact_context)
        known = self._known_from_cheap_facts(fact_context)
        # _evaluate_with_known owns last_trace: it is the one that holds the partial assignment the
        # forced-False facts are recovered from.
        return self._evaluate_with_known(fact_context, known)

    def annotate_feature_overlaps(self, overlaps: Sequence[FeatureOverlap]) -> list[list[OverlapConsequence]]:
        """Calculate consequences for a batch of overlaps, one at a time.

        Args:
            overlaps (Sequence[FeatureOverlap]): The overlaps to annotate, in order.

        Returns:
            list[list[OverlapConsequence]]: One consequence list per input overlap, same order.
        """
        return [self.annotate_feature_overlap(overlap) for overlap in overlaps]
