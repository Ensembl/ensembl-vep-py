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
"""Evaluate consequence rules in rank order.

For each overlap, the strategy computes cheap guard facts, uses their digest to cache candidate
rules, evaluates candidates in rank order, and applies subsumption.

This strategy is the conformance reference. Every other consequence strategy must produce the same
ordered ``OverlapConsequence`` terms for the same input.
"""

from __future__ import annotations

from collections.abc import Sequence
import logging
from typing import TYPE_CHECKING

from attrs import field, frozen

from ensembl.vep.annotators.consequence.strategies.base import ConsequenceAnnotator
from ensembl.vep.annotators.registry import register_annotator
from ensembl.vep.annotators.consequence.facts.base import get_fact
from ensembl.vep.annotators.consequence.context.annotator import AnnotatorContext
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.rules.rulesets import FactRule, FactRuleSet, build_fact_ruleset, prune
from ensembl.vep.trace import AnnotationTrace, TraceRecorder, describe_overlap
from ensembl.vep.trace import reconstruct
from ensembl.vep.models.overlap import FeatureOverlap, OverlapConsequence, TranscriptVariationAllele

if TYPE_CHECKING:
    from ensembl.vep.caches.codon_cache import CodonCache
    from ensembl.vep.caches.zone_map import ZoneMapCache

__all__ = ["LinearFactsAnnotator"]

logger = logging.getLogger(__name__)


@frozen(kw_only=True, auto_attribs=True)
class _GuardDigestAlphabet:
    """The sorted set of fact atoms that can be used as cheap guards.

    Attributes:
        atom_names (tuple[str, ...]): Every cheap guard atom name referenced by any rule in the
            rule set this was built from, sorted for a stable digest string ordering.
    """

    atom_names: tuple[str, ...] = field(default=())

    @classmethod
    def build(cls, ruleset: FactRuleSet) -> _GuardDigestAlphabet:
        """Collect every cheap guard atom referenced by any rule in a rule set.

        Args:
            ruleset (FactRuleSet): The rule set to inspect.

        Returns:
            _GuardDigestAlphabet: The sorted alphabet of cheap guard atom names.
        """
        atom_names = sorted({atom_name for rule in ruleset.rules for atom_name in rule.cheap_guard()})
        return cls(atom_names=tuple(atom_names))


@register_annotator("consequence.linear")
class LinearFactsAnnotator(ConsequenceAnnotator[FactRuleSet]):
    """Evaluate a fact rule set with a rank-ordered linear scan.

    There is also a cheap guard to find shorter candidate rule lists, implemented by a digest-keyed
    candidate cache.

    Attributes:
        name (str): Registry name, "consequence.linear".
        ruleset (FactRuleSet): The compiled rule set this annotator evaluates.
    """

    name = "consequence.linear"

    def __init__(
        self,
        ruleset: FactRuleSet | None = None,
        annotator_context: AnnotatorContext | None = None,
        *,
        zone_map: ZoneMapCache | None = None,
        codon_cache: CodonCache | None = None,
    ) -> None:
        """Initialize the annotator.

        Args:
            ruleset (FactRuleSet | None): The rule set to evaluate, or None to build the default
                rule set from the checked-in consequence catalogue.
            annotator_context (AnnotatorContext | None): Run-scoped caches to evaluate against,
                or None (the default) for a fresh, private one.
            zone_map (ZoneMapCache | None): Optional zone-map cache.
            codon_cache (CodonCache | None): Optional shared-codon cache.
        """
        super().__init__(
            ruleset if ruleset is not None else build_fact_ruleset(),
            annotator_context,
            zone_map=zone_map,
            codon_cache=codon_cache,
        )
        self._guard_digest_alphabet: _GuardDigestAlphabet = _GuardDigestAlphabet.build(self.ruleset)

    def _compute_guard_flags(self, fact_context: FactContext) -> dict[str, bool]:
        """Compute the guard atom values in this annotator's guard digest alphabet.

        Args:
            fact_context (FactContext): The overlap to evaluate the guard atoms against.

        Returns:
            dict[str, bool]: Guard atom name -> boolean value.
        """
        return {
            atom_name: bool(get_fact(fact_context, atom_name))
            for atom_name in self._guard_digest_alphabet.atom_names
        }

    def _compute_digest(self, fact_context: FactContext) -> str:
        """Compute digest string from the guard flags for an overlap.

        Args:
            fact_context (FactContext): The overlap to evaluate the guard atoms against.

        Returns:
            str: A digest string suitable as a candidate-rule-list cache key: identical for two
                overlaps whenever every guard atom agrees between them.
        """
        return "".join(
            f"{atom_name}{int(value)}" for atom_name, value in self._compute_guard_flags(fact_context).items()
        )

    def _candidates(self, fact_context: FactContext) -> tuple[FactRule, ...]:
        """Returns the rules whose cheap guard holds for this overlap.

        Args:
            fact_context (FactContext): The context for the overlap being evaluated.

        Returns:
            tuple[FactRule, ...]: Every rule from self.ruleset whose cheap_guard() atoms are
                all true for this overlap.
        """
        annotator_context = self.annotator_context
        candidate_cache = annotator_context.digest_candidate_cache
        digest = self._compute_digest(fact_context)
        cached = candidate_cache.get(digest)
        if cached is None:
            annotator_context.cache_misses += 1
            guard_flags = self._compute_guard_flags(fact_context)
            cached = tuple(
                rule
                for rule in self.ruleset.rules
                if all(guard_flags.get(atom) for atom in rule.cheap_guard())
            )
            candidate_cache[digest] = cached
        else:
            annotator_context.cache_hits += 1
        return cached

    def annotate_feature_overlap(self, overlap: FeatureOverlap) -> list[OverlapConsequence]:
        """Calculate consequences for an overlap.

        Args:
            overlap (FeatureOverlap): The overlap to annotate.

        Returns:
            list[OverlapConsequence]: The consequences that apply, in rank order. Empty if
                nothing matched.

        Raises:
            TypeError: If overlap is not a TranscriptVariationAllele (subclass of FeatureOverlap).
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
        candidates = self._candidates(fact_context)
        hit_terms = [rule.term for rule in candidates if rule.when.evaluate(fact_context)]
        kept_terms = prune(self.ruleset, hit_terms, apply_default=False)
        consequences = [self.ruleset.by_term[term].to_overlap_consequence() for term in kept_terms]

        self.last_trace = None
        if recorder is not None:
            # Fact events first: recomputing the guard flags below goes back through get_fact, and
            # those memo hits would otherwise be counted as reads the real path never made.
            self.last_trace = AnnotationTrace(
                fact_events=reconstruct.memo_fact_events(fact_context),
                rule_events=reconstruct.guarded_rule_events(
                    self.ruleset.rules,
                    {rule.term for rule in candidates},
                    hit_terms,
                    self._compute_guard_flags(fact_context),
                ),
                prune_events=reconstruct.subsumption_events(self.ruleset.by_term, hit_terms),
                emitted_terms=tuple(kept_terms),
                annotator_name=self.name,
                overlap_label=describe_overlap(overlap),
            )
            logger.debug(self.last_trace.format_explanation())

        return consequences

    def annotate_feature_overlaps(self, overlaps: Sequence[FeatureOverlap]) -> list[list[OverlapConsequence]]:
        """Calculate consequences for a batch of overlaps, one at a time.

        Args:
            overlaps (Sequence[FeatureOverlap]): The overlaps to annotate, in order.

        Returns:
            list[list[OverlapConsequence]]: One consequence list per input overlap, same order.
        """
        return [self.annotate_feature_overlap(overlap) for overlap in overlaps]
