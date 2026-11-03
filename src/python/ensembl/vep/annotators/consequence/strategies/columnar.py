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
"""Evaluate consequence rules bit-parallel across an overlap batch.

Each fact is represented by a Python integer whose bit ``i`` contains the value for overlap ``i``.
Rule formulas become bitwise operations, paying Python dispatch once per rule and batch rather than
once per rule and overlap.

Evaluation orders conjuncts by closure cost and carries a live mask of unresolved overlaps. A fact
column is materialised only for that mask, while each plan node guarantees correct bits within its
input mask and zero elsewhere. Because child masks and results only shrink, unpopulated bits cannot
affect a decided result.
"""

from __future__ import annotations

from collections.abc import Sequence
import logging
from typing import TYPE_CHECKING, ClassVar

from attrs import frozen

from ensembl.vep.annotators.consequence.strategies.base import ConsequenceAnnotator
from ensembl.vep.annotators.consequence.facts.base import REGISTRY, get_fact, total_cost
from ensembl.vep.annotators.consequence.context.annotator import AnnotatorContext
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.registry import register_annotator
from ensembl.vep.annotators.consequence.rules.expressions import And, Atom, Expr, Not, Or
from ensembl.vep.annotators.consequence.rules.rulesets import (
    FactRuleSet,
    build_fact_ruleset,
    prune,
)
from ensembl.vep.trace import AnnotationTrace, TraceRecorder, describe_overlap
from ensembl.vep.trace import reconstruct
from ensembl.vep.models.overlap import FeatureOverlap, OverlapConsequence, TranscriptVariationAllele

if TYPE_CHECKING:
    from ensembl.vep.caches.codon_cache import CodonCache
    from ensembl.vep.caches.zone_map import ZoneMapCache

__all__ = [
    "BitColumnBuffer",
    "ColumnarFactsAnnotator",
    "DEFAULT_BATCH_SIZE",
    "PlanNode",
    "RulePlan",
    "RulePlanSet",
]

logger = logging.getLogger(__name__)

DEFAULT_BATCH_SIZE = 256
"""Default maximum number of overlaps evaluated together as one bit-parallel batch. Increase to amortise the
plan-walking overhead over more overlaps, decrease to reduce peak memory usage. It has been observed that
increasing does not always improve speed, because of big-int arithmetic overhead."""


class BitColumnBuffer:
    """Materialise fact bit-columns lazily for one overlap batch."""

    def __init__(self, fact_contexts: Sequence[FactContext]) -> None:
        """Initialize the buffer.

        Args:
            fact_contexts (Sequence[FactContext]): The overlaps' fact contexts, batch-index order.
        """
        self._fact_contexts = fact_contexts
        self._known_bits: dict[str, int] = {}
        self._value_bits: dict[str, int] = {}

    def seed(self, fact_name: str, column: int, known_mask: int) -> None:
        """Mark a fact's column as already computed for every context in ``known_mask``.

        For a caller that produced the whole column by cheaper means than calling ``get_fact``
        per context. Unlike FactContext.prefill this does not touch the per-overlap memo, so a fact
        seeded here is visible to the rule loop but not to a ``get_fact`` call inside another fact's
        body. Therefore only seed leaf facts nothing else reads, or prefill the memo as well.

        Args:
            fact_name (str): The fact atom being seeded.
            column (int): Bit ``i`` set where the fact is true for context ``i``.
            known_mask (int): The contexts this column is exact for.
        """
        self._value_bits[fact_name] = column & known_mask
        self._known_bits[fact_name] = self._known_bits.get(fact_name, 0) | known_mask

    def column(self, fact_name: str, live_mask: int) -> int:
        """Return the bit-column for ``fact_name`` within ``live_mask``.

        Args:
            fact_name (str): The fact atom to materialise a column for.
            live_mask (int): The set of overlap indices (as bits) the caller needs an answer for.

        Returns:
            int: Bits set for every overlap index in ``live_mask`` where the fact is true. Bits
                outside ``live_mask`` are not guaranteed to mean anything.
        """
        missing_mask = live_mask & ~self._known_bits.get(fact_name, 0)
        if missing_mask:
            value_bits = self._value_bits.get(fact_name, 0)
            remaining_mask = missing_mask
            while remaining_mask:
                lowest_bit = remaining_mask & -remaining_mask
                context_index = lowest_bit.bit_length() - 1
                if get_fact(self._fact_contexts[context_index], fact_name):
                    value_bits |= lowest_bit
                remaining_mask ^= lowest_bit
            self._value_bits[fact_name] = value_bits
            self._known_bits[fact_name] = self._known_bits.get(fact_name, 0) | missing_mask
        return self._value_bits.get(fact_name, 0)


@frozen(kw_only=True, auto_attribs=True)
class _AtomPlanNode:
    """A compiled plan node for a single fact atom."""

    fact_name: str

    def evaluate(self, live_mask: int, bit_columns: BitColumnBuffer) -> int:
        """Evaluate this atom over ``live_mask``.

        Args:
            live_mask (int): The overlap indices (as bits) still undecided for this subtree.
            bit_columns (BitColumnBuffer): The batch's lazily-materialised fact columns.

        Returns:
            int: Bits set where the fact is true, within ``live_mask``.
        """
        return bit_columns.column(self.fact_name, live_mask) & live_mask


@frozen(kw_only=True, auto_attribs=True)
class _NotPlanNode:
    """A compiled plan node for a negation."""

    child: PlanNode

    def evaluate(self, live_mask: int, bit_columns: BitColumnBuffer) -> int:
        """Evaluate this negation over ``live_mask``.

        Args:
            live_mask (int): The overlap indices (as bits) still undecided for this subtree.
            bit_columns (BitColumnBuffer): The batch's lazily-materialised fact columns.

        Returns:
            int: Bits set where the child is false, within ``live_mask``.
        """
        return live_mask & ~self.child.evaluate(live_mask, bit_columns)


@frozen(kw_only=True, auto_attribs=True)
class _AndPlanNode:
    """A compiled plan node for a conjunction, children ordered cheapest-first."""

    children: tuple[PlanNode, ...]

    def evaluate(self, live_mask: int, bit_columns: BitColumnBuffer) -> int:
        """Evaluate this conjunction, narrowing the mask as each child settles.

        Args:
            live_mask (int): The overlap indices (as bits) still undecided for this subtree.
            bit_columns (BitColumnBuffer): The batch's lazily-materialised fact columns.

        Returns:
            int: Bits set where every child is true, within ``live_mask``.
        """
        still_live = live_mask
        for child in self.children:
            if not still_live:
                break
            still_live &= child.evaluate(still_live, bit_columns)
        return still_live


@frozen(kw_only=True, auto_attribs=True)
class _OrPlanNode:
    """A compiled plan node for a disjunction, children ordered cheapest-first."""

    children: tuple[PlanNode, ...]

    def evaluate(self, live_mask: int, bit_columns: BitColumnBuffer) -> int:
        """Evaluate this disjunction over ``live_mask``.

        Tracks what is still *undecided* (neither proven True nor exhausted), the dual of the
        And case, so each child is only asked about the indices no earlier child has settled.

        Args:
            live_mask (int): The overlap indices (as bits) still undecided for this subtree.
            bit_columns (BitColumnBuffer): The batch's lazily-materialised fact columns.

        Returns:
            int: Bits set where any child is true, within ``live_mask``.
        """
        undecided_mask = live_mask
        for child in self.children:
            if not undecided_mask:
                break
            undecided_mask &= ~child.evaluate(undecided_mask, bit_columns)
        return live_mask & ~undecided_mask


PlanNode = _AtomPlanNode | _NotPlanNode | _AndPlanNode | _OrPlanNode


@frozen(kw_only=True, auto_attribs=True)
class RulePlan:
    """Represent a rule as a cost-ordered bit-parallel evaluation plan.

    Attributes:
        root_plan_node (PlanNode): The root of the compiled plan.
    """

    root_plan_node: PlanNode
    _transitive_closure_cache: ClassVar[dict[str, frozenset[str]]] = {}
    """Shared across every RulePlan, deliberately, since the fact dependency graph is a property of
    REGISTRY, not of any one rule. A per-instance cache would rebuild the same closures once per rule.
    """

    @classmethod
    def build(cls, expression: Expr) -> RulePlan:
        """Compile a rule's Expr formula into a RulePlan.

        Args:
            expression (Expr): The formula to compile.

        Returns:
            RulePlan: The compiled plan, with every ``And``/``Or`` node's children sorted by
                ascending cost so the cheapest conjunct/disjunct gets the chance to settle the
                whole subtree first.
        """
        return cls(root_plan_node=cls._compile(expression))

    def evaluate(self, live_mask: int, bit_columns: BitColumnBuffer) -> int:
        """Evaluate this plan over a live mask, short-circuiting as children settle.

        Args:
            live_mask (int): The set of overlap indices (as bits) still undecided.
            bit_columns (BitColumnBuffer): The batch's lazily-materialised fact columns.

        Returns:
            int: Bits set for every overlap index in ``live_mask`` where the rule evaluates True.
        """
        return self.root_plan_node.evaluate(live_mask, bit_columns)

    @classmethod
    def _compile(cls, expression: Expr) -> PlanNode:
        """Recursively compile ``expression`` into cost-ordered plan nodes.

        Args:
            expression (Expr): The (sub)expression to compile.

        Returns:
            PlanNode: The compiled node.

        Raises:
            TypeError: If ``expression`` is a node type this annotator does not support.
        """
        if isinstance(expression, Atom):
            return _AtomPlanNode(fact_name=expression.name)
        if isinstance(expression, Not):
            return _NotPlanNode(child=cls._compile(expression.child))
        if isinstance(expression, (And, Or)):
            ordered_children = sorted(expression.children, key=cls._expression_cost)
            compiled_children = tuple(cls._compile(child) for child in ordered_children)
            return (
                _AndPlanNode(children=compiled_children)
                if isinstance(expression, And)
                else _OrPlanNode(children=compiled_children)
            )
        raise TypeError(f"unsupported expression node type: {type(expression).__name__}")

    @classmethod
    def _expression_cost(cls, expression: Expr) -> int:
        """Returns the cost of deciding ``expression``, counting each shared fact dependency once.

        Args:
            expression (Expr): The expression (typically one conjunct/disjunct of a rule's
                formula) to cost.

        Returns:
            int: The total cost of every fact ``expression`` transitively depends on.
        """
        reachable_fact_names: set[str] = set()
        for atom_name in expression.atoms():
            reachable_fact_names |= cls._transitive_closure(atom_name)
        return total_cost(reachable_fact_names)

    @classmethod
    def _transitive_closure(cls, fact_name: str) -> frozenset[str]:
        """Returns ``fact_name`` plus every fact its declared ``deps`` will pull in, transitively.

        Args:
            fact_name (str): The fact to start from.

        Returns:
            frozenset[str]: ``fact_name`` and every fact reachable from it via ``deps``.
        """
        cache = cls._transitive_closure_cache
        cached = cache.get(fact_name)
        if cached is not None:
            return cached
        cache[fact_name] = frozenset({fact_name})  # break any cycle in the declarations
        reachable = {fact_name}
        for dependency_name in REGISTRY[fact_name].deps:
            reachable |= cls._transitive_closure(dependency_name)
        result = frozenset(reachable)
        cache[fact_name] = result
        return result


@frozen(kw_only=True, auto_attribs=True)
class RulePlanSet:
    """Pair a rule set with one compiled plan per rule.

    Attributes:
        ruleset (FactRuleSet): The compiled rule set, in rank order.
        plans (tuple[RulePlan, ...]): One compiled plan per rule in ``ruleset.rules``, same order.
    """

    ruleset: FactRuleSet
    plans: tuple[RulePlan, ...]

    @classmethod
    def build(cls, ruleset: FactRuleSet) -> RulePlanSet:
        """Compile one plan per rule in ``ruleset``.

        Args:
            ruleset (FactRuleSet): The rule set to compile plans for.

        Returns:
            RulePlanSet: The rule set paired with its compiled plans.
        """
        return cls(ruleset=ruleset, plans=tuple(RulePlan.build(rule.when) for rule in ruleset.rules))

    def hit_terms_by_context(self, bit_columns: BitColumnBuffer, context_count: int) -> list[list[str]]:
        """Return the terms that evaluated as true for each fact context.

        Args:
            bit_columns (BitColumnBuffer): The batch's fact columns.
            context_count (int): Number of overlaps in this chunk.

        Returns:
            list[list[str]]: The terms that fired for each fact context, in rule rank order.
        """
        hit_terms_by_context: list[list[str]] = [[] for _ in range(context_count)]
        full_mask = (1 << context_count) - 1
        for rule, plan in zip(self.ruleset.rules, self.plans):
            fired_mask = plan.evaluate(full_mask, bit_columns) & full_mask
            remaining_mask = fired_mask
            while remaining_mask:
                lowest_bit = remaining_mask & -remaining_mask
                context_index = lowest_bit.bit_length() - 1
                hit_terms_by_context[context_index].append(rule.term)
                remaining_mask ^= lowest_bit
        return hit_terms_by_context

    def evaluate_over_columns(
        self, bit_columns: BitColumnBuffer, context_count: int
    ) -> list[list[OverlapConsequence]]:
        """Evaluate every rule bit-parallel over a column buffer according to their plan.

        The column buffer here, already built, gets populated via plain lazy materialisation.

        Args:
            bit_columns (BitColumnBuffer): The batch's fact columns.
            context_count (int): Number of overlaps in this chunk.

        Returns:
            list[list[OverlapConsequence]]: One consequence list per input context, same order.
        """
        ruleset = self.ruleset
        return [
            [
                ruleset.by_term[term].to_overlap_consequence()
                for term in prune(ruleset, hit_terms, apply_default=False)
            ]
            for hit_terms in self.hit_terms_by_context(bit_columns, context_count)
        ]


@register_annotator("consequence.columnar")
class ColumnarFactsAnnotator(ConsequenceAnnotator[FactRuleSet]):
    """Evaluate consequence rules bit-parallel over overlap batches.

    Attributes:
        name (str): Registry name, "consequence.columnar".
        ruleset (FactRuleSet): The compiled rule set this annotator evaluates.
        batch_size (int): Maximum number of overlaps evaluated together as one bit-parallel batch.
        last_traces (list[AnnotationTrace]): One derivation per overlap in the most recent
            annotate_feature_overlaps call, in the same order, or an empty list if that call did
            not trace. The inherited last_trace attribute mirrors this list's final entry.
    """

    name = "consequence.columnar"

    def __init__(
        self,
        ruleset: FactRuleSet | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
        annotator_context: AnnotatorContext | None = None,
        *,
        zone_map: ZoneMapCache | None = None,
        codon_cache: CodonCache | None = None,
    ) -> None:
        """Initialize the annotator.

        Args:
            ruleset (FactRuleSet | None): The rule set to evaluate, or None to build the default
                rule set from the checked-in consequence catalogue.
            batch_size (int): Maximum number of overlaps evaluated together as one batch.
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
        self.batch_size = batch_size
        self._planset = RulePlanSet.build(self.ruleset)
        self.last_traces: list[AnnotationTrace] = []

    def _evaluate_chunk(self, fact_contexts: Sequence[FactContext]) -> list[list[OverlapConsequence]]:
        """Evaluate every rule, bit-parallel, over one chunk of at most ``batch_size`` overlaps.

        Args:
            fact_contexts (Sequence[FactContext]): The chunk's fact contexts, batch-index order.

        Returns:
            list[list[OverlapConsequence]]: One consequence list per input context, same order.
        """
        context_count = len(fact_contexts)
        if context_count == 0:
            return []
        bit_columns = BitColumnBuffer(fact_contexts)
        return self._planset.evaluate_over_columns(bit_columns, context_count)

    def _trace_chunk(self, fact_contexts: Sequence[FactContext]) -> list[AnnotationTrace]:
        """Debug-only replay of one already-evaluated chunk, to explain each overlap.

        Deliberately a second and independent pass over the plan rather than a change to _evaluate_chunk which
        is the hot batch path. Only ever called when tracing is already known to be enabled.

        Args:
            fact_contexts (Sequence[FactContext]): The already-evaluated chunk, in the same order
                passed to the preceding ``_evaluate_chunk`` call.

        Returns:
            list[AnnotationTrace]: One derivation per context, in the same order.
        """
        context_count = len(fact_contexts)
        if context_count == 0:
            return []
        # Taken before the replay: every column below is a memo hit, and get_fact counts those,
        # so read counts gathered afterwards would include reads the real pass never made.
        fact_events_by_context = [reconstruct.memo_fact_events(context) for context in fact_contexts]

        bit_columns = BitColumnBuffer(fact_contexts)
        hit_terms_by_context = self._planset.hit_terms_by_context(bit_columns, context_count)

        traces = []
        for context_index, fact_context in enumerate(fact_contexts):
            hit_terms = hit_terms_by_context[context_index]
            trace = AnnotationTrace(
                fact_events=fact_events_by_context[context_index],
                rule_events=reconstruct.plain_rule_events(self.ruleset.rules, hit_terms),
                prune_events=reconstruct.subsumption_events(self.ruleset.by_term, hit_terms),
                emitted_terms=tuple(prune(self.ruleset, hit_terms, apply_default=False)),
                annotator_name=self.name,
                overlap_label=describe_overlap(fact_context.overlap),
            )
            traces.append(trace)
            logger.debug(trace.format_explanation())
        return traces

    def annotate_feature_overlaps(self, overlaps: Sequence[FeatureOverlap]) -> list[list[OverlapConsequence]]:
        """Calculate consequences for a batch of overlaps at once.

        Args:
            overlaps (Sequence[FeatureOverlap]): The overlaps to annotate, in order.

        Returns:
            list[list[OverlapConsequence]]: One consequence list per input overlap, same order.

        Raises:
            TypeError: If any element of ``overlaps`` is not a TranscriptVariationAllele.
        """
        transcript_variation_alleles: list[TranscriptVariationAllele] = []
        for overlap in overlaps:
            if not isinstance(overlap, TranscriptVariationAllele):
                raise TypeError(
                    f"{type(self).__name__} only annotates TranscriptVariationAllele overlaps, got "
                    f"{type(overlap).__name__}"
                )
            transcript_variation_alleles.append(overlap)

        debug_enabled = logger.isEnabledFor(logging.DEBUG)
        fact_contexts = [
            FactContext(
                overlap=overlap,
                recorder=(TraceRecorder() if debug_enabled else None),
                annotator_context=self.annotator_context,
            )
            for overlap in transcript_variation_alleles
        ]
        results: list[list[OverlapConsequence]] = []
        traces: list[AnnotationTrace] = []
        for chunk_start in range(0, len(fact_contexts), self.batch_size):
            chunk = fact_contexts[chunk_start : chunk_start + self.batch_size]
            # Caches are applied per overlap, before the chunk's columns are materialised, because
            # prefilling after expensive tier gets called would save nothing.
            for fact_context in chunk:
                self.annotator_context.prepare(fact_context)
            results.extend(self._evaluate_chunk(chunk))
            if debug_enabled:
                traces.extend(self._trace_chunk(chunk))

        self.last_traces = traces
        self.last_trace = traces[-1] if traces else None
        return results
