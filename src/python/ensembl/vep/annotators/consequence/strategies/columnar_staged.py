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
"""Combine bit-parallel rule evaluation with staged fact propagation.

Cheap fact columns are materialised for the complete chunk. A topological sweep of the acyclic
``requires`` graph then propagates false bits into forced-false masks without per-context fixpoint
work. This matches ``implied_false`` while preserving batch evaluation.

Forced-false masks may derive only from full-chunk columns or earlier forced masks. Values from a
column populated for a smaller live mask cannot establish facts outside that mask.
"""

from __future__ import annotations

from collections.abc import Sequence
import logging
from typing import TYPE_CHECKING

from ensembl.vep.annotators.consequence.strategies.base import ConsequenceAnnotator
from ensembl.vep.annotators.consequence.facts.base import REGISTRY
from ensembl.vep.annotators.consequence.context.annotator import AnnotatorContext
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.registry import register_annotator
from ensembl.vep.annotators.consequence.strategies.columnar import (
    DEFAULT_BATCH_SIZE,
    BitColumnBuffer,
    RulePlanSet,
)
from ensembl.vep.annotators.consequence.strategies.staged import CheapAtomAlphabet
from ensembl.vep.annotators.consequence.rules.rulesets import FactRuleSet, build_fact_ruleset, prune
from ensembl.vep.trace import AnnotationTrace, TraceRecorder, describe_overlap
from ensembl.vep.trace import reconstruct
from ensembl.vep.models.overlap import FeatureOverlap, OverlapConsequence, TranscriptVariationAllele

if TYPE_CHECKING:
    from ensembl.vep.caches.codon_cache import CodonCache
    from ensembl.vep.caches.zone_map import ZoneMapCache

__all__ = ["ColumnarStagedAnnotator", "topological_requires_order"]

logger = logging.getLogger(__name__)


def topological_requires_order() -> list[tuple[str, tuple[str, ...]]]:
    """Return facts with requirements in topological order.

    Ordered so every atom in a fact's ``requires`` precedes the fact itself, an explicit
    topological sort over the ``requires`` graph in REGISTRY, not a reliance on declaration order.

    Returns:
        list[tuple[str, tuple[str, ...]]]: Topologically ordered ``(fact_name, requires)`` pairs.

    Raises:
        ValueError: If the ``requires`` graph has a cycle. A cycle or a bad forward reference
            must fail loudly at annotator construction, not silently under-propagate at
            evaluation time.
    """
    gated: dict[str, tuple[str, ...]] = {
        name: registered_fact.requires
        for name, registered_fact in REGISTRY.items()
        if registered_fact.requires
    }
    order: list[str] = []
    permanently_visited: set[str] = set()
    in_progress: set[str] = set()

    def visit(fact_name: str) -> None:
        if fact_name in permanently_visited or fact_name not in gated:
            return
        if fact_name in in_progress:
            raise ValueError(f"cycle in `requires` graph at {fact_name!r}")
        in_progress.add(fact_name)
        for required_name in gated[fact_name]:
            visit(required_name)
        in_progress.discard(fact_name)
        permanently_visited.add(fact_name)
        order.append(fact_name)

    for fact_name in gated:
        visit(fact_name)
    return [(fact_name, gated[fact_name]) for fact_name in order]


class _StagedBitColumnBuffer(BitColumnBuffer):
    """BitColumnBuffer plus forced-false masks.

    The function set_forced_false must be called once, with masks derived only from populated columns
    across full-chunk and before any column call that should benefit from them. The override marks forced
    bits known, so that the base class's ``get_fact`` loop skips them entirely, without ever setting them in
    ``_value_bits`` (they are correctly False, i.e. absent, by construction).
    """

    def __init__(self, fact_contexts: Sequence[FactContext]) -> None:
        """Initialize the buffer.

        Args:
            fact_contexts (Sequence[FactContext]): The overlaps' fact contexts, batch-index order.
        """
        super().__init__(fact_contexts)
        self._forced_false: dict[str, int] = {}

    def set_forced_false(self, forced: dict[str, int]) -> None:
        """Record the per-fact forced-false masks derived by the topological sweep.

        Args:
            forced (dict[str, int]): Fact name -> mask of overlap indices forced False by
                ``requires`` propagation.
        """
        self._forced_false = forced

    def column(self, fact_name: str, live_mask: int) -> int:
        """Return the bit-column for ``fact_name``, honouring forced-false masks.

        Args:
            fact_name (str): The fact atom to materialise a column for.
            live_mask (int): The set of overlap indices (as bits) the caller needs an answer for.

        Returns:
            int: Bits set for every overlap index in ``live_mask`` where the fact is true. Bits
                outside ``live_mask`` are not guaranteed to mean anything.
        """
        forced_false_mask = self._forced_false.get(fact_name)
        if forced_false_mask:
            newly_known_mask = live_mask & forced_false_mask & ~self._known_bits.get(fact_name, 0)
            if newly_known_mask:
                self._known_bits[fact_name] = self._known_bits.get(fact_name, 0) | newly_known_mask
        return super().column(fact_name, live_mask)


@register_annotator("consequence.columnar-staged")
class ColumnarStagedAnnotator(ConsequenceAnnotator[FactRuleSet]):
    """Evaluate bit-parallel rules while avoiding unnecessary expensive facts.

    Attributes:
        name (str): Registry name, "consequence.columnar-staged".
        ruleset (FactRuleSet): The compiled rule set this annotator evaluates.
        batch_size (int): Maximum number of overlaps evaluated together as one bit-parallel batch.
        last_traces (list[AnnotationTrace]): One derivation per overlap in the most recent
            annotate_feature_overlaps call, in the same order, or an empty list if that call did
            not trace. The inherited last_trace attribute mirrors this list's final entry.
    """

    name = "consequence.columnar-staged"

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
        self._cheap_atom_alphabet = CheapAtomAlphabet.build(self.ruleset.rules)
        # Computed and validated once at construction, not per chunk.
        self._topological_requires_order = topological_requires_order()
        self.last_traces: list[AnnotationTrace] = []

    def seed_columns(self, fact_contexts: Sequence[FactContext]) -> dict[str, int]:
        """Return fact columns already known for the whole chunk.

        Extension point for a subclass that can compute a fact for every context in the chunk at
        once by some cheaper route than per-context ``get_fact`` (see ColumnarKernelsAnnotator).

        A returned column must be exact for **every** context in the chunk, since it is marked
        fully known. A partially-correct column would be read as authoritative and could lead to
        incorrect results. Returning nothing (the default) leaves every column to the ordinary
        lazy path.

        Args:
            fact_contexts (Sequence[FactContext]): The chunk's contexts, batch-index order.

        Returns:
            dict[str, int]: Fact name -> bit column, bit ``i`` being the value for context ``i``.
        """
        zone_map_cache = self.annotator_context.zone_map_cache
        if zone_map_cache is None:
            return {}
        # Evaluated and returned as columns rather than per-overlap memo prefills
        return zone_map_cache.columns(fact_contexts)

    def _evaluate_staged_chunk(self, fact_contexts: Sequence[FactContext]) -> list[list[OverlapConsequence]]:
        """Evaluate every rule, bit-parallel and ``requires``-staged, over one chunk of overlaps.

        Args:
            fact_contexts (Sequence[FactContext]): The chunk's fact contexts, batch-index order.

        Returns:
            list[list[OverlapConsequence]]: One consequence list per input context, same order.
        """
        context_count = len(fact_contexts)
        if context_count == 0:
            return []
        full_mask = (1 << context_count) - 1
        bit_columns = _StagedBitColumnBuffer(fact_contexts)

        # 0. Anything a subclass can fill for the whole chunk at once, marked fully known so the
        #    lazy path never touches it.
        for fact_name, column in self.seed_columns(fact_contexts).items():
            bit_columns.seed(fact_name, column, full_mask)

        # 1. Materialise every cheap atom for the whole chunk. The same upfront spend
        #    consequence.staged already pays, memoised per context. The seeded columns above may
        #    do not get re-computed.
        known_true_masks = {
            atom_name: bit_columns.column(atom_name, full_mask)
            for atom_name in self._cheap_atom_alphabet.atom_names
        }
        false_bits = {atom_name: full_mask & ~mask for atom_name, mask in known_true_masks.items()}

        # 2. One topological sweep of big-int ORs. A bit-parallel version of `implied_false`.
        forced: dict[str, int] = {}
        for fact_name, requires_names in self._topological_requires_order:
            forced_mask = 0
            for required_name in requires_names:
                forced_mask |= false_bits.get(required_name, 0)
            forced_mask &= full_mask
            if fact_name in false_bits:
                # Already a directly-materialised cheap column (e.g. overlap_exon, which itself has
                # requires=("overlap_feature",)): its column is already exact everywhere, so
                # this OR is a no-op that just keeps the sweep self-consistent for anything gated
                # behind it.
                false_bits[fact_name] |= forced_mask
                continue
            if forced_mask:
                forced[fact_name] = forced_mask
                false_bits[fact_name] = forced_mask  # feeds facts gated on this one
        bit_columns.set_forced_false(forced)

        # 3. Now evaluate each rule using already pre-filled columns.
        return self._planset.evaluate_over_columns(bit_columns, context_count)

    def _trace_staged_chunk(self, fact_contexts: Sequence[FactContext]) -> list[AnnotationTrace]:
        """Debug-only replay of one already-evaluated chunk, to explain each overlap.

        Deliberately a second and independent pass over the plan rather than a change to
        _evaluate_staged_chunk method which is the hot batch path. Only ever called when tracing is
        already known to be enabled.

        Every ``bit_columns.column`` call below is a memo hit -- the real pass already computed and
        charged everything this chunk could possibly need -- so this replay never charges any
        additional cost. What it rebuilds is the per-context partial assignment; the forced-False
        *facts* fall out of that assignment afterwards, so nothing here re-implements the
        propagation in order to explain it.

        Args:
            fact_contexts (Sequence[FactContext]): The already-evaluated chunk, in the same order
                passed to the preceding ``_evaluate_staged_chunk`` call.

        Returns:
            list[AnnotationTrace]: One derivation per context, in the same order.
        """
        context_count = len(fact_contexts)
        if context_count == 0:
            return []
        # Taken before the replay: every column below is a memo hit, and get_fact counts those, so
        # read counts gathered afterwards would include reads the real pass never made.
        fact_events_by_context = [reconstruct.memo_fact_events(context) for context in fact_contexts]

        full_mask = (1 << context_count) - 1
        bit_columns = _StagedBitColumnBuffer(fact_contexts)

        known_true_masks = {
            atom_name: bit_columns.column(atom_name, full_mask)
            for atom_name in self._cheap_atom_alphabet.atom_names
        }
        false_bits = {atom_name: full_mask & ~mask for atom_name, mask in known_true_masks.items()}

        forced: dict[str, int] = {}
        for fact_name, requires_names in self._topological_requires_order:
            forced_mask = 0
            for required_name in requires_names:
                forced_mask |= false_bits.get(required_name, 0)
            forced_mask &= full_mask
            if fact_name in false_bits:
                false_bits[fact_name] |= forced_mask
                continue
            if not forced_mask:
                continue
            forced[fact_name] = forced_mask
            false_bits[fact_name] = forced_mask
        bit_columns.set_forced_false(forced)

        hit_terms_by_context = self._planset.hit_terms_by_context(bit_columns, context_count)

        traces = []
        for context_index, fact_context in enumerate(fact_contexts):
            bit = 1 << context_index
            known = {atom_name: bool(mask & bit) for atom_name, mask in known_true_masks.items()}
            for atom_name, mask in forced.items():
                if mask & bit:
                    known[atom_name] = False
            hit_terms = hit_terms_by_context[context_index]
            trace = AnnotationTrace(
                fact_events=fact_events_by_context[context_index]
                + reconstruct.forced_false_events(known, fact_context.memo),
                rule_events=reconstruct.partial_rule_events(self.ruleset.rules, known, hit_terms),
                prune_events=reconstruct.subsumption_events(self.ruleset.by_term, hit_terms),
                emitted_terms=tuple(prune(self.ruleset, hit_terms, apply_default=False)),
                annotator_name=self.name,
                overlap_label=describe_overlap(fact_context.overlap),
            )
            traces.append(trace)
            logger.debug(trace.format_explanation())
        return traces

    def annotate_feature_overlaps(self, overlaps: Sequence[FeatureOverlap]) -> list[list[OverlapConsequence]]:
        """Calculate consequences for a batch of transcript/variant-allele overlaps at once.

        Args:
            overlaps (Sequence[FeatureOverlap]): The overlaps to annotate, in order. Every element
                must be a TranscriptVariationAllele.

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

        # See ColumnarFactsAnnotator.annotate_feature_overlaps: a single isEnabledFor check gates
        # everything tracing-related below.
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
            # Caches are applied per overlap, before the chunk's columns are materialised --
            # prefilling after the sequence tier has been reached would save nothing.
            for fact_context in chunk:
                # The zone map is applied per chunk as columns instead -- see seed_columns.
                self.annotator_context.prepare(fact_context, apply_zone_map=False)
            results.extend(self._evaluate_staged_chunk(chunk))
            if debug_enabled:
                traces.extend(self._trace_staged_chunk(chunk))

        self.last_traces = traces
        self.last_trace = traces[-1] if traces else None
        return results
