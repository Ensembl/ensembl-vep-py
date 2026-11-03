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
"""Columnar-staged evaluation whose cheapest fact columns are filled by native batch kernels.

Only flat integer arrays cross the native boundary; model objects remain in Python. Kernels process
whole columns because a scalar boundary crossing costs more than the operations being accelerated.

The seeded facts cover short-variant shape and feature-span comparisons that remain per-context
work after zone-map prefill. Exon and intron kernels are not seeded here because building their CSR
arrays costs more than the Python scan when a zone map has not already supplied the values.
"""

from __future__ import annotations

from collections.abc import Sequence
import logging
from typing import TYPE_CHECKING

from ensembl.vep.annotators.registry import register_annotator
from ensembl.vep.annotators.consequence.strategies.columnar import DEFAULT_BATCH_SIZE
from ensembl.vep.annotators.consequence.strategies.columnar_staged import ColumnarStagedAnnotator
from ensembl.vep.annotators.consequence.context.annotator import AnnotatorContext
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.rules.rulesets import FactRuleSet
from ensembl.vep.common.constants import EMPTY_ALLELE
from ensembl.vep.models.variation import StructuralVariationAllele
from ensembl.vep.ensembl_vep_core import (  # pylint: disable=no-name-in-module
    batch_contain_span,
    batch_overlap_span,
    batch_variant_shapes,
)

if TYPE_CHECKING:
    from ensembl.vep.caches.codon_cache import CodonCache
    from ensembl.vep.caches.zone_map import ZoneMapCache

__all__ = ["ColumnarKernelsAnnotator", "SEEDED_FACT_NAMES"]

logger = logging.getLogger(__name__)

# The short-variant shape facts, paired with the index of the kernel mask that answers each.
_SHORT_SHAPE_FACTS: tuple[tuple[str, int], ...] = (
    ("is_snv", 0),
    ("is_short_insertion", 1),
    ("is_short_deletion", 2),
)

SEEDED_FACT_NAMES: tuple[str, ...] = tuple(name for name, _ in _SHORT_SHAPE_FACTS) + (
    "is_short",
    "overlap_feature",
    "contain_feature",
)
"""Facts filled by native kernels and checked against their Python implementations in tests."""


@register_annotator("consequence.columnar-kernels")
class ColumnarKernelsAnnotator(ColumnarStagedAnnotator):
    """Evaluate staged columns with native kernels for inexpensive facts.

    Identical output to ColumnarStagedAnnotator by construction: the kernels compute the same
    formulas, and everything they do not seed falls through to the ordinary Python path.

    Attributes:
        name (str): Registry name, "consequence.columnar-kernels".
        kernel_calls (int): Number of FFI crossings made so far.
        columns_seeded (int): Number of (fact, overlap) values filled by a kernel so far.
    """

    name = "consequence.columnar-kernels"

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
            ruleset (FactRuleSet | None): The rule set to evaluate, or None for the default.
            batch_size (int): Maximum number of overlaps evaluated together as one batch.
            annotator_context (AnnotatorContext | None): Run-scoped caches to evaluate against,
                or None (the default) for a fresh, private one.
            zone_map (ZoneMapCache | None): Optional zone-map cache.
            codon_cache (CodonCache | None): Optional shared-codon cache.
        """
        super().__init__(
            ruleset,
            batch_size,
            annotator_context,
            zone_map=zone_map,
            codon_cache=codon_cache,
        )
        self.kernel_calls = 0
        self.columns_seeded = 0

    def seed_columns(self, fact_contexts: Sequence[FactContext]) -> dict[str, int]:
        """Compute every kernel-backed fact for the whole chunk, in a couple of FFI crossings.

        The short-variant shape facts are returned as columns only and not prefilled into each context's
        memo as nothing reads them through ``get_fact``. The overlap_feature and contain_feature are
        different, as every positional fact declares ``requires=("overlap_feature",)``, those two are
        additionally prefilled into each memo. This per-context fan-out is the single most expensive part of
        this strategy, which is why it is confined to the two facts that genuinely need it.

        Args:
            fact_contexts (Sequence[FactContext]): The chunk's contexts, batch-index order.

        Returns:
            dict[str, int]: Fact name -> bit column for the whole chunk.
        """
        context_count = len(fact_contexts)
        if not context_count:
            return {}

        variant_mins: list[int] = []
        variant_maxs: list[int] = []
        feature_starts: list[int] = []
        feature_ends: list[int] = []
        reference_lengths: list[int] = []
        alternate_lengths: list[int] = []
        short_column = 0
        for context_index, fact_context in enumerate(fact_contexts):
            variant_min, variant_max = fact_context.variant_span
            variant_mins.append(variant_min)
            variant_maxs.append(variant_max)
            feature = fact_context.overlap.feature
            feature_starts.append(feature.start)
            feature_ends.append(feature.end)
            variant = fact_context.overlap.overlapping_feature
            reference_lengths.append(len((variant.reference_allele or "").replace(EMPTY_ALLELE, "")))
            alternate_lengths.append(len((variant.alternate_allele or "").replace(EMPTY_ALLELE, "")))
            if not isinstance(variant, StructuralVariationAllele):
                short_column |= 1 << context_index

        shape_masks = batch_variant_shapes(reference_lengths, alternate_lengths)
        overlap_bitmap = batch_overlap_span(variant_mins, variant_maxs, feature_starts, feature_ends)
        contain_bitmap = batch_contain_span(variant_mins, variant_maxs, feature_starts, feature_ends)
        self.kernel_calls += 3

        columns: dict[str, int] = {
            fact_name: int.from_bytes(shape_masks[mask_index], "little") & short_column
            for fact_name, mask_index in _SHORT_SHAPE_FACTS
        }
        # The kernel answers the length comparison for every overlap in the chunk, including any
        # structural ones. `requires=("is_short",)` is what makes those answers False in the
        # ordinary path, but a seeded column bypasses the decorator's guard entirely. So the gate
        # has to be reapplied here. Without it a structural allele whose ref/alt strings happen to
        # differ in length might be seeded as a short insertion and the kernel path would disagree
        # with every other annotator.
        columns["is_short"] = short_column

        overlap_column = int.from_bytes(overlap_bitmap, "little")
        contain_column = int.from_bytes(contain_bitmap, "little")
        columns["contain_feature"] = contain_column

        # We also add in the zone map columns, if a zone map is attached.
        zone_map_cache = self.annotator_context.zone_map_cache
        if zone_map_cache is not None:
            columns.update(zone_map_cache.columns(fact_contexts))

        # When a zone map is attached it has already answered the `overlap_feature` fact, so the kernel's
        # answer is redundant.
        if zone_map_cache is None:
            columns["overlap_feature"] = overlap_column
            for context_index, fact_context in enumerate(fact_contexts):
                bit = 1 << context_index
                fact_context.prefill({"overlap_feature": bool(overlap_column & bit)})
        self.columns_seeded += context_count * len(columns)
        return columns
