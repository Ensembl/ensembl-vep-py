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
"""Shared consequence annotator implementation."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Generic, Sequence, TypeVar

from ensembl.vep.annotators.consequence.context import AnnotatorContext
from ensembl.vep.common import Readthrough
from ensembl.vep.models.overlap import FeatureOverlap, OverlapConsequence

if TYPE_CHECKING:
    from ensembl.vep.trace import AnnotationTrace
    from ensembl.vep.caches.codon_cache import CodonCache
    from ensembl.vep.caches.zone_map import ZoneMapCache

__all__ = ["ConsequenceAnnotator", "RulesetT"]

RulesetT = TypeVar("RulesetT")


class ConsequenceAnnotator(ABC, Generic[RulesetT]):
    """Implement shared behaviour for compiled consequence-rule annotators.

    Attributes:
        name (str): The registry name this annotator is registered under.
        ruleset (RulesetT): The compiled rule set this annotator evaluates.
        last_trace (AnnotationTrace | None): The derivation trace recorded by the most recent
            ``annotate_feature_overlap`` call, or None if that call did not trace (the annotators
            did not opt-in). The native ``rust.*`` annotators never set this.
        annotator_context (AnnotatorContext): The run-scoped caches this annotator evaluates
            against.
    """

    name: str = "base"

    def __init__(
        self,
        ruleset: RulesetT,
        annotator_context: AnnotatorContext | None = None,
        *,
        zone_map: ZoneMapCache | None = None,
        codon_cache: CodonCache | None = None,
        readthrough: Readthrough | None = None,
    ) -> None:
        """Initialize the annotator with a compiled rule set and its optional caches.

        The ``zone_map`` and ``codon_cache`` can be either coming from annotator context or supplied when
        creating the annotator.

        Args:
            ruleset (RulesetT): The compiled rule set this annotator will evaluate.
            annotator_context (AnnotatorContext | None): Run-scoped caches to evaluate against,
                or None (the default) for a fresh, empty one owned by this annotator alone.
            zone_map (ZoneMapCache | None): Zone-map cache to enable. Overrides whatever
                ``annotator_context`` carried; None leaves that untouched.
            codon_cache (CodonCache | None): Shared-codon cache to enable. Overrides whatever
                ``annotator_context`` carried; None leaves that untouched.
            readthrough (Readthrough | None): Stop-codon readthrough level for this run. Overrides
                whatever ``annotator_context`` carried; None leaves that untouched.
        """
        self.ruleset = ruleset
        self.annotator_context = annotator_context if annotator_context is not None else AnnotatorContext()
        if zone_map is not None:
            self.annotator_context.zone_map_cache = zone_map
        if codon_cache is not None:
            self.annotator_context.codon_cache = codon_cache
        if readthrough is not None:
            self.annotator_context.readthrough = readthrough
        self.last_trace: AnnotationTrace | None = None

    @property
    def cache_hits(self) -> int:
        """Return candidate-cache hits recorded on this annotator's context.

        Returns:
            int: The hit count. Shared with any other annotator given the same context.
        """
        return self.annotator_context.cache_hits

    @property
    def cache_misses(self) -> int:
        """Return candidate-cache misses recorded on this annotator's context.

        Returns:
            int: The miss count. Shared with any other annotator given the same context.
        """
        return self.annotator_context.cache_misses

    @abstractmethod
    def annotate_feature_overlaps(self, overlaps: Sequence[FeatureOverlap]) -> list[list[OverlapConsequence]]:
        """Calculate consequences for a batch of overlaps at once.

        Args:
            overlaps (Sequence[FeatureOverlap]): The overlaps to annotate, in order.

        Returns:
            list[list[OverlapConsequence]]: One consequence list per overlap, same order.
        """
        raise NotImplementedError

    def annotate_feature_overlap(self, overlap: FeatureOverlap) -> list[OverlapConsequence]:
        """Calculate consequences for a feature overlap.

        Supplied in terms of the batch method, so implementing ``annotate_feature_overlaps`` is
        enough to get both halves of the contract.

        **If your batch method is a loop over this one, you must override this one too.** A
        strategy that evaluates per overlap (``linear``/``staged``) implements the real work here
        and a plain loop there. Implementing *only* a looping batch method leaves this delegation
        pointing back at it, and the pair recurses until RecursionError.

        Args:
            overlap (FeatureOverlap): The overlap to annotate.

        Returns:
            list[OverlapConsequence]: The consequences that apply. May be empty.
        """
        return self.annotate_feature_overlaps([overlap])[0]
