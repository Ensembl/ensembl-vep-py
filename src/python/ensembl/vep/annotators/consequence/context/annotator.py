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
"""Run-level state shared across overlap annotations.

``FactContext`` belongs to one variant-feature overlap. ``AnnotatorContext`` spans the annotation
run and enables cross-overlap cache reuse.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from attrs import define, field

from ensembl.vep.common import DEFAULT_DISTANCE, Readthrough

if TYPE_CHECKING:
    from ensembl.vep.models.overlap import TranscriptVariationAllele
    from ensembl.vep.annotators.consequence.context.fact import FactContext
    from ensembl.vep.annotators.consequence.rules.rulesets import FactRule
    from ensembl.vep.caches.codon_cache import CodonCache
    from ensembl.vep.caches.mapper_cache import MapperIndexCache
    from ensembl.vep.caches.zone_map import ZoneMapCache

__all__ = ["AnnotatorContext"]


@define(kw_only=True, auto_attribs=True)
class AnnotatorContext:
    """Manage run-scoped caches shared across overlap annotations.

    Attributes:
        zone_map_cache (ZoneMapCache | None): Per-transcript positional-fact zone maps or None to
            disables zone-map prefill.
        codon_cache (CodonCache | None): Shared-codon edit cache or None to disable codon reuse; the
            relevant codons are then computed per overlap.
        readthrough (Readthrough): How far translation may run past a destroyed stop codon, for
            the whole run.
        upstream_distance (int): Bases before a feature's 5' end that still count as upstream.
        downstream_distance (int): Bases past a feature's 3' end that still count as downstream.
        mapper_index_cache (MapperIndexCache | None): Byte-budget LRU of per-transcript
            genomic<->CDS offset accelerators.
        digest_candidate_cache (dict[str, tuple[FactRule, ...]]): Cheap-guard digest string ->
            the rules that survive it. Used by the annotators (e.g., linear) that build a guard digest per
            overlap. Always present -- it costs one dict and is never a
            correctness risk, since the key is a total function of the guard values.
        bitmask_candidate_cache (dict[int, tuple[FactRule, ...]]): Cheap-guard bitmask -> the
            rules that survive it. The same idea keyed by an int rather than a string, for the
            annotators that already compute guard flags as bits.
        cache_hits (int): Candidate-cache hits, summed across both candidate caches.
        cache_misses (int): Candidate-cache misses, summed across both candidate caches.
    """

    zone_map_cache: ZoneMapCache | None = field(default=None)
    codon_cache: CodonCache | None = field(default=None)
    readthrough: Readthrough = field(default=Readthrough.NONE)
    upstream_distance: int = field(default=DEFAULT_DISTANCE)
    downstream_distance: int = field(default=DEFAULT_DISTANCE)
    mapper_index_cache: MapperIndexCache | None = field(default=None)
    digest_candidate_cache: dict[str, tuple[FactRule, ...]] = field(factory=dict, init=False, repr=False)
    bitmask_candidate_cache: dict[int, tuple[FactRule, ...]] = field(factory=dict, init=False, repr=False)
    cache_hits: int = field(default=0, init=False)
    cache_misses: int = field(default=0, init=False)

    def prepare(
        self, fact_context: FactContext[TranscriptVariationAllele], *, apply_zone_map: bool = True
    ) -> None:
        """Apply every enabled cache to one overlap, before its facts are computed.

        The single place the caches' ordering is written down. Annotators call this once per overlap and
        are otherwise cache-agnostic. Turning a cache on or off is a change to the context a caller will
        make, not to the annotator.

        Order is load-bearing and not arbitrary:

           assigning it invalidates the overlap's cached CDS alleles -- so it has to land before
        1. **Readthrough level first.** It changes what the sequence tier will compute, and
           anything reads them, for the same reason the codon cache has to land before the
           sequence tier are evaluated.
        2. **Zone map second.** It prefills the positional tier, including ``overlap_feature`` and
           ``overlap_exon``, which is exactly what the codon cache's gate reads next, so running
           it first makes that gate free rather than evaluated.
        3. **Codon cache third, and still before any rule is evaluated.** Injecting once the
           sequence tier has already been reached is too late: the evaluation cost has been charged
           and the cache saves nothing.

        Args:
            fact_context (FactContext): The overlap about to be evaluated.
            apply_zone_map (bool): False for a caller that will take the zone map's answers over a whole
                batch of overlaps instead (ZoneMapCache.columns), for example, a batch annotator (columnar).
                Doing both would compute every zone fact twice.
        """
        # Stamp the run's readthrough level onto the overlap before anything reads its sequence.
        # Guarded because assigning invalidates the overlap's cached CDS/peptide alleles, and the
        # default-level case must stay a true no-op.
        transcript_overlap = fact_context.overlap
        if transcript_overlap.readthrough is not self.readthrough:
            transcript_overlap.readthrough = self.readthrough

        # Both caches expose the same prefill(fact_context) shape, so applying them is symmetric;
        # each one's own gating and keying rules live with the cache rather than here.
        if apply_zone_map and self.zone_map_cache is not None:
            self.zone_map_cache.prefill(fact_context)
        if self.codon_cache is not None:
            self.codon_cache.prefill(fact_context)

    def reset_candidate_caches(self) -> None:
        """Drop both candidate caches and their counters.

        The candidate caches are keyed on cheap-guard values *for a particular rule set*, so a
        caller that swaps the rule set under a shared context must clear them or the second rule
        set will be served the first one's surviving rules.
        """
        self.digest_candidate_cache.clear()
        self.bitmask_candidate_cache.clear()
        self.cache_hits = 0
        self.cache_misses = 0
