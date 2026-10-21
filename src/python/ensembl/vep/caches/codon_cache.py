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
"""Cache codon and peptide edits shared across consequence annotations.

Caches are instance-scoped and keyed by the complete codon-window edit. Key construction uses the
same mapping functions as the ``_codon_alleles`` fact so cache identity cannot drift from fact
evaluation.
"""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Callable
from threading import RLock
from typing import TYPE_CHECKING, NamedTuple, Protocol, TypeAlias, cast, runtime_checkable

from attrs import define, field, frozen

from ensembl.vep.annotators.consequence.facts.base import get_fact
from ensembl.vep.annotators.consequence.facts.helper import variant_genomic_span
from ensembl.vep.trace.events import FactOrigin

if TYPE_CHECKING:
    from ensembl.vep.models.overlap import TranscriptVariationAllele
    from ensembl.vep.annotators.consequence.context.fact import FactContext

__all__ = [
    "CodonCache",
    "CodonEdit",
    "CodonEditKey",
    "CodonEditLRUCache",
    "CodonEditLoader",
]


class CodonEditKey(NamedTuple):
    """Identify one codon-window edit.

    The key is deliberately over-inclusive: extra fields can reduce the hit rate, while omitted
    fields can reuse an incorrect consequence. A ``NamedTuple`` keeps hot-path construction and
    hashing inexpensive while retaining named fields.

    Attributes:
        seq_region_name (str): Region the codon window sits on.
        codon_window_genomic_positions (tuple[int, ...]): Genomic positions of the whole codon(s).
        padding (int): Window bases with no genomic position -- the ``N`` before a CDS that starts
            part-way through a codon. Without it that first codon's key could equal a partial last
            codon's on another transcript over the same bases.
        strand (int): Reading direction, as ``Strand.value``.
        variant_min (int): Low genomic bound of the edit, fixing where in the window it falls.
        variant_max (int): High genomic bound of the edit.
        reference_allele (str): Reference bases.
        alternate_allele (str): Alternate bases.
    """

    seq_region_name: str
    codon_window_genomic_positions: tuple[int, ...]
    padding: int
    strand: int
    variant_min: int
    variant_max: int
    reference_allele: str
    alternate_allele: str

    @classmethod
    def build(cls, fact_context: FactContext[TranscriptVariationAllele]) -> CodonEditKey | None:
        """Build a key for the codon window touched by a variant.

        Args:
            fact_context (FactContext): Overlap being evaluated. Reads and memoises the
                ``_cds_footprint`` fact.

        Returns:
            CodonEditKey | None: Cache key, or ``None`` if the variant touches no codon window.
        """
        transcript = fact_context.overlap.feature
        footprint = cast("tuple[int, int] | None", get_fact(fact_context, "_cds_footprint"))
        mapper = transcript.cds_mapper
        if mapper is None:
            # No CDS mapper means no codon coordinates to key on. The caller's gate makes this
            # unreachable for a coding transcript, but the type allows it and an AttributeError
            # here would surface as a cache failure rather than a missing key.
            return None
        window = mapper.codon_window(footprint)
        if window is None:
            return None

        accelerator = None
        annotator_context = fact_context.annotator_context
        if (
            annotator_context is not None
            and annotator_context.mapper_index_cache is not None
            and transcript.cds_mapper is not None
        ):
            accelerator = annotator_context.mapper_index_cache.get_index(
                transcript, transcript.cds_mapper.build_index
            )

        window_start, window_end = window
        codon_window_genomic_positions = tuple(
            genomic_position
            for genomic_position in (
                transcript.cds_offset_to_genomic(offset, accelerator=accelerator)
                for offset in range(window_start, window_end)
            )
            if genomic_position is not None
        )

        variant = fact_context.overlap.overlapping_feature
        variant_min, variant_max = variant_genomic_span(fact_context)
        seq_region_name = transcript.slice.seq_region_name if transcript.slice is not None else ""
        reference_allele = variant.reference_allele or ""
        alternate_allele = variant.alternate_allele or ""
        return cls(
            seq_region_name,
            codon_window_genomic_positions,
            (window_end - window_start) - len(codon_window_genomic_positions),
            transcript.strand.value,
            variant_min,
            variant_max,
            reference_allele,
            alternate_allele,
        )


CodonEditLoader: TypeAlias = Callable[[], "CodonEdit"]


@frozen(kw_only=True, auto_attribs=True)
class CodonEdit:
    """Store the mapping and sequence result for one codon-window edit.

    Attributes:
        reference_codon (str): The reference codon window, transcript-sense.
        alternate_codon (str): ``reference_codon`` with the variant allele spliced in.
        reference_peptide (str): The amino acid translation of ``reference_codon``.
        alternate_peptide (str): The amino acid translation of ``alternate_codon``.
    """

    reference_codon: str
    alternate_codon: str
    reference_peptide: str
    alternate_peptide: str


@runtime_checkable
class CodonCache(Protocol):
    """Define storage used to reuse codon and peptide edits.

    Implementations provide ``get_codon_edit`` and ``clear``. The default ``prefill`` method
    contains the storage-independent fact integration.
    """

    def get_codon_edit(self, key: CodonEditKey, loader: CodonEditLoader) -> CodonEdit:
        """Return the cached codon edit, loading it on a miss.

        Args:
            key (CodonEditKey): Full-edit cache key.
            loader (CodonEditLoader): Loader used when ``key`` is not cached.

        Returns:
            CodonEdit: Cached or newly loaded codon edit.
        """
        raise NotImplementedError

    def clear(self) -> None:
        """Clear cached state."""
        raise NotImplementedError

    def prefill(self, fact_context: FactContext[TranscriptVariationAllele]) -> None:
        """Prefill codon and peptide facts when the overlap is coding.

        Run after zone-map prefill and before fact evaluation so cache hits avoid mapping and
        sequence work. Non-coding and non-exonic overlaps are ignored.

        Args:
            fact_context (FactContext): Overlap to prefill.
        """
        if fact_context.overlap.feature.translation is None:
            return
        if not get_fact(fact_context, "overlap_feature"):
            return
        if not get_fact(fact_context, "overlap_exon"):
            return

        key = CodonEditKey.build(fact_context)
        if key is None:
            return

        def _load_codon_edit() -> CodonEdit:
            reference_codon, alternate_codon = cast(
                "tuple[str, str]", get_fact(fact_context, "_codon_alleles")
            )
            reference_peptide, alternate_peptide = cast(
                "tuple[str, str]", get_fact(fact_context, "_peptide_alleles")
            )
            return CodonEdit(
                reference_codon=reference_codon,
                alternate_codon=alternate_codon,
                reference_peptide=reference_peptide,
                alternate_peptide=alternate_peptide,
            )

        codon_edit = self.get_codon_edit(key, _load_codon_edit)
        # On a hit, _load_codon_edit above never ran, so _codon_alleles/_peptide_alleles were never
        # computed (or charged) for this context -- this prefill is the entire mechanism, and the
        # only thing that makes either name known, correctly attributed to FactOrigin.CACHED.
        # On a miss, _load_codon_edit already populated the memo via get_fact; prefilling the same
        # values again here is redundant but harmless, and passing origin=CACHED never mislabels it,
        # because prefill only attributes the names it actually introduced to the memo.
        fact_context.prefill(
            {
                "_codon_alleles": (codon_edit.reference_codon, codon_edit.alternate_codon),
                "_peptide_alleles": (codon_edit.reference_peptide, codon_edit.alternate_peptide),
            },
            origin=FactOrigin.CACHED,
        )


@define(kw_only=True, auto_attribs=True, eq=False)
class CodonEditLRUCache(CodonCache):
    """Store codon edits in an instance-scoped LRU cache.

    ``eq=False`` preserves identity semantics for mutable shared state.

    Attributes:
        enabled (bool): Whether caching is enabled.
        max_entries (int): Maximum number of entries to retain.
        hit_count (int): Number of lookups served from the cache.
        miss_count (int): Number of lookups that invoked ``loader``.
    """

    enabled: bool = True
    max_entries: int = 65536
    hit_count: int = field(default=0, init=False)
    miss_count: int = field(default=0, init=False)
    _cache: OrderedDict[CodonEditKey, CodonEdit] = field(factory=OrderedDict, init=False, repr=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    def get_codon_edit(self, key: CodonEditKey, loader: CodonEditLoader) -> CodonEdit:
        """Return the cached codon edit, loading and storing it on a miss.

        Args:
            key (CodonEditKey): Full-edit cache key.
            loader (CodonEditLoader): Loader used when ``key`` is not cached.

        Returns:
            CodonEdit: Cached or newly loaded codon edit.
        """
        if not self.enabled:
            self.miss_count += 1
            return loader()

        with self._lock:
            cached = self._cache.get(key)
            if cached is not None:
                self._cache.move_to_end(key)
                self.hit_count += 1
                return cached

        codon_edit = loader()

        with self._lock:
            self._cache[key] = codon_edit
            self._cache.move_to_end(key)
            self.miss_count += 1
            while len(self._cache) > self.max_entries:
                self._cache.popitem(last=False)

        return codon_edit

    def clear(self) -> None:
        """Clear the cache contents."""
        with self._lock:
            self._cache.clear()
