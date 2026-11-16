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
"""Shared-codon LRU cache: amortises the expensive fact tier across variants that edit the
same codon window with the same edit.

Instead of module-level ``_GLOBAL_CACHE``, **instance-scoped**:
This repository's convention is a cache owned by whatever constructs it, not a process-global
singleton. Program-lifetime sharing is recovered the same way sharing across calls is here -
by constructing one CodonEditLRUCache at the pipeline level and passing it to every
annotator (or future cache-aware annotator) that should share it. The full-edit key below is what
keeps that sharing order-insensitive.

Read this before changing the key:
CodonEditKey.build and the fact computation it keys (``_codon_alleles`` fact) must never drift
apart, which is why the key is built from the same functions the fact reads (``Mapper.codon_window`` and
``variant_genomic_span``) rather than re-derived from scratch here.
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
    """The full-edit cache key for one codon window.

    A NamedTuple, not a value class: this is a **dict key** (see ``CodonEditLRUCache._cache``), so
    it has to be hashable and compare by value, and it is built once per coding overlap on the hot
    path. Staying a tuple keeps that cheap, a frozen attrs class measured ~5x the cost to build,
    hash and look up, while the field names here document what a bare tuple could not.

    Deliberately over-inclusive: Over-inclusion only ever costs hit rate; under-inclusion is the
    silent-wrong-answer failure mode.

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
        """Build the key for the codon window ``fact_context``'s variant touches.

        Args:
            fact_context (FactContext): The overlap being evaluated. Reads (and memoises, at the
                usual MAPPING cost) the ``_cds_footprint`` fact as a side effect.

        Returns:
            CodonEditKey | None: The cache key, or None if the variant touches no whole codon (no
                footprint, or the footprint maps to an empty/out-of-range codon window).
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
    """The cached MAPPING+SEQUENCE result for one codon-window edit.

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
    """Interface for codon/peptide-allele caches.

    ``get_codon_edit`` and ``clear`` are the storage contract an implementation must supply.
    ``prefill`` is not. Its body is generic, so it is provided here as a default, and an implementation
    that *inherits* this Protocol gets it for free. It stays ``@runtime_checkable``, so a purely
    structural implementation still satisfies ``isinstance``. Such an implementation just has to supply
    its own ``prefill``.
    """

    def get_codon_edit(self, key: CodonEditKey, loader: CodonEditLoader) -> CodonEdit:
        """Returns the cached codon edit for ``key``, computing it via ``loader`` on a miss.

        Args:
            key (CodonEditKey): The full-edit cache key (see CodonEditKey.build).
            loader (CodonEditLoader): Loader used when ``key`` is not cached.

        Returns:
            CodonEdit: The cached or newly loaded codon edit.
        """
        raise NotImplementedError

    def clear(self) -> None:
        """Clear cached state."""
        raise NotImplementedError

    def prefill(self, fact_context: FactContext[TranscriptVariationAllele]) -> None:
        """Gate, key, and inject ``_codon_alleles``/``_peptide_alleles`` from this cache.

        The body is generic: it reaches storage only through ``get_codon_edit``, which is why it is
        a default on the Protocol rather than an LRU-specific method.

        A no-op when the overlap is not coding, when the variant does not overlap the transcript's
        exons (``translation is not None`` + ``overlap_feature`` + ``overlap_exon`` - free after a
        zone-map prefill, harmlessly charged the normal way otherwise), or when the variant touches no
        whole codon.

        Must run before the cheap facts are computed and after any zone-map prefill; injecting once
        the sequence tier has already been reached is too late, the charge has already landed.
        AnnotatorContext.prepare is what guarantees that ordering.

        Args:
            fact_context (FactContext): The overlap being evaluated.
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
    """An in-memory LRU cache for codon/peptide alleles, keyed on the full edit.

    Instance-scoped: construct one and pass it to every annotator that should share it.

    ``eq=False`` keeps identity semantics: this is a shared, mutable instance, not a value object,
    so two caches with equal counters must not compare equal.

    Attributes:
        enabled (bool): Whether the cache is enabled. Default is True.
        max_entries (int): The maximum number of entries to keep in the cache. Default is 65536.
        hit_count (int): Number of get_codon_edit calls served from the cache so far.
        miss_count (int): Number of get_codon_edit calls that invoked ``loader`` so far.
    """

    enabled: bool = True
    max_entries: int = 65536
    hit_count: int = field(default=0, init=False)
    miss_count: int = field(default=0, init=False)
    _cache: OrderedDict[CodonEditKey, CodonEdit] = field(factory=OrderedDict, init=False, repr=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    def get_codon_edit(self, key: CodonEditKey, loader: CodonEditLoader) -> CodonEdit:
        """Returns the cached codon edit for ``key``, computing and storing it via ``loader`` on a miss.

        Args:
            key (CodonEditKey): The full-edit cache key (see CodonEditKey.build).
            loader (CodonEditLoader): Loader used when ``key`` is not cached. Not called on a hit.

        Returns:
            CodonEdit: The cached or newly loaded codon edit.
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
