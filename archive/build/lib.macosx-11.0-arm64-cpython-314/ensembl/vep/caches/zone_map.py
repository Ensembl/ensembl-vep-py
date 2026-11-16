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
"""Zone maps: compile a transcript's zone.

Most of what facts computes is positional (``overlap_exon``, ``overlap_intron``,
``overlap_intron_interior``, the six splice-offset facts, the UTR tests, the start/stop codon overlaps,
``overlap_cds``, ``is_upstream``/``is_downstream``.
For a fixed transcript each of those is a fixed function of a genomic coordinate and a transcript is reused
across every variant that touches it. So the interval scans re-derive, once per variant, something that
only ever depended on the transcript.

The fix is to precompute it. Every positional fact is a union of genomic intervals, so the endpoints of
all of them cut the transcript's span (plus its upstream/downstream flanks) into a small number of
zones inside which every positional fact is constant. Each zone gets a bitmask; the zones are
stored as a sorted array of start coordinates (ZoneMap). Answering the positional tier for a
variant is then a ``bisect`` plus a short OR, instead of a dozen interval scans.

The answers are injected via FactContext.prefill, which seeds the memo table without charging the
seeded facts' cost.

This lives under ``caches/`` because that is what it is. A zone map is built once per
transcript and reused across every variant that touches it, keyed by transcript, with hit/miss
counters. ``ZoneMapCache`` is a byte-budgeted LRU (``DEFAULT_ZONE_MAP_BYTES``): ``starts`` is a
packed ``array('q')`` and ``masks`` a tuple of shared boxed ints (see ``_MASK_INTERN_CACHE``), and
the least-recently-used map is evicted first once the budget is exceeded. The zone map itself
stays opt-in (``--zone_map``); it pays for itself once a transcript is touched by more than one
variant.

A variant covers ``min_position..max_position`` (spans, and not points), which may straddle several
zones. Every positional fact is a "does the variant *overlap* X" test, so the variant's value is
the union of the covered zones' masks, hence the OR in ZoneMap.lookup. Facts that are not pure unions
are handled explicitly rather than approximately, for example, ``is_upstream``/``is_downstream`` return
False whenever the variant touches the transcript at all, so the flank bit is ANDed with
``not overlap_feature`` in ZoneMapCache.prefill.

Limits: this does not touch the expensive tier. ``_cds_footprint``, ``_codon_alleles``,
``_peptide_alleles`` and everything downstream depend on the alt allele, so nothing about them
can be folded into the transcript. It is also a loss for a transcript touched by exactly one
variant - the build is a sort over a few dozen interval endpoints plus a mask per zone, cheaper
than one full positional pass but not free, and it only amortises from the second variant onward.
"""

from __future__ import annotations

import sys
import time
from array import array
from typing import ClassVar
from bisect import bisect_right
from collections import OrderedDict, defaultdict
from collections.abc import Sequence
from threading import RLock

from attrs import define, field, frozen

from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.facts.helper import (
    intron_offsets_to_genomic_span,
    variant_genomic_span,
)
from ensembl.vep.common import DEFAULT_DISTANCE, Strand
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele

__all__ = [
    "DEFAULT_ZONE_MAP_BYTES",
    "MEMO_FACT_NAMES",
    "SERVED_FACT_NAMES",
    "ZoneMap",
    "ZoneMapCache",
]

DEFAULT_ZONE_MAP_BYTES = 256 * 1024 * 1024
"""Default byte budget for ``ZoneMapCache``, shared by the config default and every call site."""

POSITIONAL_FACT_NAMES: tuple[str, ...] = (
    "overlap_feature",
    "overlap_exon",
    "overlap_intron",
    "overlap_intron_interior",
    "overlap_intron_boundary",
    "overlap_donor_site",
    "overlap_acceptor_site",
    "overlap_donor_5th_base",
    "overlap_donor_region",
    "overlap_polypyrimidine_tract",
    "overlap_splice_region",
    "is_upstream",
    "is_downstream",
)
"""Every fact whose value a variant's genomic span determines against a fixed transcript, purely
from exon/intron/flank geometry."""


_ZONE_FACT_NAMES: tuple[str, ...] = POSITIONAL_FACT_NAMES[:-2] + (
    "overlap_five_prime_utr",
    "overlap_three_prime_utr",
    "overlap_start_codon",
    "overlap_stop_codon",
    "overlap_cds",
)
"""Every fact that is served straight out of the zone mask, as a plain union over the zones the
variant covers. Bit order is this tuple's order. The five facts added here come from the
coding-region tier, but reduce to a plain interval union just like the positional facts above
them."""

# Flank bits: a union over zones, but only meaningful when the variant misses the transcript
# entirely, so they are not served directly out of the zone mask (see prefill()).
_UPSTREAM_FACT_NAME = "is_upstream"
_DOWNSTREAM_FACT_NAME = "is_downstream"

_MASK_FACT_NAMES: tuple[str, ...] = _ZONE_FACT_NAMES + (_UPSTREAM_FACT_NAME, _DOWNSTREAM_FACT_NAME)
_FACT_BIT: dict[str, int] = {
    fact_name: 1 << bit_index for bit_index, fact_name in enumerate(_MASK_FACT_NAMES)
}

# Everything this annotator prefills, i.e. every fact computation it removes from the bill.
SERVED_FACT_NAMES: tuple[str, ...] = _MASK_FACT_NAMES

MEMO_FACT_NAMES: tuple[str, ...] = (
    "overlap_feature",
    "overlap_exon",
    "overlap_cds",
    "overlap_start_codon",
    "overlap_stop_codon",
)
"""Served facts that must also reach each overlap's memo, not just a column.

Exactly the served facts that appear in some other fact's ``requires``. The guard the ``fact``
decorator installs (consequence/facts/base.py) reads them through ``get_fact``, so a column alone
would leave the guard recomputing them. The rest are read only by rule formulas, which a columnar
annotator answers from the column directly."""

_MASK_INTERN_CACHE: dict[int, int] = {}
"""Process-lifetime cache of boxed zone-mask ints, shared by every ``ZoneMap.build`` call.

Real input has very few distinct masks -- a median of 2, max 8, per neighbouring-variant chunk
(see ``ZoneMapCache._implied_facts_by_mask``) -- so pointing every zone map's ``masks`` tuple at
the same handful of shared int objects keeps memory close to array-packed storage while a tuple
lookup pays no per-access unboxing cost, unlike ``array.array``. A microbenchmark over a 500-2000
zone map backs this: array/array lookups took 34-41% longer than array-starts/deduped-masks-tuple,
for +4 B/zone (~16 B/zone total, well under a plain tuple's ~48 B/zone). Bounded by the number of
distinct masks the ``_FACT_BIT`` layout can produce (2 ** len(_MASK_FACT_NAMES) at worst), so it
never needs clearing.
"""


_ESTIMATED_ENTRY_OVERHEAD_BYTES = 250
"""Fixed allowance ``ZoneMap.estimated_bytes`` adds on top of its containers' own sizes.

Covers what a per-container ``sys.getsizeof`` cannot see: the ``ZoneMap`` object itself, the
cache's ``(stable_id, version)`` key tuple, and its ``OrderedDict``/``_sizes`` entries. Calibrated
against ``tracemalloc``-measured insertions into a real ``ZoneMapCache`` (see
``tests/caches/test_zone_map.py``); the gap is close to size-invariant, since none of those four
things grow with zone count.
"""


def _int64_array(values: Sequence[int]) -> array[int]:
    """Pack zone start positions into a signed 64-bit array."""
    return array("q", values)


def _dedup_masks(values: Sequence[int]) -> tuple[int, ...]:
    """Return ``values`` as a tuple, each slot pointing at one shared boxed int.

    See ``_MASK_INTERN_CACHE``.
    """
    return tuple(_MASK_INTERN_CACHE.setdefault(value, value) for value in values)


@frozen(kw_only=True, auto_attribs=True, hash=False)
class ZoneMap:
    """A sorted partition of one transcript's span (plus flanks) into positional-fact zones.

    ``starts[i]`` is the first genomic position of zone ``i`` and ``masks[i]`` its fact bitmask;
    zone ``i`` runs to ``starts[i + 1] - 1`` (the last zone runs to highest_position).
    Positions outside ``lowest_position..highest_position`` carry no positional facts at all.

    Not hashable (``hash=False``): ``starts`` is an ``array.array``, which is not hashable either,
    so leaving attrs' default frozen-class ``__hash__`` in place would advertise a hash that raises
    the moment it is called.

    Attributes:
        lowest_position (int): Lowest genomic position covered by any zone.
        highest_position (int): Highest genomic position covered by any zone.
        starts (array[int]): Ascending zone start positions, packed as a signed 64-bit array (~8
            bytes a zone instead of a boxed-int tuple's ~40+).
        masks (tuple[int, ...]): Fact bitmask for the zone starting at the same index in
            ``starts``. A tuple, not a packed array: masks repeat heavily across zones and maps, so
            ``_dedup_masks`` points every slot at one shared boxed int, which keeps memory close to
            array-packed storage while avoiding the per-lookup unboxing cost an array pays.
        inverted (tuple[tuple[int, int, int], ...]): ``(bit, start, end)`` triples with
            ``start > end``, a degenerate interval a union cannot represent, kept out of the
            zone array and matched instead by a variant that spans ``[end, start]``. None of
            ZoneMap.build's own intervals produce one today; kept so a future guard relaxation
            cannot silently produce a wrong zone map rather than a wrong answer here.
    """

    lowest_position: int
    highest_position: int
    starts: array[int] = field(converter=_int64_array, default=())
    masks: tuple[int, ...] = field(converter=_dedup_masks, default=())
    inverted: tuple[tuple[int, int, int], ...] = field(default=())

    def lookup(self, min_position: int, max_position: int) -> int:
        """Returns the union of the masks of every zone the span overlaps.

        Args:
            min_position (int): Lower bound (inclusive) of the genomic span to look up.
            max_position (int): Upper bound (inclusive) of the genomic span to look up.

        Returns:
            int: The OR of every overlapping zone's mask, plus any inverted interval the span
                fully covers.
        """
        mask = 0
        if max_position >= self.lowest_position and min_position <= self.highest_position:
            clipped_min = max(min_position, self.lowest_position)
            clipped_max = min(max_position, self.highest_position)
            start_index = bisect_right(self.starts, clipped_min) - 1
            end_index = bisect_right(self.starts, clipped_max) - 1
            for zone_index in range(start_index, end_index + 1):
                mask |= self.masks[zone_index]
        for bit, interval_start, interval_end in self.inverted:
            if min_position <= interval_end and max_position >= interval_start:
                mask |= bit
        return mask

    def estimated_bytes(self) -> int:
        """Estimate this zone map's heap footprint, for LRU byte accounting.

        ``sys.getsizeof`` on ``starts``/``masks``/``inverted`` is O(1) -- it reports each
        container's own buffer size without walking its elements -- and exact for ``starts``
        (a packed array) and close to exact for ``masks`` (a tuple of pointers to shared ints, so
        its own boxed-int cost is not charged here; see ``_MASK_INTERN_CACHE``). What it cannot see
        is added as one fixed allowance: the ``ZoneMap`` object itself, the cache's
        ``(stable_id, version)`` key tuple, and its ``OrderedDict``/``_sizes`` entries.

        Returns:
            int: The estimated size in bytes.
        """
        return (
            sys.getsizeof(self.starts)
            + sys.getsizeof(self.masks)
            + sys.getsizeof(self.inverted)
            + _ESTIMATED_ENTRY_OVERHEAD_BYTES
        )

    @classmethod
    def build(
        cls,
        transcript: Transcript,
        upstream_distance: int = DEFAULT_DISTANCE,
        downstream_distance: int = DEFAULT_DISTANCE,
    ) -> ZoneMap:
        """Compile a transcript into its positional-fact zone partition.

        Decomposes every POSITIONAL_FACT_NAMES fact directly from the transcript's exon/intron/flank
        geometry. Uses intron_offsets_to_genomic_span for the splice-offset arithmetic, the same function
        intron_offset_hit uses for the per-variant tests, so the two cannot diverge. Adds the coding-region
        intervals directly from the transcript's own CDS/codon mapping.

        Args:
            transcript (Transcript): The transcript to compile.
            upstream_distance (int): Bases before the 5' end that ``is_upstream`` covers.
            downstream_distance (int): Bases past the 3' end that ``is_downstream`` covers.

        Returns:
            ZoneMap: The compiled zone partition.
        """
        spans: list[tuple[int, int, int]] = []  # (bit, start, end)
        inverted: list[tuple[int, int, int]] = []

        def add(fact_name: str, start: int, end: int) -> None:
            bit = _FACT_BIT[fact_name]
            if start <= end:
                spans.append((bit, start, end))
            else:
                inverted.append((bit, start, end))

        add("overlap_feature", transcript.start, transcript.end)

        for exon in transcript.spliced_exons:
            add("overlap_exon", exon.start, exon.end)

        for intron in transcript.introns:
            add("overlap_intron", intron.start, intron.end)
            # Inverted for an intron shorter than five bases, which add() keeps out of the zones:
            # the fact's own genomic_overlap test then matches only a span covering the whole window.
            add("overlap_intron_interior", intron.start + 2, intron.end - 2)
            # Genomic, strand-agnostic windows -- no strand resolution needed, unlike every
            # offset-based span below. See overlap_intron_boundary's own docstring (position.py) for
            # the known coverage gap relative to overlap_polypyrimidine_tract/overlap_splice_region.
            add("overlap_intron_boundary", intron.start - 4, intron.start + 8)
            add("overlap_intron_boundary", intron.end - 8, intron.end + 4)
            strand = transcript.strand
            add("overlap_donor_site", *intron_offsets_to_genomic_span(intron, strand, 1, 2, True))
            add("overlap_acceptor_site", *intron_offsets_to_genomic_span(intron, strand, 1, 2, False))
            add("overlap_donor_5th_base", *intron_offsets_to_genomic_span(intron, strand, 5, 5, True))
            add("overlap_donor_region", *intron_offsets_to_genomic_span(intron, strand, 3, 6, True))
            add(
                "overlap_polypyrimidine_tract",
                *intron_offsets_to_genomic_span(intron, strand, 3, 17, False),
            )
            add("overlap_splice_region", *intron_offsets_to_genomic_span(intron, strand, 3, 8, True))
            add("overlap_splice_region", *intron_offsets_to_genomic_span(intron, strand, 3, 8, False))
            # The exonic side of each boundary is strand-independent, unlike the offset spans above.
            add("overlap_splice_region", intron.start - 3, intron.start - 1)
            add("overlap_splice_region", intron.end + 1, intron.end + 3)

        for fact_name, flank_start, flank_end in _flank_windows(
            transcript, upstream_distance, downstream_distance
        ):
            add(fact_name, flank_start, flank_end)

        # overlap_cds/overlap_five_prime_utr/overlap_three_prime_utr are not part of the geometry
        # decomposition above -- they reduce to a plain interval union over an already exon-clipped
        # genomic spans list (Transcript.cds_genomic_spans/five_prime_utr_genomic_spans/
        # three_prime_utr_genomic_spans), added directly from that list instead.
        for block_start, block_end in transcript.cds_genomic_spans or ():
            add("overlap_cds", block_start, block_end)
        for block_start, block_end in transcript.five_prime_utr_genomic_spans or ():
            add("overlap_five_prime_utr", block_start, block_end)
        for block_start, block_end in transcript.three_prime_utr_genomic_spans or ():
            add("overlap_three_prime_utr", block_start, block_end)

        for span_start, span_end in transcript.start_codon_genomic_spans() or ():
            add("overlap_start_codon", span_start, span_end)
        for span_start, span_end in transcript.stop_codon_genomic_spans() or ():
            add("overlap_stop_codon", span_start, span_end)

        if not spans:
            raise ValueError(
                f"no positional spans to build a zone map from for transcript {transcript.stable_id!r}"
            )

        lowest_position = min(start for _bit, start, _end in spans)
        highest_position = max(end for _bit, _start, end in spans)

        cut_points = {lowest_position, highest_position + 1}
        for _bit, start, end in spans:
            cut_points.add(start)
            cut_points.add(end + 1)
        zone_boundaries = sorted(
            point for point in cut_points if lowest_position <= point <= highest_position + 1
        )

        # Sweep the cut points once, counting open spans per bit, instead of testing every span
        # against every zone: a many-exon transcript has thousands of both. Every span start and
        # end + 1 is itself a cut point, so each event lands on a zone boundary.
        opening_bits: dict[int, list[int]] = defaultdict(list)
        closing_bits: dict[int, list[int]] = defaultdict(list)
        for bit, start, end in spans:
            opening_bits[start].append(bit)
            closing_bits[end + 1].append(bit)

        open_counts: dict[int, int] = defaultdict(int)
        mask = 0
        zone_starts: list[int] = []
        zone_masks: list[int] = []
        for zone_start in zone_boundaries[:-1]:
            for bit in closing_bits.get(zone_start, ()):
                open_counts[bit] -= 1
                if open_counts[bit] == 0:
                    mask &= ~bit
            for bit in opening_bits.get(zone_start, ()):
                open_counts[bit] += 1
                mask |= bit
            zone_starts.append(zone_start)
            zone_masks.append(mask)

        return cls(
            lowest_position=lowest_position,
            highest_position=highest_position,
            starts=zone_starts,
            masks=zone_masks,
            inverted=tuple(inverted),
        )


_FLANK_ONLY_MARGIN = 16
"""Distance outside a transcript beyond which only flank facts can hold.

Every non-flank span a zone map records lies within the transcript's own start and end; the margin
is slack, not a measured bound.
"""


def _flank_windows(
    transcript: Transcript, upstream_distance: int, downstream_distance: int
) -> tuple[tuple[str, int, int], ...]:
    """Return the ``is_upstream`` and ``is_downstream`` genomic windows of a transcript.

    The two distances need not be equal, so each window is measured from the end it belongs to:
    on the reverse strand the 5' end is the transcript's highest coordinate.
    """

    before, after = upstream_distance, downstream_distance
    if transcript.strand == Strand.REVERSE:
        before, after = after, before
    lower_flank = (transcript.start - before, transcript.start - 1)
    upper_flank = (transcript.end + 1, transcript.end + after)
    if transcript.strand == Strand.FORWARD:
        return (("is_upstream", *lower_flank), ("is_downstream", *upper_flank))
    return (("is_upstream", *upper_flank), ("is_downstream", *lower_flank))


def _flank_mask(
    transcript: Transcript,
    min_position: int,
    max_position: int,
    upstream_distance: int,
    downstream_distance: int,
) -> int:
    """Return the zone mask of a span that does not reach the transcript itself."""

    mask = 0
    for fact_name, flank_start, flank_end in _flank_windows(
        transcript, upstream_distance, downstream_distance
    ):
        if min_position <= flank_end and max_position >= flank_start:
            mask |= _FACT_BIT[fact_name]
    return mask


@define(kw_only=True, auto_attribs=True, eq=False)
class ZoneMapCache:
    """Per-transcript zone maps, byte-budgeted and evicted least-recently-used.

    The unit of reuse is the transcript, so one instance shared across a whole run is what makes
    the map pay for itself. Building it costs more than one positional pass and only amortises
    from the second variant touching that transcript onward. Construct one and hand it to every
    strategy that should share it.

    ``max_bytes`` follows ``MapperIndexLRUCache``'s shape: build outside the lock, insert and evict
    under it. ``eq=False`` keeps identity semantics (two caches with the same counters are not the
    same cache); attrs' default field-wise ``__eq__`` would otherwise ignore that this class is a
    shared, mutable instance, not a value object.

    Attributes:
        max_bytes (int): Byte budget for cached zone maps. Default is ``DEFAULT_ZONE_MAP_BYTES``.
        upstream_distance (int): The run's upstream pairing distance, which the flank zones cover.
        downstream_distance (int): The run's downstream pairing distance.
        hit_count (int): Number of ``zone_map`` calls served from the cache so far.
        miss_count (int): Number of ``zone_map`` calls that built a map so far.
        maps_built (int): Alias for ``miss_count``, kept for existing callers/counters.
        eviction_count (int): Number of entries evicted for being over budget so far.
        build_seconds (float): Total wall-clock time spent building zone maps so far.
        facts_served (int): Total number of fact values prefilled (uncharged) so far.
    """

    max_bytes: int = DEFAULT_ZONE_MAP_BYTES
    upstream_distance: int = DEFAULT_DISTANCE
    downstream_distance: int = DEFAULT_DISTANCE
    hit_count: int = field(default=0, init=False)
    miss_count: int = field(default=0, init=False)
    eviction_count: int = field(default=0, init=False)
    build_seconds: float = field(default=0.0, init=False)
    facts_served: int = field(default=0, init=False)
    _zone_maps: OrderedDict[tuple[str, int | None], ZoneMap] = field(
        factory=OrderedDict, init=False, repr=False
    )
    _sizes: dict[tuple[str, int | None], int] = field(factory=dict, init=False, repr=False)
    _total_bytes: int = field(default=0, init=False, repr=False)
    _lock: RLock = field(factory=RLock, init=False, repr=False)

    _implied_facts_by_mask: ClassVar[dict[int, tuple[str, ...]]] = {}
    """Memo: zone mask -> the facts that mask implies, flank suppression already applied.

    ``ZoneMap.lookup`` returns one packed mask per overlap, and the same mask recurs constantly --
    a 256-overlap chunk of the heavy-locus corpus has a *median of 2* distinct masks (max 8 over
    all 198 chunks), because neighbouring variants on one transcript land in the same zone.
    Resolving a mask once per distinct value instead of once per overlap takes the transpose in
    ``columns`` from ~156us to ~43us per chunk.

    Shared across every instance, deliberately: the mapping is a pure function of the module's
    ``_FACT_BIT`` layout, not of any one cache's state, so a per-instance table would just be
    rebuilt identically. Bounded by 2 ** len(_MASK_FACT_NAMES) entries in the worst case, and it
    degrades gracefully. On input with far more distinct masks than this corpus it converges on
    the old per-overlap cost plus one dict lookup, never worse by more than that.
    """

    @staticmethod
    def _implied_fact_names(mask: int) -> tuple[str, ...]:
        """Returns the fact names a zone mask implies, memoised on the mask.

        Args:
            mask (int): A packed zone mask, as returned by ``ZoneMap.lookup``.

        Returns:
            tuple[str, ...]: Every fact the mask sets, with the flank facts already suppressed if
                the variant overlaps the feature at all.
        """
        cache = ZoneMapCache._implied_facts_by_mask
        implied_fact_names = cache.get(mask)
        if implied_fact_names is None:
            implied = [fact_name for fact_name in _ZONE_FACT_NAMES if mask & _FACT_BIT[fact_name]]
            # Flanks are suppressed by any overlap with the feature, so they are not a plain
            # union over the covered zones.
            if not mask & _FACT_BIT["overlap_feature"]:
                implied.extend(
                    fact_name
                    for fact_name in (_UPSTREAM_FACT_NAME, _DOWNSTREAM_FACT_NAME)
                    if mask & _FACT_BIT[fact_name]
                )
            implied_fact_names = cache[mask] = tuple(implied)
        return implied_fact_names

    @property
    def maps_built(self) -> int:
        """Alias for ``miss_count``, kept for callers that count zone maps built."""
        return self.miss_count

    def zone_map(self, transcript: Transcript) -> ZoneMap:
        """Returns the (cached) zone map for a transcript, building it on first use.

        Keyed by ``(stable_id, version)`` rather than object identity, so two model objects
        describing the same transcript share one map. The build happens outside the lock, so one
        slow build cannot block other threads' hits; only the dict/counter mutation is serialised.
        A concurrent miss on the same key builds twice (both count), but the size accounting is
        exact: the losing build's bytes are subtracted before the winner's are added.

        Args:
            transcript (Transcript): The transcript to look up.

        Returns:
            ZoneMap: The transcript's zone partition.
        """
        cache_key = (transcript.stable_id, transcript.version)
        with self._lock:
            cached = self._zone_maps.get(cache_key)
            if cached is not None:
                self._zone_maps.move_to_end(cache_key)
                self.hit_count += 1
                return cached

        build_started = time.perf_counter()
        zone_map = ZoneMap.build(transcript, self.upstream_distance, self.downstream_distance)
        build_elapsed = time.perf_counter() - build_started
        size = zone_map.estimated_bytes()

        with self._lock:
            if cache_key in self._sizes:
                self._total_bytes -= self._sizes[cache_key]
            self._zone_maps[cache_key] = zone_map
            self._zone_maps.move_to_end(cache_key)
            self._sizes[cache_key] = size
            self._total_bytes += size
            self.miss_count += 1
            self.build_seconds += build_elapsed
            while self._total_bytes > self.max_bytes and len(self._zone_maps) > 1:
                evicted_key, _ = self._zone_maps.popitem(last=False)
                self._total_bytes -= self._sizes.pop(evicted_key, 0)
                self.eviction_count += 1
            return zone_map

    def clear(self) -> None:
        """Clear the cache contents."""
        with self._lock:
            self._zone_maps.clear()
            self._sizes.clear()
            self._total_bytes = 0

    def columns(self, fact_contexts: Sequence[FactContext]) -> dict[str, int]:
        """Returns every zone-map-served fact as a bit column over a whole chunk.

        The column-shaped counterpart to prefill, and the one that is worth using. prefill has to
        build a 16-entry dict per overlap and write it into that overlap's memo, after which the
        columnar annotators still call ``get_fact`` once per (fact, overlap) to read it back. But the
        answer is already in the shape a bit-parallel annotator wants.

        The five facts in MEMO_FACT_NAMES are additionally written to each overlap's memo, because
        they appear in some other fact's ``requires`` and the decorator's guard reads them through
        ``get_fact``. The other eleven never need to exist per overlap at all.

        Args:
            fact_contexts (Sequence[FactContext]): The chunk's contexts, batch-index order.

        Returns:
            dict[str, int]: Fact name -> bit column, bit ``i`` being the value for context ``i``.
        """
        columns: dict[str, int] = {fact_name: 0 for fact_name in SERVED_FACT_NAMES}

        for context_index, fact_context in enumerate(fact_contexts):
            variant_min, variant_max = fact_context.variant_span
            transcript = fact_context.overlap.feature
            if (
                variant_max < transcript.start - _FLANK_ONLY_MARGIN
                or variant_min > transcript.end + _FLANK_ONLY_MARGIN
            ):
                # Most overlaps on a real genome reach the transcript only through its flank.
                # Only the flank facts can hold there, so answer them without building a map.
                mask = _flank_mask(
                    transcript,
                    variant_min,
                    variant_max,
                    self.upstream_distance,
                    self.downstream_distance,
                )
            else:
                mask = self.zone_map(transcript).lookup(variant_min, variant_max)
            bit = 1 << context_index

            for fact_name in self._implied_fact_names(mask):
                columns[fact_name] |= bit

        # Only the facts some other fact's `requires` guard reads need to exist per overlap.
        for context_index, fact_context in enumerate(fact_contexts):
            bit = 1 << context_index
            fact_context.prefill({fact_name: bool(columns[fact_name] & bit) for fact_name in MEMO_FACT_NAMES})

        self.facts_served += len(fact_contexts) * len(columns)
        return columns

    def prefill(self, fact_context: FactContext[TranscriptVariationAllele]) -> None:
        """Write every zone-map-served fact into the context memo, uncharged.

        Args:
            fact_context (FactContext): The overlap to prefill facts for.
        """
        variant_min, variant_max = variant_genomic_span(fact_context)
        zone_map = self.zone_map(fact_context.overlap.feature)
        mask = zone_map.lookup(variant_min, variant_max)

        served_values: dict[str, object] = {
            fact_name: bool(mask & _FACT_BIT[fact_name]) for fact_name in _ZONE_FACT_NAMES
        }

        # The flanking facts are suppressed by any overlap with the transcript, so they are not a
        # plain union over the covered zones.
        variant_misses_transcript = not served_values["overlap_feature"]
        served_values[_UPSTREAM_FACT_NAME] = variant_misses_transcript and bool(
            mask & _FACT_BIT[_UPSTREAM_FACT_NAME]
        )
        served_values[_DOWNSTREAM_FACT_NAME] = variant_misses_transcript and bool(
            mask & _FACT_BIT[_DOWNSTREAM_FACT_NAME]
        )

        fact_context.prefill(served_values)
        self.facts_served += len(served_values)
