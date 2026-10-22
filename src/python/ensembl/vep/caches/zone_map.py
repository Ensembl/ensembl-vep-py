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
"""Cache transcript geometry as positional-fact zones.

Zone boundaries partition a transcript and its flanks into spans with constant positional facts.
Each zone stores a fact bitmask, allowing variant spans to be evaluated with a binary search and
short bitwise union instead of repeated interval scans. Allele-dependent sequence facts remain
outside the zone map.

``ZoneMapCache`` reuses maps across variants and evicts least-recently-used maps by estimated byte
size. Computed facts are inserted through ``FactContext.prefill`` without charging their normal
evaluation cost.
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
"""Default byte budget for a zone-map cache."""

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
"""Facts determined by a genomic span and fixed transcript geometry."""


_ZONE_FACT_NAMES: tuple[str, ...] = POSITIONAL_FACT_NAMES[:-2] + (
    "overlap_five_prime_utr",
    "overlap_three_prime_utr",
    "overlap_start_codon",
    "overlap_stop_codon",
    "overlap_cds",
)
"""Facts served as a union of the zone masks covered by a variant."""

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
"""Facts written to each overlap memo as well as returned in columns.

These facts appear in another fact's requirements and must remain available through ``get_fact``.
The remaining served facts are read directly from columns by columnar annotators.
"""

_MASK_INTERN_CACHE: dict[int, int] = {}
"""Interned zone masks shared by all zone maps.

Zone masks repeat frequently. Sharing boxed integers keeps tuple lookup fast while reducing the
memory overhead of repeated values. The cache is bounded by the masks expressible by
``_FACT_BIT``.
"""


_ESTIMATED_ENTRY_OVERHEAD_BYTES = 250
"""Fixed per-entry allowance used by ``ZoneMap.estimated_bytes``.

The allowance covers the zone-map object, cache key, and cache bookkeeping that container sizes do
not include. It is calibrated against ``tracemalloc`` measurements in the zone-map tests.
"""


_FLANK_ONLY_MARGIN = 16
"""Distance beyond a transcript where only flank facts can hold.

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
    """Partition a transcript and its flanks into positional-fact zones.

    ``starts[i]`` is the first genomic position of zone ``i`` and ``masks[i]`` its fact bitmask;
    each zone ends immediately before the next start. ``hash=False`` reflects the unhashable
    ``array`` used for starts.

    Attributes:
        lowest_position (int): Lowest genomic position covered by any zone.
        highest_position (int): Highest genomic position covered by any zone.
        starts (array[int]): Ascending zone starts packed as signed 64-bit integers.
        masks (tuple[int, ...]): Fact mask for each corresponding zone start.
        inverted (tuple[tuple[int, int, int], ...]): Degenerate ``(bit, start, end)`` intervals
            matched when a variant spans both bounds.
    """

    lowest_position: int
    highest_position: int
    starts: array[int] = field(converter=_int64_array, default=())
    masks: tuple[int, ...] = field(converter=_dedup_masks, default=())
    inverted: tuple[tuple[int, int, int], ...] = field(default=())

    def lookup(self, min_position: int, max_position: int) -> int:
        """Return the union of masks for every zone overlapped by a span.

        Args:
            min_position (int): Inclusive lower bound of the genomic span.
            max_position (int): Inclusive upper bound of the genomic span.

        Returns:
            int: Combined mask for overlapping zones and covered inverted intervals.
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
        """Estimate the heap space used by this zone map.

        Container sizes are combined with a fixed allowance for the zone-map object and cache
        bookkeeping. Interned mask integers are not charged per map.

        Returns:
            int: Estimated size in bytes.
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
        """Build the positional-fact zones for a transcript.

        Splice windows use ``intron_offsets_to_genomic_span``, matching the per-variant facts.
        Coding intervals come from the transcript's CDS, UTR, and codon mappings.

        Args:
            transcript (Transcript): Transcript to compile.
            upstream_distance (int): Bases before the 5' end that ``is_upstream`` covers.
            downstream_distance (int): Bases past the 3' end that ``is_downstream`` covers.

        Returns:
            ZoneMap: Compiled zone partition.
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


@define(kw_only=True, auto_attribs=True, eq=False)
class ZoneMapCache:
    """Store per-transcript zone maps in a byte-budget LRU cache.

    Maps are built outside the lock and inserted or evicted while holding it. ``eq=False``
    preserves identity semantics for mutable shared state.

    Attributes:
        max_bytes (int): Byte budget for cached zone maps.
        upstream_distance (int): Upstream distance covered by flank zones.
        downstream_distance (int): Downstream distance covered by flank zones.
        hit_count (int): Number of cache hits.
        miss_count (int): Number of maps built after cache misses.
        maps_built (int): Alias for ``miss_count``.
        eviction_count (int): Number of entries evicted.
        build_seconds (float): Total time spent building zone maps.
        facts_served (int): Number of fact values prefilled without evaluation.
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
    """Cache fact names implied by each zone mask.

    The mapping depends only on ``_FACT_BIT`` and is shared by all instances. It is bounded by the
    masks expressible by ``_MASK_FACT_NAMES``.
    """

    @staticmethod
    def _implied_fact_names(mask: int) -> tuple[str, ...]:
        """Return the fact names implied by a zone mask.

        Args:
            mask (int): Packed zone mask returned by ``ZoneMap.lookup``.

        Returns:
            tuple[str, ...]: Set facts with incompatible flank facts suppressed.
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
        """Return the number of zone maps built."""
        return self.miss_count

    def zone_map(self, transcript: Transcript) -> ZoneMap:
        """Return a transcript's cached zone map, building it on a miss.

        Entries are keyed by stable ID and version. Concurrent misses may build the same map more
        than once, but insertion keeps byte accounting exact.

        Args:
            transcript (Transcript): Transcript to look up.

        Returns:
            ZoneMap: Transcript zone partition.
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
        """Return zone-map facts as bit columns for a context batch.

        Facts in ``MEMO_FACT_NAMES`` are also written to each context because other facts read
        them through ``get_fact``.

        Args:
            fact_contexts (Sequence[FactContext]): Contexts in batch-index order.

        Returns:
            dict[str, int]: Fact names mapped to bit columns.
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
        """Prefill every zone-map fact without charging evaluation cost.

        Args:
            fact_context (FactContext): Overlap context to prefill.
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
