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
"""Genomic <-> CDS <-> protein coordinate mapping for a coding region of a transcript.

Built from a transcript's CDS genomic spans and strand only. Stays reusable outside a
``Transcript`` entirely.

Always builds a small (O(exon count)) bisect index at construction:
- prefix sums over the spans, walked in transcript order (reversed on the reverse strand)
- the span in genomic ascending order and related prefix sums.

Can optionally use a fully-expanded offset<->position index for O(1) lookups, at the cost of higher
memory usage. Measured ~15MB for a titin-scale (109 kb, 363-exon) CDS, against ~52 KiB for using bisect. A
byte-budget LRU ``caches/mapper_cache.py`` can limit memory usage for the index.
"""

from __future__ import annotations

import bisect
import logging

from attrs import define, field

from ensembl.vep.common import Strand

logger = logging.getLogger(__name__)

__all__ = ["Mapper", "MapperIndex"]


@define(kw_only=True, auto_attribs=True)
class MapperIndex:
    """Index holding the offset<->genomic position mapping.

    Turns O(log exon count) bisect lookup of Mapper into an O(1) lookup, at the cost of holding this
    O(length) structure in memory.

    Attributes:
        offset_to_position (list[int]): 0-based CDS offset -> Genomic position, index ``i`` is the
            genomic position of CDS offset ``i``.
        position_to_offset (dict[int, int]): Genomic position -> 0-based CDS offset, inverse of
            ``offset_to_position``.
    """

    offset_to_position: list[int]
    position_to_offset: dict[int, int]


@define(kw_only=True, auto_attribs=True)
class Mapper:
    """Genomic <-> offset coordinate mapping over an ordered list of genomic spans.

    Deliberately generic: the spans can be a transcript's CDS blocks (``Transcript.cds_mapper``),
    its whole spliced exon set (``Transcript.cdna_mapper``), or any other ordered block list.
    "Offset" therefore means whatever the caller built the mapper from, CDS-relative in the first
    case, cDNA-relative in the second.

    The exception is the codon and protein methods (``codon_window``, ``codon_genomic_spans``,
    ``offset_to_protein_offset``, ``protein_offset_to_span``). Those divide offsets by three and so
    only mean anything on a mapper built over CDS spans; calling them on a cDNA mapper is a
    caller error, not something this class can detect.

    Attributes:
        genomic_spans (list[tuple[int, int]]): The genomic spans, in transcript (5' to 3') order.
            Index 0 is always the first span. Genomic-descending on the reverse strand.
        strand (Strand): The owning transcript's strand.
        coding_start_offset (int): 0-based offset to the first coding base
        coding_end_offset (int): 0-based offset to the last coding base
    """

    genomic_spans: list[tuple[int, int]]
    strand: Strand
    coding_start_offset: int
    coding_end_offset: int
    _prefix_lengths: list[int] = field(init=False, repr=False)
    _ordered_genomic_spans: list[tuple[int, int]] = field(init=False, repr=False)
    _ordered_span_starts: list[int] = field(init=False, repr=False)
    _ordered_prefix_lengths: list[int] = field(init=False, repr=False)

    def __attrs_post_init__(self) -> None:
        """Build both bisect views: transcript order (for offset->genomic) and genomic order (for
        genomic->offset)."""
        prefix_lengths = [0]
        for span_start, span_end in self.genomic_spans:
            prefix_lengths.append(prefix_lengths[-1] + (span_end - span_start + 1))
        self._prefix_lengths = prefix_lengths

        # Order by genomic start position ascending
        genomic_order = sorted(range(len(self.genomic_spans)), key=lambda index: self.genomic_spans[index][0])

        ordered_genomic_spans: list[tuple[int, int]] = []
        ordered_span_starts: list[int] = []
        ordered_prefix_lengths: list[int] = []
        for index in genomic_order:
            span = self.genomic_spans[index]
            ordered_genomic_spans.append(span)
            ordered_span_starts.append(span[0])
            ordered_prefix_lengths.append(prefix_lengths[index])
        self._ordered_genomic_spans = ordered_genomic_spans
        self._ordered_span_starts = ordered_span_starts
        self._ordered_prefix_lengths = ordered_prefix_lengths

    @property
    def length(self) -> int:
        """Total number of nucleotides covered by ``genomic_spans``.

        Returns:
            int: Length in bases (0 if ``genomic_spans`` is empty). CDS length on a CDS mapper,
                cDNA length on a cDNA mapper.
        """
        return self._prefix_lengths[-1]

    def genomic_to_offset(
        self, genomic_position: int, *, accelerator: MapperIndex | None = None
    ) -> int | None:
        """Map a genomic position to its 0-based CDS/cDNA offset.

        Args:
            genomic_position (int): A genomic position, 1-based.
            accelerator (MapperIndex | None): If given and it covers this CDS/cDNA, answers in O(1)
                via dict lookup; otherwise falls back to an O(log exon count) bisect search.

        Returns:
            int | None: The 0-based offset into the CDS/cDNA, or None if the position is not in the CDS/cDNA.
        """
        if self.strand == Strand.UNKNOWN:
            logger.warning(
                "Strand is unknown, cannot map genomic positions with confidence. "
                "Assuming the strand is forward."
            )

        if accelerator is not None:
            return accelerator.position_to_offset.get(genomic_position)

        span_index = bisect.bisect_right(self._ordered_span_starts, genomic_position) - 1
        if span_index < 0:
            return None

        span_start, span_end = self._ordered_genomic_spans[span_index]
        if not span_start <= genomic_position <= span_end:
            return None

        offset_in_span = (
            (span_end - genomic_position)
            if self.strand == Strand.REVERSE
            else (genomic_position - span_start)
        )
        return self._ordered_prefix_lengths[span_index] + offset_in_span

    def offset_to_genomic(self, offset: int, *, accelerator: MapperIndex | None = None) -> int | None:
        """Map a 0-based CDS/cDNA offset to its genomic position.

        Args:
            offset (int): 0-based offset into the CDS/cDNA.
            accelerator (MapperIndex | None): If given, answers in O(1) via list index; otherwise
                falls back to an O(log exon count) bisect search.

        Returns:
            int | None: The genomic position, or None if the offset is out of range.
        """
        if self.strand == Strand.UNKNOWN:
            logger.warning(
                "Strand is unknown, cannot map genomic positions with confidence. "
                "Assuming the strand is forward."
            )

        if accelerator is not None:
            index = accelerator.offset_to_position
            return index[offset] if 0 <= offset < len(index) else None

        if not 0 <= offset < self.length:
            return None

        span_index = bisect.bisect_right(self._prefix_lengths, offset) - 1
        span_start, span_end = self.genomic_spans[span_index]
        offset_in_span = offset - self._prefix_lengths[span_index]
        return (span_end - offset_in_span) if self.strand == Strand.REVERSE else (span_start + offset_in_span)

    def genomic_span_to_offsets(self, min_position: int, max_position: int) -> tuple[int, int] | None:
        """Map a genomic span ``[min_position, max_position]`` to CDS/cDNA offsets ``[low, high]``.

        Equivalent to calling ``genomic_to_offset`` on every genomic position in the span and
        taking the min/max of whatever is overlapping the CDS/cDNA. But instead of calling for every genomic
        position in the span, this uses a bisect search to find the overlapping coding spans and then
        calculates the offsets from the span's prefix sums and the portion overlapping with the span.

        Suitable for a variant's genomic span, which may be very large if it is structural.

        Args:
            min_position (int): Lower bound (inclusive) of the genomic span, 1-based.
            max_position (int): Upper bound (inclusive) of the genomic span, 1-based.

        Returns:
            tuple[int, int] | None: (low, high) 0-based CDS/cDNA offsets, or None if the span touches
                no CDS/cDNA spans at all.
        """
        if self.strand == Strand.UNKNOWN:
            logger.warning(
                "Strand is unknown, cannot map genomic positions with confidence. "
                "Assuming the strand is forward."
            )

        if max_position < min_position:
            return None

        last_index = bisect.bisect_right(self._ordered_span_starts, max_position) - 1
        if last_index < 0:
            return None

        # The candidate span whose end <= min_position ended before min_position means the genomic span starts
        # right after it), in which case the first span that can possibly overlap is the next one.
        first_index = bisect.bisect_right(self._ordered_span_starts, min_position) - 1
        if first_index < 0 or self._ordered_genomic_spans[first_index][1] < min_position:
            first_index += 1

        if first_index > last_index:
            return None

        # Get the offsets from the first and last overlapping CDS/cDNA spans
        offsets: list[int] = []
        # dict.fromkeys, not a set: it dedupes the same way when first_index == last_index but
        # iterates in a defined order rather than hash order.
        for index in dict.fromkeys((first_index, last_index)):
            span_start, span_end = self._ordered_genomic_spans[index]
            overlap_start = max(span_start, min_position)
            overlap_end = min(span_end, max_position)
            prefix_length = self._ordered_prefix_lengths[index]
            if self.strand == Strand.REVERSE:
                offsets.append(prefix_length + (span_end - overlap_end))
                offsets.append(prefix_length + (span_end - overlap_start))
            else:
                offsets.append(prefix_length + (overlap_start - span_start))
                offsets.append(prefix_length + (overlap_end - span_start))

        if not offsets:
            return None
        return min(offsets), max(offsets)

    def codon_genomic_spans(self, offset: int) -> list[tuple[int, int]] | None:
        """Genomic spans of the whole codon containing ``offset``.

        A single ``(min, max)`` envelope silently excludes an intron the codon straddles. Therefore the
        three individual genomic positions are grouped into maximal runs of consecutive positions, so a
        codon split across an exon-exon junction comes back as two (or in an extreme case, three) disjoint
        spans.

        Args:
            offset (int): Any 0-based CDS offset within the codon.

        Returns:
            list[tuple[int, int]] | None: Ascending genomic spans, or None if none of the three bases
                map to a genomic position.
        """
        if offset < self.coding_start_offset or offset > self.coding_end_offset:
            return None

        # Get the offset of the first base of the codon containing offset
        codon_start = max((offset // 3) * 3, self.coding_start_offset)

        # Then determine the genomic positions of all three bases of that codon
        genomic_positions = sorted(
            position
            for position in (self.offset_to_genomic(codon_start + base_index) for base_index in range(3))
            if position is not None
        )
        if not genomic_positions:
            return None

        # Group the three genomic positions into minimal number of spans
        spans: list[tuple[int, int]] = []
        run_start = run_end = genomic_positions[0]
        for position in genomic_positions[1:]:
            if position == run_end + 1:
                run_end = position
            else:
                spans.append((run_start, run_end))
                run_start = run_end = position
        spans.append((run_start, run_end))
        return spans

    def codon_window(self, offsets: tuple[int, int] | None) -> tuple[int, int] | None:
        """0-based CDS half-open offsets ``[start, end)`` extending to whole codon(s) window.

        Args:
            offsets (tuple[int, int] | None): ``(low, high)`` 0-based CDS offsets.

        Returns:
            tuple[int, int] | None: The codon-aligned window, or None if ``offsets`` is None or
                the window is empty (e.g. the offsets range is entirely past ``length``).
        """
        if offsets is None or offsets[0] > self.coding_end_offset or offsets[1] < self.coding_start_offset:
            return None
        low_offset, high_offset = offsets
        codon_start = max((low_offset // 3) * 3, self.coding_start_offset)
        codon_end = min((high_offset // 3) * 3 + 3, self.coding_end_offset + 1)
        return (codon_start, codon_end) if codon_end > codon_start else None

    def offset_to_protein_offset(self, offset: int) -> int | None:
        """Map a 0-based CDS offset to its 0-based protein (residue) offset.

        Args:
            offset (int): 0-based offset into the CDS.

        Returns:
            int | None: 0-based residue index, or None if ``offset`` is out of range.
        """
        if not 0 <= offset < self.length:
            return None
        return offset // 3

    def protein_offset_to_span(self, protein_offset: int) -> tuple[int, int] | None:
        """Map a 0-based protein (residue) offset to its codon's CDS-offset span.

        Args:
            protein_offset (int): 0-based residue index.

        Returns:
            tuple[int, int] | None: Inclusive 0-based ``(low, high)`` CDS-offset span of the
                residue's codon (clipped to ``length`` for a partial final codon), or None if
                ``protein_offset`` is negative or starts beyond the CDS.
        """
        low = protein_offset * 3
        if protein_offset < 0 or low >= self.length:
            return None
        return low, min(low + 2, self.length - 1)

    def build_index(self) -> MapperIndex:
        """Build the O(length) accelerator.

        Returns:
            MapperIndex: The fully-expanded offset<->genomic position structures.
        """
        genomic_positions: list[int] = []
        for span_start, span_end in self.genomic_spans:
            span_positions = (
                range(span_end, span_start - 1, -1)
                if self.strand == Strand.REVERSE
                else range(span_start, span_end + 1)
            )
            genomic_positions.extend(span_positions)
        position_to_offset = {position: offset for offset, position in enumerate(genomic_positions)}
        return MapperIndex(offset_to_position=genomic_positions, position_to_offset=position_to_offset)
