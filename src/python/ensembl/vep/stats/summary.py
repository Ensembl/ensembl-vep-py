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
"""Term-level annotation run summaries.

Most of it is already on the model: ``Variation.class_so_term`` is the variant-classes chart and
``Variation.most_severe_consequence`` the most-severe chart. The collector adds only what is
genuinely cross-variant -- every consequence, the coding ones, variants by chromosome and where
along it they fall -- and the run's own metadata. Cheap enough to be on by default: about a
microsecond per overlap, against a fact layer called ninety times as often.
"""

from __future__ import annotations

from collections.abc import Mapping

from attrs import define, field, frozen

from ensembl.vep.annotators.consequence.rules.rulesets import coding_terms
from ensembl.vep.common.constants import UNKNOWN_VALUE
from ensembl.vep.models.variation import Variation
from ensembl.vep.render.chart import Slice
from ensembl.vep.render.report import ChartKind, ReportSection

__all__ = ["POSITION_BINS", "RunSummary", "RunSummaryCollector", "RunSummaryDelta"]

POSITION_BINS = 10
"""How many equal bins a chromosome's observed span is divided into for the position chart."""

_METADATA_HEADERS = ("Key", "Value")
_CHROMOSOME_HEADERS = ("Chromosome", "Variants")
_POSITION_HEADERS = ("Chromosome", "Bin", "From", "To", "Variants")


def _counts(totals: Mapping[str, int]) -> tuple[tuple[str, int], ...]:
    """Return counts largest first, then by name, so a chart reads top-down."""
    return tuple(sorted(totals.items(), key=lambda entry: (-entry[1], entry[0])))


def _term_slices(counts: tuple[tuple[str, int], ...]) -> tuple[Slice, ...]:
    """Return term counts as slices, each taking its colour from the chart's default palette."""
    return tuple(Slice(label=term, value=count) for term, count in counts)


@frozen(kw_only=True, auto_attribs=True)
class RunSummary:
    """Term-level summary of an annotation run.

    Attributes:
        variants (int): Variations observed.
        alleles (int): Alternate alleles across them.
        overlaps (int): Feature overlaps across those alleles.
        variant_classes (tuple[tuple[str, int], ...]): Class SO term -> variants.
        most_severe (tuple[tuple[str, int], ...]): Most severe term -> variants.
        consequences (tuple[tuple[str, int], ...]): Every term -> overlaps carrying it.
        coding_consequences (tuple[tuple[str, int], ...]): The coding subset of the above.
        by_chromosome (tuple[tuple[str, int], ...]): Chromosome -> variants.
        positions (tuple[tuple[str, int, int, int, int], ...]): Chromosome, bin, from, to,
            variants -- the distribution along each chromosome's observed span. Streaming merges
            retain exact spans and totals; bin allocation is estimated if a later batch expands a
            chromosome's span.
        metadata (tuple[tuple[str, str], ...]): The run's own facts: command, counts, timings.
    """

    variants: int = field(default=0)
    alleles: int = field(default=0)
    overlaps: int = field(default=0)
    variant_classes: tuple[tuple[str, int], ...] = field(default=())
    most_severe: tuple[tuple[str, int], ...] = field(default=())
    consequences: tuple[tuple[str, int], ...] = field(default=())
    coding_consequences: tuple[tuple[str, int], ...] = field(default=())
    by_chromosome: tuple[tuple[str, int], ...] = field(default=())
    positions: tuple[tuple[str, int, int, int, int], ...] = field(default=())
    metadata: tuple[tuple[str, str], ...] = field(default=())

    @property
    def title(self) -> str:
        """Return the document title."""
        return "Run summary"

    @property
    def subtitle(self) -> str:
        """Return the headline counts."""
        return f"{self.variants} variants, {self.alleles} alleles, {self.overlaps} overlaps"

    def sections(self) -> tuple[ReportSection, ...]:
        """Return the summary as report sections, in the order the Perl summary shows them."""
        return (
            ReportSection(
                identifier="metadata", title="Run", headers=_METADATA_HEADERS, rows=tuple(self.metadata)
            ),
            ReportSection(
                identifier="classes",
                title="Variant classes",
                slices=tuple(Slice(label=name, value=count) for name, count in self.variant_classes),
                chart=ChartKind.PIE,
            ),
            ReportSection(
                identifier="most_severe",
                title="Most severe consequences",
                slices=_term_slices(self.most_severe),
                chart=ChartKind.PIE,
            ),
            ReportSection(
                identifier="consequences",
                title="Consequences (all)",
                slices=_term_slices(self.consequences),
                chart=ChartKind.BARS,
            ),
            ReportSection(
                identifier="coding",
                title="Coding consequences",
                slices=_term_slices(self.coding_consequences),
                chart=ChartKind.BARS,
            ),
            ReportSection(
                identifier="chromosomes",
                title="Variants by chromosome",
                slices=tuple(Slice(label=name, value=count) for name, count in self.by_chromosome),
                chart=ChartKind.BARS,
            ),
            ReportSection(
                identifier="positions",
                title="Position in chromosome",
                headers=_POSITION_HEADERS,
                rows=tuple(self.positions),
            ),
        )


@frozen(kw_only=True, auto_attribs=True)
class RunSummaryDelta:
    """Mergeable summary counts from one worker buffer.

    Attributes:
        variants (int): Number of variants.
        alleles (int): Number of alleles.
        overlaps (int): Number of feature overlaps.
        classes (tuple[tuple[str, int], ...]): Counts by variant class.
        most_severe (tuple[tuple[str, int], ...]): Counts by most severe consequence.
        consequences (tuple[tuple[str, int], ...]): Counts by consequence term.
        chromosomes (tuple[tuple[str, int], ...]): Counts by sequence region.
        positions (tuple[tuple[str, int], ...]): Observed positions by sequence region.
    """

    variants: int = 0
    alleles: int = 0
    overlaps: int = 0
    classes: tuple[tuple[str, int], ...] = field(default=())
    most_severe: tuple[tuple[str, int], ...] = field(default=())
    consequences: tuple[tuple[str, int], ...] = field(default=())
    chromosomes: tuple[tuple[str, int], ...] = field(default=())
    positions: tuple[tuple[str, int], ...] = field(default=())


@define(kw_only=True, auto_attribs=True)
class _StreamingPositionBins:
    """Bounded position distribution for one recurrent sequence region.

    Reproject existing counts by bin midpoint when the observed span expands. This keeps the
    region's endpoints and total exact while bounding storage to ``POSITION_BINS`` counters.
    """

    low: int | None = None
    high: int | None = None
    counts: list[int] = field(factory=lambda: [0] * POSITION_BINS)

    def add(self, positions: tuple[int, ...]) -> None:
        """Add exact positions to the current bounded distribution."""
        if not positions:
            return

        new_low = min(positions) if self.low is None else min(self.low, *positions)
        new_high = max(positions) if self.high is None else max(self.high, *positions)
        if self.low is not None and self.high is not None:
            self._resize(new_low, new_high)
        self.low = new_low
        self.high = new_high
        for position in positions:
            self.counts[self._index(position)] += 1

    def rows(self, chromosome: str) -> tuple[tuple[str, int, int, int, int], ...]:
        """Return up to ten non-empty coordinate ranges for this distribution."""
        assert self.low is not None
        assert self.high is not None
        width = self._width(self.low, self.high)
        span = self.high - self.low + 1
        bin_count = min(POSITION_BINS, (span + width - 1) // width)
        rows = []
        for index, count in enumerate(self.counts[:bin_count]):
            bin_low = self.low + index * width
            rows.append((chromosome, index + 1, bin_low, min(self.high, bin_low + width - 1), count))
        return tuple(rows)

    def _resize(self, low: int, high: int) -> None:
        """Reproject retained counts if new positions expand the observed span."""
        assert self.low is not None
        assert self.high is not None
        if low == self.low and high == self.high:
            return

        old_low = self.low
        old_high = self.high
        old_width = self._width(old_low, old_high)
        old_counts = self.counts
        self.low = low
        self.high = high
        self.counts = [0] * POSITION_BINS
        for index, count in enumerate(old_counts):
            if not count:
                continue
            bin_low = old_low + index * old_width
            bin_high = min(old_high, bin_low + old_width - 1)
            midpoint = bin_low + (bin_high - bin_low) // 2
            self.counts[self._index(midpoint)] += count

    def _index(self, position: int) -> int:
        """Return the bin index for a position in the current span."""
        assert self.low is not None
        assert self.high is not None
        width = self._width(self.low, self.high)
        return min((position - self.low) // width, POSITION_BINS - 1)

    @staticmethod
    def _width(low: int, high: int) -> int:
        """Return the integral width of one report bin."""
        return max(1, (high - low + 1 + POSITION_BINS - 1) // POSITION_BINS)


@define(kw_only=True, auto_attribs=True)
class RunSummaryCollector:
    """Incremental collector for a ``RunSummary``.

    Attributes:
        metadata (dict[str, str]): The run's own facts, written into the summary as given.
    """

    metadata: dict[str, str] = field(factory=dict)
    _variants: int = field(default=0, init=False)
    _alleles: int = field(default=0, init=False)
    _overlaps: int = field(default=0, init=False)
    _classes: dict[str, int] = field(factory=dict, init=False)
    _most_severe: dict[str, int] = field(factory=dict, init=False)
    _consequences: dict[str, int] = field(factory=dict, init=False)
    _by_chromosome: dict[str, int] = field(factory=dict, init=False)
    _positions: dict[str, list[int]] = field(factory=dict, init=False)
    _streaming_positions: dict[str, _StreamingPositionBins] = field(factory=dict, init=False)

    def observe(self, variation: Variation) -> None:
        """Count one variation and everything hanging off it.

        Args:
            variation (Variation): The variation, with its alleles' consequences already known --
                whether assigned by an annotator or read back from a file.
        """
        uploaded_locus = variation.uploaded_locus
        uploaded_name = uploaded_locus.seq_region_name if uploaded_locus is not None else None
        chromosome = uploaded_name or UNKNOWN_VALUE
        self._assert_not_ordered_region(chromosome)
        self._variants += 1
        self._by_chromosome[chromosome] = self._by_chromosome.get(chromosome, 0) + 1
        self._positions.setdefault(chromosome, []).append(variation.start)

        class_term = variation.class_so_term
        class_name = class_term.value if class_term is not None else "unknown"
        self._classes[class_name] = self._classes.get(class_name, 0) + 1

        most_severe = variation.most_severe_consequence
        if most_severe is not None:
            self._most_severe[most_severe.term] = self._most_severe.get(most_severe.term, 0) + 1

        for allele in variation.variation_alleles:
            self._alleles += 1
            self._overlaps += len(allele.transcript_variation_alleles)
            for consequence in allele.overlap_consequences:
                self._consequences[consequence.term] = self._consequences.get(consequence.term, 0) + 1

    def merge(self, delta: RunSummaryDelta) -> None:
        """Merge a worker-buffer summary delta."""

        for chromosome, _position in delta.positions:
            self._assert_not_ordered_region(chromosome)
        self._merge_counts(delta)
        for chromosome, position in delta.positions:
            self._positions.setdefault(chromosome, []).append(position)

    def merge_ordered_delta(self, delta: RunSummaryDelta) -> None:
        """Merge one sequence-region worker delta in arbitrary source order.

        Every position in a delta must belong to one sequence region. Deltas for
        a region may recur after other regions. A bounded histogram per region
        releases every raw position after merging. Expanding a region's observed
        span reprojects its existing counts by bin midpoint, so endpoints and
        totals stay exact while individual bin counts become an estimate.

        Raises:
            ValueError: If a delta mixes sequence regions or the region has unordered positions.
        """

        regions = {chromosome for chromosome, _position in delta.positions}
        if not regions:
            self._merge_counts(delta)
            return
        if len(regions) != 1:
            raise ValueError("an ordered summary delta must contain one sequence region")

        region = next(iter(regions))
        if region in self._positions:
            raise ValueError(f"ordered summary region {region!r} already has unordered positions")
        self._merge_counts(delta)
        positions = tuple(position for _chromosome, position in delta.positions)
        self._streaming_positions.setdefault(region, _StreamingPositionBins()).add(positions)

    @property
    def completed_position_regions(self) -> tuple[str, ...]:
        """Return regions whose raw positions have been released into bounded bins."""

        return tuple(sorted(self._streaming_positions))

    @property
    def retained_ordered_position_count(self) -> int:
        """Return the number of raw positions retained by streaming merges."""

        return 0

    def _merge_counts(self, delta: RunSummaryDelta) -> None:
        """Merge non-positional counts from a worker delta."""

        self._variants += delta.variants
        self._alleles += delta.alleles
        self._overlaps += delta.overlaps
        for totals, values in (
            (self._classes, delta.classes),
            (self._most_severe, delta.most_severe),
            (self._consequences, delta.consequences),
            (self._by_chromosome, delta.chromosomes),
        ):
            for name, count in values:
                totals[name] = totals.get(name, 0) + count

    def delta(self) -> RunSummaryDelta:
        """Return the mergeable counts accumulated so far."""
        positions: list[tuple[str, int]] = []
        for chromosome, values in self._positions.items():
            for position in values:
                positions.append((chromosome, position))
        return RunSummaryDelta(
            variants=self._variants,
            alleles=self._alleles,
            overlaps=self._overlaps,
            classes=tuple(self._classes.items()),
            most_severe=tuple(self._most_severe.items()),
            consequences=tuple(self._consequences.items()),
            chromosomes=tuple(self._by_chromosome.items()),
            positions=tuple(positions),
        )

    def finish(self) -> RunSummary:
        """Return the summary of everything observed so far."""
        coding = coding_terms()
        consequences = _counts(self._consequences)
        return RunSummary(
            variants=self._variants,
            alleles=self._alleles,
            overlaps=self._overlaps,
            variant_classes=_counts(self._classes),
            most_severe=_counts(self._most_severe),
            consequences=consequences,
            coding_consequences=tuple((term, count) for term, count in consequences if term in coding),
            by_chromosome=_counts(self._by_chromosome),
            positions=self._all_binned_positions(),
            metadata=tuple(self.metadata.items()),
        )

    def _binned_positions(self) -> tuple[tuple[str, int, int, int, int], ...]:
        """Return each chromosome's variants binned along its observed span."""
        rows: list[tuple[str, int, int, int, int]] = []
        for chromosome in sorted(self._positions):
            positions = self._positions[chromosome]
            distribution = _StreamingPositionBins()
            distribution.add(tuple(positions))
            rows.extend(distribution.rows(chromosome))
        return tuple(rows)

    def _all_binned_positions(self) -> tuple[tuple[str, int, int, int, int], ...]:
        """Return unordered and streaming position rows in report order."""

        rows = list(self._binned_positions())
        for chromosome, distribution in self._streaming_positions.items():
            rows.extend(distribution.rows(chromosome))
        return tuple(sorted(rows, key=lambda row: (row[0], row[1])))

    def _assert_not_ordered_region(self, chromosome: str) -> None:
        """Reject unordered positions that would invalidate streaming bins."""

        if chromosome in self._streaming_positions:
            raise ValueError(f"summary region {chromosome!r} is already managed by ordered deltas")
