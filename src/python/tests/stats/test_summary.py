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
"""Tests for the term-level run summary."""

from __future__ import annotations

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import
import pytest

# ^ keeps the annotators-before-caches import order.
from ensembl.vep.common import CoordSystemLevel, Locus, Strand
from ensembl.vep.models.overlap import OverlapConsequence
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.models.variation import Variation
from ensembl.vep.stats.summary import POSITION_BINS, RunSummaryCollector, RunSummaryDelta


def _variation(chromosome: str, position: int, alleles: str, *terms: tuple[str, int]) -> Variation:
    """Return a variation whose single allele carries the given (term, rank) consequences."""
    built = Variation(
        start=position,
        end=position + len(alleles.split("/")[0]) - 1,
        strand=Strand.FORWARD,
        slice=Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name=chromosome,
            start=1,
            end=10**6,
            strand=Strand.FORWARD,
        ),
        allele_string=alleles,
        uploaded_locus=Locus(seq_region_name=chromosome, start=position, end=position),
    )
    consequences = [OverlapConsequence(term=term, accession="", rank=rank) for term, rank in terms]
    built.variation_alleles[0].seed_attrs(overlap_consequences=consequences)
    return built


class TestCollector:
    """What one pass over some variations accumulates."""

    def test_counts_variants_alleles_and_classes(self) -> None:
        """The class chart reads straight off Variation.class_so_term."""
        collector = RunSummaryCollector()
        collector.observe(_variation("1", 100, "A/T", ("missense_variant", 13)))
        collector.observe(_variation("1", 200, "AC/A", ("frameshift_variant", 5)))
        summary = collector.finish()
        assert summary.variants == 2
        assert summary.alleles == 2
        assert dict(summary.variant_classes) == {"SNV": 1, "indel": 1}

    def test_most_severe_is_the_lowest_rank(self) -> None:
        """Most severe comes off Variation.most_severe_consequence, rank 1 first."""
        collector = RunSummaryCollector()
        collector.observe(_variation("1", 100, "A/T", ("intron_variant", 28), ("stop_gained", 4)))
        summary = collector.finish()
        assert summary.most_severe == (("stop_gained", 1),)
        assert dict(summary.consequences) == {"intron_variant": 1, "stop_gained": 1}

    def test_coding_consequences_are_the_catalogue_derived_subset(self) -> None:
        """Coding means the rule reads overlap_cds; intron_variant does not."""
        collector = RunSummaryCollector()
        collector.observe(_variation("1", 100, "A/T", ("intron_variant", 28), ("missense_variant", 13)))
        summary = collector.finish()
        assert dict(summary.coding_consequences) == {"missense_variant": 1}

    def test_variants_are_counted_per_chromosome_and_binned_along_it(self) -> None:
        """The positional distribution spans each chromosome's observed range."""
        collector = RunSummaryCollector()
        for position in (100, 200, 300, 1000):
            collector.observe(_variation("2", position, "A/T"))
        collector.observe(_variation("X", 5, "A/T"))
        summary = collector.finish()
        assert dict(summary.by_chromosome) == {"2": 4, "X": 1}
        rows_for_2 = [row for row in summary.positions if row[0] == "2"]
        assert len(rows_for_2) == POSITION_BINS
        assert sum(row[4] for row in rows_for_2) == 4
        assert rows_for_2[0][2] == 100 and rows_for_2[-1][3] == 1000

    def test_metadata_is_carried_through(self) -> None:
        """The run's own facts appear as given."""
        summary = RunSummaryCollector(metadata={"input": "x.vcf"}).finish()
        assert summary.metadata == (("input", "x.vcf"),)
        assert summary.sections()[0].rows == (("input", "x.vcf"),)

    def test_an_unannotated_variation_still_counts(self) -> None:
        """No consequences means no most-severe entry, not an error."""
        collector = RunSummaryCollector()
        collector.observe(_variation("1", 100, "A/T"))
        summary = collector.finish()
        assert summary.variants == 1
        assert summary.most_severe == ()

    def test_the_uploaded_region_name_is_what_is_counted(self) -> None:
        """Let run summaries follow the input's own spelling, not the reference's."""
        collector = RunSummaryCollector()
        variation = _variation("1", 100, "A/T")
        variation.uploaded_locus = Locus(seq_region_name="chr1", start=100, end=100)

        collector.observe(variation)

        summary = collector.finish()
        assert summary.by_chromosome == (("chr1", 1),)
        assert {row[0] for row in summary.positions} == {"chr1"}


class TestOrderedDeltaMerge:
    """Ordered worker-buffer position merging."""

    def test_matches_the_unordered_summary_across_sequence_regions(self) -> None:
        """Preserve exact position rows while buffers arrive in region order."""

        first = _delta(_variation("1", 100, "A/T"), _variation("1", 200, "A/T"))
        second = _delta(_variation("1", 300, "A/T"))
        third = _delta(_variation("2", 50, "AC/A"), _variation("2", 900, "A/T"))

        unordered = RunSummaryCollector()
        ordered = RunSummaryCollector()
        for delta in (first, second, third):
            unordered.merge(delta)
            ordered.merge_ordered_delta(delta)

        assert ordered.finish() == unordered.finish()

    def test_releases_all_merged_region_positions(self) -> None:
        """Keep bounded bin counts instead of raw positions for every region."""

        collector = RunSummaryCollector()
        collector.merge_ordered_delta(_delta(_variation("1", 100, "A/T"), _variation("1", 200, "A/T")))
        collector.merge_ordered_delta(_delta(_variation("2", 50, "A/T")))

        assert collector.completed_position_regions == ("1", "2")
        assert collector.retained_ordered_position_count == 0

    def test_rejects_a_delta_containing_multiple_regions(self) -> None:
        """Require one sequence region per ordered worker payload."""

        collector = RunSummaryCollector()
        delta = _delta(_variation("1", 100, "A/T"), _variation("2", 200, "A/T"))

        with pytest.raises(ValueError, match="one sequence region"):
            collector.merge_ordered_delta(delta)

        assert collector.finish().variants == 0

    def test_accepts_a_region_reappearing_later_in_source_order(self) -> None:
        """A valid unsorted VCF may contain region order A, B, A."""

        collector = RunSummaryCollector()
        collector.merge_ordered_delta(_delta(_variation("1", 100, "A/T")))
        collector.merge_ordered_delta(_delta(_variation("2", 50, "A/T")))
        collector.merge_ordered_delta(_delta(_variation("1", 200, "A/T")))

        summary = collector.finish()
        rows_for_1 = [row for row in summary.positions if row[0] == "1"]
        assert summary.variants == 3
        assert dict(summary.by_chromosome) == {"1": 2, "2": 1}
        assert rows_for_1[0][2] == 100
        assert rows_for_1[-1][3] == 200
        assert sum(row[4] for row in rows_for_1) == 2

    @pytest.mark.parametrize(
        "positions",
        [
            pytest.param((100,), id="single_position"),
            pytest.param((100, 103, 107), id="sub_ten_base_span"),
            pytest.param((100, 110), id="eleven_base_span"),
            pytest.param((100, 111), id="twelve_base_span"),
            pytest.param((100, 120), id="twenty_one_base_span"),
        ],
    )
    def test_short_spans_have_valid_position_bin_ranges(self, positions: tuple[int, ...]) -> None:
        """Do not emit synthetic ranges whose start lies beyond their end."""
        collector = RunSummaryCollector()
        for position in positions:
            collector.merge_ordered_delta(_delta(_variation("1", position, "A/T")))

        rows = collector.finish().positions

        assert rows[0][2] == min(positions)
        assert rows[-1][3] == max(positions)
        assert all(low <= high for _chromosome, _bin, low, high, _count in rows)
        assert sum(row[4] for row in rows) == len(positions)


def _delta(*variations: Variation) -> RunSummaryDelta:
    """Return a worker-style delta for variations in one buffer."""

    collector = RunSummaryCollector()
    for variation in variations:
        collector.observe(variation)
    return collector.delta()
