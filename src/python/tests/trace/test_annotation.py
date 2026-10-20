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
# pylint: disable=redefined-outer-name
"""Tests for the finished AnnotationTrace and how it renders."""

from __future__ import annotations

import pytest

from ensembl.vep.trace.annotation import AnnotationTrace, describe_overlap
from ensembl.vep.trace.events import FactEvent, FactOrigin, PruneEvent, RuleEvent


@pytest.fixture
def sample_trace() -> AnnotationTrace:
    """A small, hand-built trace exercising every fact origin and a pruning decision."""
    return AnnotationTrace(
        fact_events=(
            FactEvent(name="overlap_feature", value=True, cost=1, origin=FactOrigin.COMPUTED, reads=3),
            FactEvent(name="overlap_exon", value=True, cost=0, origin=FactOrigin.PREFILLED, reads=1),
            FactEvent(name="_codon_alleles", value=("ATG", "ACG"), cost=0, origin=FactOrigin.CACHED, reads=2),
            FactEvent(
                name="overlap_cds",
                value=False,
                cost=20,
                origin=FactOrigin.FORCED_FALSE,
                forced_by=("overlap_exon",),
            ),
        ),
        rule_events=(
            RuleEvent(term="missense_variant", rank=13, fired=True, settled_early=False, reason="atom"),
            RuleEvent(
                term="coding_sequence_variant", rank=23, fired=True, settled_early=False, reason="atom"
            ),
        ),
        prune_events=(PruneEvent(dropped_term="coding_sequence_variant", subsumed_by="missense_variant"),),
        emitted_terms=("missense_variant",),
        annotator_name="consequence.linear",
        overlap_label="ENST0001 / chrX:1051 A>C",
    )


class TestTraceTotals:
    """The cost figures a trace derives from its own events."""

    def test_computed_cost_counts_only_computed_facts(self, sample_trace: AnnotationTrace) -> None:
        """Only a COMPUTED fact charged anything."""
        assert sample_trace.computed_cost == 1

    def test_avoided_cost_counts_forced_false_facts(self, sample_trace: AnnotationTrace) -> None:
        """A forced-False fact carries the cost it did not charge."""
        assert sample_trace.avoided_cost == 20

    def test_events_by_origin_selects_in_trace_order(self, sample_trace: AnnotationTrace) -> None:
        """Selecting one origin preserves the order the derivation materialised facts in."""
        names = [event.name for event in sample_trace.events_by_origin(FactOrigin.COMPUTED)]
        assert names == ["overlap_feature"]

    def test_an_empty_trace_has_no_cost(self) -> None:
        """Every field defaults to empty, so an unfilled trace still answers."""
        empty = AnnotationTrace()
        assert empty.computed_cost == 0
        assert empty.avoided_cost == 0
        assert not empty.events_by_origin(FactOrigin.CACHED)


class TestDescribeOverlap:
    """describe_overlap's duck-typed formatting and graceful degradation."""

    def test_describe_overlap_names_the_feature_and_the_allele(self) -> None:
        """An overlap description pairs the feature with the allele's own name."""

        class _Allele:
            display_name = "chrX:1051:A>C"

        class _Feature:
            stable_id = "ENST0001"

        class _Overlap:
            feature = _Feature()
            overlapping_feature = _Allele()

        assert describe_overlap(_Overlap()) == "ENST0001 / chrX:1051:A>C"

    def test_missing_overlapping_feature_degrades_gracefully(self) -> None:
        """A missing allele degrades to a question mark rather than raising."""

        class _Feature:
            stable_id = "ENST0001"

        class _Overlap:
            feature = _Feature()
            overlapping_feature = None

        assert describe_overlap(_Overlap()) == "ENST0001 / ?"

    def test_bare_object_does_not_raise(self) -> None:
        """A completely unrelated object does not raise: every attribute lookup degrades."""
        assert describe_overlap(object()) is not None


class TestFormatExplanation:
    """The rendering: header, fact rows, rule rows, pruning, emitted line."""

    def test_rendering_contains_every_origin_and_the_emitted_term(
        self, sample_trace: AnnotationTrace
    ) -> None:
        """A rendered explanation names the overlap, each origin, the pruning and the result."""
        rendered = sample_trace.format_explanation()
        assert "annotator: consequence.linear" in rendered
        assert "ENST0001 / chrX:1051 A>C" in rendered
        assert "prefilled by zone map" in rendered
        assert "codon cache hit" in rendered
        assert "forced False via overlap_exon" in rendered
        assert "coding_sequence_variant <- subsumed by missense_variant" in rendered
        assert "emitted: missense_variant" in rendered

    def test_rendering_reports_read_counts(self, sample_trace: AnnotationTrace) -> None:
        """Read counts are what replaces the old per-lookup memo-hit events."""
        rendered = sample_trace.format_explanation()
        assert "6 reads" in rendered
        assert "x3" in rendered

    def test_rendering_reports_avoided_cost(self, sample_trace: AnnotationTrace) -> None:
        """The summary line says what the derivation did not have to charge."""
        assert "cost 20 avoided" in sample_trace.format_explanation()

    def test_empty_trace_reports_none_pruned_and_none_emitted(self) -> None:
        """An empty trace still renders, reporting no pruning and nothing emitted."""
        rendered = AnnotationTrace().format_explanation()
        assert "pruned: (none)" in rendered
        assert "emitted: (none)" in rendered

    def test_a_computed_fact_shows_its_cost_tier(self, sample_trace: AnnotationTrace) -> None:
        """The tier name comes from CostTier rather than a table kept here."""
        assert "BOUNDS" in sample_trace.format_explanation()
