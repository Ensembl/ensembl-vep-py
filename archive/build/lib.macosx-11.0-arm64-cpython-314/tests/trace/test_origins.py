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
"""Integration tests for fact origins recorded by consequence strategies."""

# pylint: disable=redefined-outer-name

from __future__ import annotations

import logging

import pytest

from ensembl.vep.annotators.registry import get_annotator
from ensembl.vep.caches.codon_cache import CodonEditLRUCache
from ensembl.vep.caches.zone_map import ZoneMapCache
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.trace import FactOrigin

_ALL_STRATEGIES: tuple[tuple[str, str], ...] = (
    ("consequence.linear", "ensembl.vep.annotators.consequence.strategies.linear"),
    ("consequence.staged", "ensembl.vep.annotators.consequence.strategies.staged"),
    ("consequence.columnar", "ensembl.vep.annotators.consequence.strategies.columnar"),
    ("consequence.columnar-staged", "ensembl.vep.annotators.consequence.strategies.columnar_staged"),
    ("consequence.columnar-kernels", "ensembl.vep.annotators.consequence.strategies.columnar_staged"),
)

_STAGED_STRATEGIES: tuple[tuple[str, str], ...] = (
    ("consequence.staged", "ensembl.vep.annotators.consequence.strategies.staged"),
    ("consequence.columnar-staged", "ensembl.vep.annotators.consequence.strategies.columnar_staged"),
    ("consequence.columnar-kernels", "ensembl.vep.annotators.consequence.strategies.columnar_staged"),
)

_CODING_OVERLAP_KEY = "ENST0001::ENST0001_cds3_del1"
_DEEP_INTRON_OVERLAP_KEY = "ENST0001::ENST0001_intron0_deep"


@pytest.fixture(params=_ALL_STRATEGIES, ids=[name for name, _logger_name in _ALL_STRATEGIES])
def strategy_name_and_logger(request: pytest.FixtureRequest) -> tuple[str, str]:
    """Parametrize over every consequence strategy and its tracing logger."""
    return request.param


@pytest.fixture(params=_STAGED_STRATEGIES, ids=[name for name, _logger_name in _STAGED_STRATEGIES])
def staged_strategy_name_and_logger(request: pytest.FixtureRequest) -> tuple[str, str]:
    """Parametrize over strategies that derive forced-false facts."""
    return request.param


class TestForcedFalseOrigin:
    """Tests for facts settled by staged dependency propagation."""

    def test_overlap_cds_is_forced_false_by_overlap_exon(
        self,
        staged_strategy_name_and_logger: tuple[str, str],
        boundary_overlaps: dict[str, TranscriptVariationAllele],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Test deep intronic overlaps avoid computing the CDS overlap fact."""
        annotator_name, logger_name = staged_strategy_name_and_logger
        caplog.set_level(logging.DEBUG, logger=logger_name)
        annotator = get_annotator(annotator_name)

        annotator.annotate_feature_overlap(boundary_overlaps[_DEEP_INTRON_OVERLAP_KEY])
        trace = getattr(annotator, "last_trace", None)

        assert trace is not None
        events_by_name = {event.name: event for event in trace.fact_events}
        overlap_cds = events_by_name["overlap_cds"]
        assert overlap_cds.origin is FactOrigin.FORCED_FALSE
        assert overlap_cds.forced_by == ("overlap_exon",)
        assert overlap_cds.value is False

    def test_zone_map_prefill_still_allows_downstream_forcing(
        self,
        staged_strategy_name_and_logger: tuple[str, str],
        boundary_overlaps: dict[str, TranscriptVariationAllele],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Test a prefilled false fact can force dependent facts false."""
        annotator_name, logger_name = staged_strategy_name_and_logger
        caplog.set_level(logging.DEBUG, logger=logger_name)
        annotator = get_annotator(annotator_name, zone_map=ZoneMapCache())

        annotator.annotate_feature_overlap(boundary_overlaps[_DEEP_INTRON_OVERLAP_KEY])
        trace = getattr(annotator, "last_trace", None)

        assert trace is not None
        events_by_name = {event.name: event for event in trace.fact_events}
        assert events_by_name["overlap_cds"].origin is FactOrigin.PREFILLED
        forced_events = trace.events_by_origin(FactOrigin.FORCED_FALSE)
        assert any("overlap_cds" in event.forced_by for event in forced_events)


class TestPrefilledOrigin:
    """Tests for facts served by an enabled zone map."""

    def test_overlap_exon_is_recorded_as_prefilled(
        self,
        strategy_name_and_logger: tuple[str, str],
        boundary_overlaps: dict[str, TranscriptVariationAllele],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Test every strategy attributes a zone-map fact to prefill."""
        annotator_name, logger_name = strategy_name_and_logger
        caplog.set_level(logging.DEBUG, logger=logger_name)
        annotator = get_annotator(annotator_name, zone_map=ZoneMapCache())

        annotator.annotate_feature_overlap(boundary_overlaps[_CODING_OVERLAP_KEY])
        trace = getattr(annotator, "last_trace", None)

        assert trace is not None
        overlap_exon = next(event for event in trace.fact_events if event.name == "overlap_exon")
        assert overlap_exon.origin is FactOrigin.PREFILLED
        assert overlap_exon.cost == 0


class TestCachedOrigin:
    """Tests for facts served by an enabled codon cache."""

    def test_codon_alleles_are_recorded_as_cached_on_a_hit(
        self,
        strategy_name_and_logger: tuple[str, str],
        boundary_overlaps: dict[str, TranscriptVariationAllele],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Test every strategy attributes a shared codon edit to the cache."""
        annotator_name, logger_name = strategy_name_and_logger
        caplog.set_level(logging.DEBUG, logger=logger_name)
        codon_cache = CodonEditLRUCache()
        overlap = boundary_overlaps[_CODING_OVERLAP_KEY]
        warm_up_annotator = get_annotator(annotator_name, codon_cache=codon_cache)
        hit_annotator = get_annotator(annotator_name, codon_cache=codon_cache)

        warm_up_annotator.annotate_feature_overlap(overlap)
        hit_annotator.annotate_feature_overlap(overlap)
        trace = getattr(hit_annotator, "last_trace", None)

        assert trace is not None
        codon_alleles = next(event for event in trace.fact_events if event.name == "_codon_alleles")
        assert codon_alleles.origin is FactOrigin.CACHED
