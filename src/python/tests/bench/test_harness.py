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
"""Tests for the conformance gate and the timing harness."""

from __future__ import annotations

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the annotators and the facts they read.
from ensembl.vep.annotators.consequence.facts.base import REGISTRY
from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.bench.harness import (
    REFERENCE_STRATEGY,
    Timing,
    TimingPolicy,
    annotator_extras,
    conformance_gate,
    consequence_strategies,
    fact_cost,
    time_best,
)
from ensembl.vep.corpus.synthetic.default import default_corpus
from ensembl.vep.models.overlap import OverlapConsequence, TranscriptVariationAllele


class _SilentAnnotator:
    """Emits nothing for any overlap: wrong wherever the reference emits something."""

    def annotate_feature_overlap(self, _overlap: TranscriptVariationAllele) -> list[OverlapConsequence]:
        """Return no consequences."""
        return []

    def annotate_feature_overlaps(
        self, overlaps: list[TranscriptVariationAllele]
    ) -> list[list[OverlapConsequence]]:
        """Return no consequences for every overlap."""
        return [[] for _ in overlaps]


@pytest.fixture(name="overlaps", scope="module")
def fixture_overlaps() -> list[TranscriptVariationAllele]:
    """Forty boundary overlaps."""
    return [overlap for _name, overlap in default_corpus(fuzz_count=0).named_overlaps[:40]]


class TestConformanceGate:
    """The gate that runs before any timing."""

    def test_agreeing_candidates_pass(self, overlaps: list[TranscriptVariationAllele]) -> None:
        """A conforming strategy is reported as such, with what was checked."""
        report = conformance_gate(
            overlaps,
            reference=build_annotator(annotator=REFERENCE_STRATEGY),
            candidates={"consequence.staged": build_annotator(annotator="consequence.staged")},
        )
        assert report.ok
        assert report.overlaps == len(overlaps)
        assert report.candidates == ("consequence.staged",)

    def test_a_wrong_candidate_is_caught_and_named(self, overlaps: list[TranscriptVariationAllele]) -> None:
        """The mismatch names the strategy, the reference and the feature."""
        report = conformance_gate(
            overlaps,
            reference=build_annotator(annotator=REFERENCE_STRATEGY),
            candidates={"silent": _SilentAnnotator()},  # type: ignore[dict-item]
        )
        assert not report.ok
        assert report.mismatches[0].startswith(f"silent disagrees with {REFERENCE_STRATEGY} on ")

    def test_the_gate_stops_at_its_limit(self, overlaps: list[TranscriptVariationAllele]) -> None:
        """One broken strategy cannot flood the report."""
        report = conformance_gate(
            overlaps,
            reference=build_annotator(annotator=REFERENCE_STRATEGY),
            candidates={"silent": _SilentAnnotator()},  # type: ignore[dict-item]
            limit=3,
        )
        assert len(report.mismatches) == 3
        assert report.truncated

    def test_the_reference_strategy_is_listed_first(self) -> None:
        """So a speed-up column always has its reference in the first row."""
        strategies = consequence_strategies()
        assert strategies[0] == REFERENCE_STRATEGY
        assert all(name.startswith("consequence.") for name in strategies)


class TestTiming:
    """Warmup, repeats, and the best sample."""

    def test_warmup_passes_are_run_and_not_timed(self) -> None:
        """Warmup plus repeats calls, repeats samples."""
        calls = []
        timing = time_best(lambda: calls.append(1), TimingPolicy(repeats=3, warmup=2))
        assert len(calls) == 5
        assert len(timing.samples) == 3
        assert timing.best == min(timing.samples)

    def test_ms_per_thousand_contexts(self) -> None:
        """Half a second over 500 contexts is a millisecond per context, a second per thousand."""
        assert Timing(samples=(0.5,)).ms_per_thousand(500) == pytest.approx(1000.0)


class TestFactCost:
    """Counting what an annotator computes."""

    def test_it_counts_and_restores_the_registry(self, overlaps: list[TranscriptVariationAllele]) -> None:
        """Facts and cost come back positive, and every compute is its original afterwards."""
        before = {name: registered.compute for name, registered in REGISTRY.items()}
        facts_per_context, cost_per_context = fact_cost(
            build_annotator(annotator=REFERENCE_STRATEGY), overlaps
        )
        assert facts_per_context > 0
        assert cost_per_context > 0
        assert {name: registered.compute for name, registered in REGISTRY.items()} == before


class TestExtras:
    """Instrumentation read from where it actually lives."""

    def test_cache_counters_are_read_through_the_context(
        self, overlaps: list[TranscriptVariationAllele]
    ) -> None:
        """The zone map and codon cache count on themselves, not on the annotator."""
        annotator = build_annotator(annotator="consequence.staged", zone_map=True, codon_cache=True)
        for overlap in overlaps:
            annotator.annotate_feature_overlap(overlap)
        extras = annotator_extras(annotator)
        assert int(extras["zone map maps built"]) > 0  # type: ignore[call-overload]
        assert "zone map hits" in extras
        assert "zone map evictions" in extras
        assert "codon cache hits" in extras
        assert "settled early" in extras

    def test_an_uninstrumented_object_has_no_extras(self) -> None:
        """Nothing to report is an empty mapping, not an error."""
        assert not annotator_extras(object())
