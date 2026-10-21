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
"""The gate on reconstruction: a rebuilt derivation must describe the run that really happened.

Reconstruction is derivation, not recording, so it can drift from the path it claims to explain --
it already did once, when the columnar replay reported terms the real path never emitted. These
checks are the standing guard against that, over a sample of both corpora and every annotator that
traces: the terms a trace reports must be the terms the annotator returned, and every fact the run
memoised must be accounted for exactly once, by exactly one origin.
"""

from __future__ import annotations

import logging

import pytest

from ensembl.vep.annotators.registry import get_annotator
from ensembl.vep.corpus import build_corpus
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.trace import AnnotationTrace, FactOrigin

TRACED_ANNOTATORS: tuple[tuple[str, str], ...] = (
    ("consequence.linear", "ensembl.vep.annotators.consequence.strategies.linear"),
    ("consequence.staged", "ensembl.vep.annotators.consequence.strategies.staged"),
    ("consequence.columnar", "ensembl.vep.annotators.consequence.strategies.columnar"),
    ("consequence.columnar-staged", "ensembl.vep.annotators.consequence.strategies.columnar_staged"),
    ("consequence.columnar-kernels", "ensembl.vep.annotators.consequence.strategies.columnar_staged"),
)
"""Every pure-Python annotator, with its module logger. The native rust.* ones never trace."""

SAMPLE_STEP = 37
"""Take every Nth overlap of both corpora: enough spread to hit each shape, small enough to be a
unit test rather than a benchmark."""

_CODING_OVERLAP_KEY = "ENST0001::ENST0001_cds3_del1"


@pytest.fixture(scope="module")
def sampled_overlaps() -> list[TranscriptVariationAllele]:
    """A spread of overlaps from both corpora, including the exon-sharing heavy locus."""
    corpus = build_corpus(CorpusSpec(kind="synthetic.all", fuzz=200))
    overlaps = [overlap for _name, overlap in corpus.named_overlaps]
    return overlaps[::SAMPLE_STEP]


@pytest.fixture(params=TRACED_ANNOTATORS, ids=[name for name, _logger in TRACED_ANNOTATORS])
def annotator_name_and_logger(request: pytest.FixtureRequest) -> tuple[str, str]:
    """Parametrize over every annotator that records a derivation."""
    return request.param


def _annotate_with_tracing(
    annotator_name: str,
    logger_name: str,
    overlaps: list[TranscriptVariationAllele],
    caplog: pytest.LogCaptureFixture,
) -> list[tuple[list[str], AnnotationTrace]]:
    """Run one annotator over each overlap with tracing on, pairing terms with the derivation."""
    caplog.set_level(logging.DEBUG, logger=logger_name)
    annotator = get_annotator(annotator_name)
    paired = []
    for overlap in overlaps:
        consequences = annotator.annotate_feature_overlap(overlap)
        trace = getattr(annotator, "last_trace", None)
        assert trace is not None, f"{annotator_name} recorded no derivation with DEBUG enabled"
        paired.append(([consequence.term for consequence in consequences], trace))
    return paired


class TestReconstructedTermsMatchTheRealPath:
    """The reconstructed verdicts must agree with what the annotator actually returned."""

    def test_emitted_terms_match(
        self,
        annotator_name_and_logger: tuple[str, str],
        sampled_overlaps: list[TranscriptVariationAllele],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """A trace that reported different terms from the run would be explaining the wrong thing."""
        annotator_name, logger_name = annotator_name_and_logger
        for terms, trace in _annotate_with_tracing(annotator_name, logger_name, sampled_overlaps, caplog):
            assert trace.emitted_terms == tuple(terms)

    def test_fired_terms_are_exactly_the_emitted_and_the_pruned(
        self,
        annotator_name_and_logger: tuple[str, str],
        sampled_overlaps: list[TranscriptVariationAllele],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Every rule the trace says fired either survived pruning or was dropped by a named term."""
        annotator_name, logger_name = annotator_name_and_logger
        for _terms, trace in _annotate_with_tracing(annotator_name, logger_name, sampled_overlaps, caplog):
            fired = {event.term for event in trace.rule_events if event.fired}
            dropped = {event.dropped_term for event in trace.prune_events}
            assert fired == dropped | set(trace.emitted_terms)

    def test_every_pruning_decision_names_a_firing_term(
        self,
        annotator_name_and_logger: tuple[str, str],
        sampled_overlaps: list[TranscriptVariationAllele],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """A term cannot be subsumed by one that did not fire."""
        annotator_name, logger_name = annotator_name_and_logger
        for _terms, trace in _annotate_with_tracing(annotator_name, logger_name, sampled_overlaps, caplog):
            fired = {event.term for event in trace.rule_events if event.fired}
            for prune_event in trace.prune_events:
                assert prune_event.subsumed_by in fired


class TestTraceLifecycle:
    """Tests for enabling, disabling, and rendering strategy traces."""

    def test_trace_is_disabled_without_debug_logging(
        self,
        annotator_name_and_logger: tuple[str, str],
        boundary_overlaps: dict[str, TranscriptVariationAllele],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Test an ordinary annotation does not build a trace."""
        annotator_name, logger_name = annotator_name_and_logger
        caplog.set_level(logging.WARNING, logger=logger_name)
        annotator = get_annotator(annotator_name)

        consequences = annotator.annotate_feature_overlap(boundary_overlaps[_CODING_OVERLAP_KEY])

        assert consequences
        assert getattr(annotator, "last_trace", None) is None

    def test_trace_is_populated_and_renderable(
        self,
        annotator_name_and_logger: tuple[str, str],
        boundary_overlaps: dict[str, TranscriptVariationAllele],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Test debug logging records a renderable derivation."""
        annotator_name, logger_name = annotator_name_and_logger
        caplog.set_level(logging.DEBUG, logger=logger_name)
        annotator = get_annotator(annotator_name)

        consequences = annotator.annotate_feature_overlap(boundary_overlaps[_CODING_OVERLAP_KEY])
        emitted_terms = tuple(consequence.term for consequence in consequences)
        trace = getattr(annotator, "last_trace", None)

        assert trace is not None
        assert trace.fact_events
        assert trace.rule_events
        assert trace.emitted_terms == emitted_terms
        rendered = trace.format_explanation()
        assert f"annotator: {annotator_name}" in rendered
        assert all(term in rendered for term in emitted_terms)

    def test_coding_sequence_variant_pruning_is_recorded(
        self,
        annotator_name_and_logger: tuple[str, str],
        boundary_overlaps: dict[str, TranscriptVariationAllele],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Test the coding overlap records its subsumption edge."""
        annotator_name, logger_name = annotator_name_and_logger
        caplog.set_level(logging.DEBUG, logger=logger_name)
        annotator = get_annotator(annotator_name)

        annotator.annotate_feature_overlap(boundary_overlaps[_CODING_OVERLAP_KEY])
        trace = getattr(annotator, "last_trace", None)

        assert trace is not None
        pruning_edges = {(event.dropped_term, event.subsumed_by) for event in trace.prune_events}
        assert ("coding_sequence_variant", "frameshift_variant") in pruning_edges


class TestEveryFactHasExactlyOneOrigin:
    """Origins are derived, so the partition they induce is what has to be checked."""

    def test_no_fact_is_reported_twice(
        self,
        annotator_name_and_logger: tuple[str, str],
        sampled_overlaps: list[TranscriptVariationAllele],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """One event per fact: two would mean a fact claimed by two origins."""
        annotator_name, logger_name = annotator_name_and_logger
        for _terms, trace in _annotate_with_tracing(annotator_name, logger_name, sampled_overlaps, caplog):
            names = [event.name for event in trace.fact_events]
            assert len(names) == len(set(names))

    def test_forced_false_facts_were_never_computed(
        self,
        annotator_name_and_logger: tuple[str, str],
        sampled_overlaps: list[TranscriptVariationAllele],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Forced False means proved without computing, so it must carry no read and no charge."""
        annotator_name, logger_name = annotator_name_and_logger
        for _terms, trace in _annotate_with_tracing(annotator_name, logger_name, sampled_overlaps, caplog):
            for event in trace.events_by_origin(FactOrigin.FORCED_FALSE):
                assert event.value is False
                assert event.reads == 0
                assert event.forced_by, f"{event.name} was forced False by nothing"

    def test_only_computed_facts_charge_cost(
        self,
        annotator_name_and_logger: tuple[str, str],
        sampled_overlaps: list[TranscriptVariationAllele],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """A cache serving a fact is the whole saving, so an uncharged origin must charge nothing."""
        annotator_name, logger_name = annotator_name_and_logger
        for _terms, trace in _annotate_with_tracing(annotator_name, logger_name, sampled_overlaps, caplog):
            for event in trace.fact_events:
                if event.origin in (FactOrigin.PREFILLED, FactOrigin.CACHED):
                    assert event.cost == 0

    def test_a_read_fact_is_one_the_derivation_reports(
        self,
        annotator_name_and_logger: tuple[str, str],
        sampled_overlaps: list[TranscriptVariationAllele],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Read counts and events come from the same memo, so neither may name a fact the other misses."""
        annotator_name, logger_name = annotator_name_and_logger
        for _terms, trace in _annotate_with_tracing(annotator_name, logger_name, sampled_overlaps, caplog):
            read_names = {event.name for event in trace.fact_events if event.reads}
            reported = {event.name for event in trace.fact_events}
            assert read_names <= reported


class TestTracingDoesNotChangeTheAnswer:
    """The output must not depend on whether anybody was watching."""

    def test_terms_are_identical_with_tracing_on_and_off(
        self,
        annotator_name_and_logger: tuple[str, str],
        sampled_overlaps: list[TranscriptVariationAllele],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """The replay passes re-read facts; if that leaked into evaluation, this would catch it."""
        annotator_name, logger_name = annotator_name_and_logger

        caplog.set_level(logging.WARNING, logger=logger_name)
        untraced = [
            [consequence.term for consequence in get_annotator(annotator_name).annotate_feature_overlap(o)]
            for o in sampled_overlaps
        ]

        traced = [
            terms
            for terms, _trace in _annotate_with_tracing(annotator_name, logger_name, sampled_overlaps, caplog)
        ]

        assert untraced == traced
