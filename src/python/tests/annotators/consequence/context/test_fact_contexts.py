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
"""Tests for the ``FactContext`` class."""

from __future__ import annotations

from ensembl.vep.annotators.consequence.facts.base import get_fact
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.trace import FactOrigin, TraceRecorder
from ensembl.vep.trace import reconstruct
from ensembl.vep.common import Strand
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import VariationAllele


def _create_transcript_variation_allele(
    transcript: Transcript, *, position: int
) -> TranscriptVariationAllele:
    """Create a single-nucleotide transcript variation allele."""
    variation_allele = VariationAllele(
        start=position,
        end=position,
        strand=Strand.FORWARD,
        slice=transcript.slice,
        reference_allele="A",
        alternate_allele="T",
    )
    return TranscriptVariationAllele(feature=transcript, overlapping_feature=variation_allele)


class TestContextConstruction:
    """What a freshly built context exposes."""

    def test_overlap_halves_are_reached_directly(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
    ) -> None:
        """Both halves of the overlap are read straight off ``overlap``."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        overlap = _create_transcript_variation_allele(transcript, position=transcript.start)
        context = FactContext(overlap=overlap)

        assert context.overlap is overlap
        assert context.overlap.feature is transcript
        assert context.overlap.overlapping_feature is overlap.overlapping_feature

    def test_recorder_defaults_to_none(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Tracing is opt-in: a FactContext built without a recorder has none."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        context = FactContext(
            overlap=_create_transcript_variation_allele(transcript, position=transcript.start)
        )
        assert context.recorder is None


class TestMemoisation:
    """The memo table: computed at most once, private per context, cleared on reset."""

    def test_memoisation_computes_a_fact_at_most_once(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
    ) -> None:
        """A fact already in the memo is served from it rather than recomputed.

        Observed through the read count and the single memo entry: two lookups, one
        materialisation.
        """
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        overlap = _create_transcript_variation_allele(transcript, position=transcript.start)
        context = FactContext(overlap=overlap, recorder=TraceRecorder())

        first_value = get_fact(context, "overlap_feature")
        second_value = get_fact(context, "overlap_feature")

        assert first_value == second_value
        assert context.recorder is not None
        assert context.recorder.read_counts["overlap_feature"] == 2
        events = [event for event in reconstruct.memo_fact_events(context) if event.name == "overlap_feature"]
        assert len(events) == 1, "must not be computed a second time"
        assert events[0].origin is FactOrigin.COMPUTED

    def test_two_contexts_over_the_same_overlap_do_not_share_memo(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
    ) -> None:
        """Two contexts over the same overlap do not share memo."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        overlap = _create_transcript_variation_allele(transcript, position=transcript.start)

        first_context = FactContext(overlap=overlap)
        second_context = FactContext(overlap=overlap)

        get_fact(first_context, "overlap_feature")

        assert first_context.memo
        assert second_context.memo == {}


class TestReset:
    """Everything ``reset`` clears, so a context can be reused for another overlap.

    The following tests are checked by giving the context a *new* overlap and reading the value back,
    rather than by inspecting the protected attribute that holds it: the stale read before the
    reset is the failure this guards against, and it is only visible through the public surface.
    """

    def test_reset_clears_the_memo(self, synthetic_transcripts_by_stable_id: dict[str, Transcript]) -> None:
        """Reset clears the memo table, so the next read recomputes."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        overlap = _create_transcript_variation_allele(transcript, position=transcript.start)
        context = FactContext(overlap=overlap)

        get_fact(context, "overlap_feature")
        assert context.memo

        context.reset()

        assert context.memo == {}

    def test_reset_clears_the_memoised_variant_span(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """The span is memoised per overlap, so it must be dropped along with the memo."""
        first = synthetic_transcripts_by_stable_id["ENST0001"]
        second = synthetic_transcripts_by_stable_id["ENST0003"]
        context = FactContext(overlap=_create_transcript_variation_allele(first, position=first.start))
        assert context.variant_span == (first.start, first.start)

        context.overlap = _create_transcript_variation_allele(second, position=second.start)
        assert context.variant_span == (first.start, first.start), "stale until reset"

        context.reset()

        assert context.variant_span == (second.start, second.start)

    def test_reset_clears_the_memoised_biotype(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """The biotype is resolved once per overlap, so it must be dropped along with the memo."""
        first = synthetic_transcripts_by_stable_id["ENST0001"]
        second = synthetic_transcripts_by_stable_id["ENST0003"]
        assert first.biotype != second.biotype
        context = FactContext(overlap=_create_transcript_variation_allele(first, position=first.start))
        assert context.biotype == first.biotype

        context.overlap = _create_transcript_variation_allele(second, position=second.start)
        assert context.biotype == first.biotype, "stale until reset"

        context.reset()

        assert context.biotype == second.biotype


class TestPrefill:
    """Seeding the memo from a cache, and the origin each seeded fact is credited with."""

    def test_prefill_seeds_the_memo(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
    ) -> None:
        """Prefill puts a value straight into the memo, and get_fact then serves it without
        recomputing -- the mechanism the zone-map and codon caches inject through."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        overlap = _create_transcript_variation_allele(transcript, position=transcript.start)
        context = FactContext(overlap=overlap, recorder=TraceRecorder())

        context.prefill({"overlap_feature": True})

        assert context.memo["overlap_feature"] is True
        assert get_fact(context, "overlap_feature") is True
        events = list(reconstruct.memo_fact_events(context))
        assert len(events) == 1
        assert events[0].origin is FactOrigin.PREFILLED

    def test_prefill_without_a_recorder_does_not_raise(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
    ) -> None:
        """prefill is a no-op with respect to tracing when self.recorder is None."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        context = FactContext(
            overlap=_create_transcript_variation_allele(transcript, position=transcript.start)
        )
        context.prefill({"overlap_feature": True})
        assert context.memo["overlap_feature"] is True

    def test_prefill_records_prefilled_origin_by_default(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
    ) -> None:
        """prefill's default origin is FactOrigin.PREFILLED (the zone map's origin)."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        context = FactContext(
            overlap=_create_transcript_variation_allele(transcript, position=transcript.start),
            recorder=TraceRecorder(),
        )
        context.prefill({"overlap_exon": True, "overlap_cds": False})

        events_by_name = {event.name: event for event in reconstruct.memo_fact_events(context)}
        assert events_by_name["overlap_exon"].origin is FactOrigin.PREFILLED
        assert events_by_name["overlap_exon"].value is True
        assert events_by_name["overlap_exon"].cost == 0
        assert events_by_name["overlap_cds"].origin is FactOrigin.PREFILLED
        assert events_by_name["overlap_cds"].value is False

    def test_prefill_records_cached_origin_when_requested(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
    ) -> None:
        """prefill(..., origin=FactOrigin.CACHED) is how the shared-codon cache tags its
        injections."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        context = FactContext(
            overlap=_create_transcript_variation_allele(transcript, position=transcript.start),
            recorder=TraceRecorder(),
        )
        context.prefill({"_codon_alleles": ("ATG", "ACG")}, origin=FactOrigin.CACHED)

        event = reconstruct.memo_fact_events(context)[0]
        assert event.name == "_codon_alleles"
        assert event.origin is FactOrigin.CACHED
        assert event.value == ("ATG", "ACG")

    def test_prefill_does_not_claim_a_fact_the_memo_already_held(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
    ) -> None:
        """A codon-cache miss computes the facts and then prefills them back; they stay COMPUTED.

        This is what lets an origin be derived rather than recorded: only the names a prefill call
        actually introduced to the memo are the ones a cache genuinely served.
        """
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        context = FactContext(
            overlap=_create_transcript_variation_allele(transcript, position=transcript.start),
            recorder=TraceRecorder(),
        )
        get_fact(context, "overlap_feature")

        context.prefill({"overlap_feature": True}, origin=FactOrigin.CACHED)

        events = list(reconstruct.memo_fact_events(context))
        assert len(events) == 1
        assert events[0].origin is FactOrigin.COMPUTED
