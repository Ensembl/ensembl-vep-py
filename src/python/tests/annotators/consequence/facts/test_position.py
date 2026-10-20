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
"""Tests for consequence position facts."""

from __future__ import annotations

from collections.abc import Callable

import pytest

from ensembl.vep.annotators.consequence.facts.base import get_fact
from ensembl.vep.annotators.consequence.context.annotator import AnnotatorContext
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.common import Strand
from ensembl.vep.models.genic import Exon, Intron, Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import VariationAllele


def _create_transcript_variation_allele(
    transcript: Transcript, *, start: int, end: int, reference: str, alternate: str
) -> TranscriptVariationAllele:
    """Create a transcript variation allele at explicit genomic coordinates."""
    variation_allele = VariationAllele(
        start=start,
        end=end,
        strand=Strand.FORWARD,
        slice=transcript.slice,
        reference_allele=reference,
        alternate_allele=alternate,
    )
    return TranscriptVariationAllele(feature=transcript, overlapping_feature=variation_allele)


def _context_at(transcript: Transcript, start: int, end: int | None = None) -> FactContext:
    """A fact context for a substitution spanning ``[start, end]``, one base when ``end`` is None.

    Always a substitution of equal length, so that nothing under test can be swayed by the
    variant's shape; every fact here is a function of its coordinates alone.
    """
    span_end = start if end is None else end
    span_length = span_end - start + 1
    overlap = _create_transcript_variation_allele(
        transcript,
        start=start,
        end=span_end,
        reference="A" * span_length,
        alternate="T" * span_length,
    )
    return FactContext(overlap=overlap)


def _exons_in_genomic_order(transcript: Transcript) -> list[Exon]:
    """The transcript's exons in ascending genomic order.

    ``spliced_exons`` is in transcript order, so on the reverse strand index 0 is the
    highest-coordinate exon; sorting keeps the placement tables below strand-agnostic.
    """
    return sorted(transcript.spliced_exons, key=lambda exon: exon.start)


def _intron_offset_position(intron: Intron, strand: Strand, *, offset: int, from_donor: bool) -> int:
    """Genomic position of a 1-based, transcript-sense intron offset.

    Restated here rather than imported from ``facts.helper`` so that the expectations below do
    not inherit the implementation's own view of which intron end is the donor.

    Args:
        intron (Intron): The intron to measure from.
        strand (Strand): The transcript's strand.
        offset (int): 1-based offset into the intron, from the donor end when ``from_donor``
            is True and from the acceptor end otherwise.
        from_donor (bool): Measure from the donor (5') end rather than the acceptor (3') end.

    Returns:
        int: The genomic position of that offset.
    """
    if strand == Strand.FORWARD:
        donor_anchor, acceptor_anchor, step = intron.start, intron.end, 1
    else:
        donor_anchor, acceptor_anchor, step = intron.end, intron.start, -1
    if from_donor:
        return donor_anchor + step * (offset - 1)
    return acceptor_anchor - step * (offset - 1)


def _first_exon(transcript: Transcript) -> Exon:
    return _exons_in_genomic_order(transcript)[0]


def _second_exon(transcript: Transcript) -> Exon:
    return _exons_in_genomic_order(transcript)[1]


_SPLICE_FACTS = (
    "overlap_donor_site",
    "overlap_acceptor_site",
    "overlap_donor_5th_base",
    "overlap_donor_region",
    "overlap_polypyrimidine_tract",
    "overlap_splice_region",
)


class TestRelativeToFeatureFacts:
    """Position facts relative to the whole feature."""

    @pytest.mark.parametrize("stable_id", ["ENST0001", "ENST0002"])
    @pytest.mark.parametrize(
        ("span", "expected"),
        [
            pytest.param(lambda start, end: (start - 20, start - 11), False, id="entirely_before"),
            pytest.param(lambda start, end: (start - 10, start - 1), False, id="abuts_start"),
            pytest.param(lambda start, end: (start - 10, start + 10), True, id="straddles_start"),
            pytest.param(lambda start, end: (start, start), True, id="at_first_base"),
            pytest.param(lambda start, end: (start + 10, start + 20), True, id="entirely_within"),
            pytest.param(lambda start, end: (start, end), True, id="exactly_the_feature"),
            pytest.param(lambda start, end: (end, end), True, id="at_last_base"),
            pytest.param(lambda start, end: (end - 10, end + 10), True, id="straddles_end"),
            pytest.param(lambda start, end: (start - 10, end + 10), True, id="contains_the_feature"),
            pytest.param(lambda start, end: (end + 1, end + 10), False, id="abuts_end"),
            pytest.param(lambda start, end: (end + 11, end + 20), False, id="entirely_after"),
        ],
    )
    def test_overlap_feature(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        stable_id: str,
        span: Callable[[int, int], tuple[int, int]],
        expected: bool,
    ) -> None:
        """A variant overlaps the feature whenever the two spans share at least one base."""
        transcript = synthetic_transcripts_by_stable_id[stable_id]
        start, end = span(transcript.start, transcript.end)
        assert get_fact(_context_at(transcript, start, end), "overlap_feature") is expected

    @pytest.mark.parametrize("stable_id", ["ENST0001", "ENST0002"])
    @pytest.mark.parametrize(
        ("span", "expected"),
        [
            pytest.param(lambda start, end: (start - 20, start - 11), False, id="entirely_before"),
            pytest.param(lambda start, end: (start - 10, start - 1), False, id="abuts_start"),
            pytest.param(lambda start, end: (start - 10, start + 10), False, id="straddles_start"),
            pytest.param(lambda start, end: (start, start), False, id="at_first_base"),
            pytest.param(lambda start, end: (start + 10, start + 20), False, id="entirely_within"),
            pytest.param(lambda start, end: (start, end), True, id="exactly_the_feature"),
            pytest.param(lambda start, end: (end, end), False, id="at_last_base"),
            pytest.param(lambda start, end: (end - 10, end + 10), False, id="straddles_end"),
            pytest.param(lambda start, end: (start - 10, end + 10), True, id="contains_the_feature"),
            pytest.param(lambda start, end: (end + 1, end + 10), False, id="abuts_end"),
            pytest.param(lambda start, end: (end + 11, end + 20), False, id="entirely_after"),
        ],
    )
    def test_contain_feature(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        stable_id: str,
        span: Callable[[int, int], tuple[int, int]],
        expected: bool,
    ) -> None:
        """A variant contains the feature only when it covers every base of it."""
        transcript = synthetic_transcripts_by_stable_id[stable_id]
        start, end = span(transcript.start, transcript.end)
        assert get_fact(_context_at(transcript, start, end), "contain_feature") is expected

    @pytest.mark.parametrize("stable_id", ["ENST0001", "ENST0002"])
    @pytest.mark.parametrize(
        ("span", "expected"),
        [
            pytest.param(lambda start, end: (start - 20, start - 11), False, id="entirely_before"),
            pytest.param(lambda start, end: (start - 10, start - 1), False, id="abuts_start"),
            pytest.param(lambda start, end: (start - 10, start + 10), False, id="straddles_start"),
            pytest.param(lambda start, end: (start, start), True, id="at_first_base"),
            pytest.param(lambda start, end: (start + 10, start + 20), True, id="entirely_within"),
            pytest.param(lambda start, end: (start, end), True, id="exactly_the_feature"),
            pytest.param(lambda start, end: (end, end), True, id="at_last_base"),
            pytest.param(lambda start, end: (end - 10, end + 10), False, id="straddles_end"),
            pytest.param(lambda start, end: (start - 10, end + 10), False, id="contains_the_feature"),
            pytest.param(lambda start, end: (end + 1, end + 10), False, id="abuts_end"),
            pytest.param(lambda start, end: (end + 11, end + 20), False, id="entirely_after"),
        ],
    )
    def test_within_feature(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        stable_id: str,
        span: Callable[[int, int], tuple[int, int]],
        expected: bool,
    ) -> None:
        """A variant is within the feature only when every base of it is inside the feature."""
        transcript = synthetic_transcripts_by_stable_id[stable_id]
        start, end = span(transcript.start, transcript.end)
        assert get_fact(_context_at(transcript, start, end), "within_feature") is expected

    def test_contain_and_within_are_both_true_only_for_an_exact_span_match(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """The two facts meet only where the variant and the feature have identical spans."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        exact = _context_at(transcript, transcript.start, transcript.end)
        larger = _context_at(transcript, transcript.start - 1, transcript.end + 1)
        smaller = _context_at(transcript, transcript.start + 1, transcript.end - 1)

        assert (get_fact(exact, "contain_feature"), get_fact(exact, "within_feature")) == (True, True)
        assert (get_fact(larger, "contain_feature"), get_fact(larger, "within_feature")) == (True, False)
        assert (get_fact(smaller, "contain_feature"), get_fact(smaller, "within_feature")) == (False, True)

    @pytest.mark.parametrize("stable_id", ["ENST0001", "ENST0002"])
    @pytest.mark.parametrize(
        ("span", "forward_upstream", "forward_downstream"),
        [
            pytest.param(lambda start, end: (start - 1, start - 1), True, False, id="one_base_before_start"),
            pytest.param(
                lambda start, end: (start - 5000, start - 5000), True, False, id="5000_before_start"
            ),
            pytest.param(
                lambda start, end: (start - 5001, start - 5001), False, False, id="5001_before_start"
            ),
            pytest.param(
                lambda start, end: (start - 5001, start - 5000), True, False, id="spans_the_5000_boundary"
            ),
            pytest.param(lambda start, end: (end + 1, end + 1), False, True, id="one_base_after_end"),
            pytest.param(lambda start, end: (end + 5000, end + 5000), False, True, id="5000_after_end"),
            pytest.param(lambda start, end: (end + 5001, end + 5001), False, False, id="5001_after_end"),
            pytest.param(lambda start, end: (start + 10, start + 20), False, False, id="inside_the_feature"),
            pytest.param(lambda start, end: (start - 10, start + 10), False, False, id="straddles_the_start"),
        ],
    )
    def test_upstream_and_downstream(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        stable_id: str,
        span: Callable[[int, int], tuple[int, int]],
        forward_upstream: bool,
        forward_downstream: bool,
    ) -> None:
        """Flanking placement, with the 5kb window inclusive at both ends of the feature.

        The table states the forward-strand answer; on the reverse strand the very same genomic
        placement swaps roles, which is the whole of what strand means to these two facts.
        """
        transcript = synthetic_transcripts_by_stable_id[stable_id]
        start, end = span(transcript.start, transcript.end)
        context = _context_at(transcript, start, end)

        if transcript.strand == Strand.FORWARD:
            expected_upstream, expected_downstream = forward_upstream, forward_downstream
        else:
            expected_upstream, expected_downstream = forward_downstream, forward_upstream

        assert get_fact(context, "is_upstream") is expected_upstream
        assert get_fact(context, "is_downstream") is expected_downstream

    def test_upstream_and_downstream_are_never_both_true(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """No single placement is on both sides of a feature at once."""
        placements = [
            lambda start, end: (start - 1, start - 1),
            lambda start, end: (start - 5000, start - 5000),
            lambda start, end: (start - 5001, start - 5001),
            lambda start, end: (start - 5001, start - 5000),
            lambda start, end: (end + 1, end + 1),
            lambda start, end: (end + 5000, end + 5000),
            lambda start, end: (end + 5001, end + 5001),
            lambda start, end: (start + 10, start + 20),
            lambda start, end: (start - 10, start + 10),
        ]
        for stable_id in ["ENST0001", "ENST0002"]:
            transcript = synthetic_transcripts_by_stable_id[stable_id]
            for span in placements:
                start, end = span(transcript.start, transcript.end)
                context = _context_at(transcript, start, end)
                assert not (get_fact(context, "is_upstream") and get_fact(context, "is_downstream"))


class TestConfiguredDistances:
    """``--distance`` reaches the two flanking facts through the annotator context."""

    @staticmethod
    def _context(transcript: Transcript, start: int, *, upstream: int, downstream: int) -> FactContext:
        """A one-base context annotated under the given pairing distances."""
        context = _context_at(transcript, start)
        return FactContext(
            overlap=context.overlap,
            annotator_context=AnnotatorContext(upstream_distance=upstream, downstream_distance=downstream),
        )

    def test_a_narrower_upstream_distance_excludes_a_placement_the_default_includes(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """What counts as upstream is the run's distance, not the packaged default."""
        transcript = next(
            candidate
            for candidate in synthetic_transcripts_by_stable_id.values()
            if candidate.strand == Strand.FORWARD
        )
        position = transcript.start - 2000

        assert get_fact(_context_at(transcript, position), "is_upstream") is True
        narrowed = self._context(transcript, position, upstream=1000, downstream=5000)
        assert get_fact(narrowed, "is_upstream") is False

    def test_the_two_sides_are_independent(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """A narrow upstream distance leaves the downstream side at its own value."""
        transcript = next(
            candidate
            for candidate in synthetic_transcripts_by_stable_id.values()
            if candidate.strand == Strand.FORWARD
        )

        upstream = self._context(transcript, transcript.start - 2000, upstream=1000, downstream=5000)
        downstream = self._context(transcript, transcript.end + 2000, upstream=1000, downstream=5000)

        assert get_fact(upstream, "is_upstream") is False
        assert get_fact(downstream, "is_downstream") is True

    def test_the_distances_follow_the_strand(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """On the reverse strand the upstream distance applies past the transcript's end."""
        transcript = next(
            candidate
            for candidate in synthetic_transcripts_by_stable_id.values()
            if candidate.strand == Strand.REVERSE
        )

        above = self._context(transcript, transcript.end + 2000, upstream=5000, downstream=1000)
        below = self._context(transcript, transcript.start - 2000, upstream=5000, downstream=1000)

        assert get_fact(above, "is_upstream") is True
        assert get_fact(below, "is_downstream") is False


class TestWithinFeatureFacts:
    """Facts that depend on variant position within a transcript's exons and introns."""

    @pytest.mark.parametrize("stable_id", ["ENST0001", "ENST0002"])
    @pytest.mark.parametrize(
        ("span", "expected_exon", "expected_within_exon", "expected_intron"),
        [
            pytest.param(
                lambda transcript: (_first_exon(transcript).start + 10, _first_exon(transcript).start + 20),
                True,
                True,
                False,
                id="inside_an_exon",
            ),
            pytest.param(
                lambda transcript: (_first_exon(transcript).start, _first_exon(transcript).end),
                True,
                True,
                False,
                id="the_whole_first_exon",
            ),
            pytest.param(
                lambda transcript: (_first_exon(transcript).end, _first_exon(transcript).end),
                True,
                True,
                False,
                id="last_base_of_an_exon",
            ),
            pytest.param(
                lambda transcript: (transcript.introns[0].start + 10, transcript.introns[0].start + 20),
                False,
                False,
                True,
                id="inside_an_intron",
            ),
            pytest.param(
                lambda transcript: (transcript.introns[0].start, transcript.introns[0].start),
                False,
                False,
                True,
                id="first_base_of_an_intron",
            ),
            pytest.param(
                lambda transcript: (_first_exon(transcript).end - 2, transcript.introns[0].start + 2),
                True,
                False,
                True,
                id="straddles_the_exon_intron_boundary",
            ),
            pytest.param(
                lambda transcript: (_first_exon(transcript).start, _second_exon(transcript).end),
                True,
                False,
                True,
                id="spans_two_exons_and_the_intron_between",
            ),
            pytest.param(
                lambda transcript: (transcript.start - 100, transcript.start - 50),
                False,
                False,
                False,
                id="outside_the_transcript",
            ),
        ],
    )
    def test_exon_and_intron_placement(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        stable_id: str,
        span: Callable[[Transcript], tuple[int, int]],
        expected_exon: bool,
        expected_within_exon: bool,
        expected_intron: bool,
    ) -> None:
        """Exon and intron placement, including the spans that straddle a boundary."""
        transcript = synthetic_transcripts_by_stable_id[stable_id]
        start, end = span(transcript)
        context = _context_at(transcript, start, end)

        assert get_fact(context, "overlap_exon") is expected_exon
        assert get_fact(context, "within_any_exon") is expected_within_exon
        assert get_fact(context, "overlap_intron") is expected_intron

    @pytest.mark.parametrize("stable_id", ["ENST0001", "ENST0002"])
    @pytest.mark.parametrize(
        ("offset", "expected_true"),
        [
            pytest.param(1, frozenset({"overlap_donor_site"}), id="donor+1"),
            pytest.param(2, frozenset({"overlap_donor_site"}), id="donor+2"),
            pytest.param(3, frozenset({"overlap_donor_region", "overlap_splice_region"}), id="donor+3"),
            pytest.param(4, frozenset({"overlap_donor_region", "overlap_splice_region"}), id="donor+4"),
            pytest.param(
                5,
                frozenset({"overlap_donor_5th_base", "overlap_donor_region", "overlap_splice_region"}),
                id="donor+5",
            ),
            pytest.param(6, frozenset({"overlap_donor_region", "overlap_splice_region"}), id="donor+6"),
            pytest.param(7, frozenset({"overlap_splice_region"}), id="donor+7"),
            pytest.param(8, frozenset({"overlap_splice_region"}), id="donor+8"),
            pytest.param(9, frozenset(), id="donor+9"),
        ],
    )
    def test_splice_facts_at_donor_offsets(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        stable_id: str,
        offset: int,
        expected_true: frozenset[str],
    ) -> None:
        """Each intron offset into the donor end fires exactly the splice facts that claim it."""
        transcript = synthetic_transcripts_by_stable_id[stable_id]
        intron = transcript.introns[0]
        position = _intron_offset_position(intron, transcript.strand, offset=offset, from_donor=True)
        context = _context_at(transcript, position)

        for fact_name in _SPLICE_FACTS:
            assert get_fact(context, fact_name) is (fact_name in expected_true), fact_name

    @pytest.mark.parametrize("stable_id", ["ENST0001", "ENST0002"])
    @pytest.mark.parametrize(
        ("offset", "expected_true"),
        [
            pytest.param(1, frozenset({"overlap_acceptor_site"}), id="acceptor-1"),
            pytest.param(2, frozenset({"overlap_acceptor_site"}), id="acceptor-2"),
            pytest.param(
                3,
                frozenset({"overlap_polypyrimidine_tract", "overlap_splice_region"}),
                id="acceptor-3",
            ),
            pytest.param(
                4,
                frozenset({"overlap_polypyrimidine_tract", "overlap_splice_region"}),
                id="acceptor-4",
            ),
            pytest.param(
                5,
                frozenset({"overlap_polypyrimidine_tract", "overlap_splice_region"}),
                id="acceptor-5",
            ),
            pytest.param(
                6,
                frozenset({"overlap_polypyrimidine_tract", "overlap_splice_region"}),
                id="acceptor-6",
            ),
            pytest.param(
                7,
                frozenset({"overlap_polypyrimidine_tract", "overlap_splice_region"}),
                id="acceptor-7",
            ),
            pytest.param(
                8,
                frozenset({"overlap_polypyrimidine_tract", "overlap_splice_region"}),
                id="acceptor-8",
            ),
            pytest.param(9, frozenset({"overlap_polypyrimidine_tract"}), id="acceptor-9"),
            pytest.param(10, frozenset({"overlap_polypyrimidine_tract"}), id="acceptor-10"),
            pytest.param(16, frozenset({"overlap_polypyrimidine_tract"}), id="acceptor-16"),
            pytest.param(17, frozenset({"overlap_polypyrimidine_tract"}), id="acceptor-17"),
            pytest.param(18, frozenset(), id="acceptor-18"),
        ],
    )
    def test_splice_facts_at_acceptor_offsets(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        stable_id: str,
        offset: int,
        expected_true: frozenset[str],
    ) -> None:
        """Each intron offset into the acceptor end fires exactly the splice facts that claim it."""
        transcript = synthetic_transcripts_by_stable_id[stable_id]
        intron = transcript.introns[0]
        position = _intron_offset_position(intron, transcript.strand, offset=offset, from_donor=False)
        context = _context_at(transcript, position)

        for fact_name in _SPLICE_FACTS:
            assert get_fact(context, fact_name) is (fact_name in expected_true), fact_name

    @pytest.mark.parametrize("stable_id", ["ENST0001", "ENST0002"])
    @pytest.mark.parametrize(
        ("exon_offset", "expected"),
        [
            pytest.param(1, True, id="1bp"),
            pytest.param(2, True, id="2bp"),
            pytest.param(3, True, id="3bp"),
            pytest.param(4, False, id="4bp"),
        ],
    )
    def test_splice_region_reaches_three_bases_into_the_flanking_exons(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        stable_id: str,
        exon_offset: int,
        expected: bool,
    ) -> None:
        """The exonic half of the splice region is 1-3 bases out from either intron edge.

        Measured in genomic coordinates from the intron's own edges, since the fact treats both
        edges alike -- which is itself the claim, on both strands.
        """
        transcript = synthetic_transcripts_by_stable_id[stable_id]
        intron = transcript.introns[0]

        before_intron = _context_at(transcript, intron.start - exon_offset)
        after_intron = _context_at(transcript, intron.end + exon_offset)

        assert get_fact(before_intron, "overlap_splice_region") is expected
        assert get_fact(after_intron, "overlap_splice_region") is expected

    def test_donor_is_the_low_coordinate_intron_end_on_the_forward_strand(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """On the forward strand the donor (+1) is the intron's low genomic coordinate."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        intron = transcript.introns[0]

        low_end = _context_at(transcript, intron.start)
        high_end = _context_at(transcript, intron.end)

        assert get_fact(low_end, "overlap_donor_site") is True
        assert get_fact(low_end, "overlap_acceptor_site") is False
        assert get_fact(high_end, "overlap_acceptor_site") is True
        assert get_fact(high_end, "overlap_donor_site") is False

    def test_donor_is_the_high_coordinate_intron_end_on_the_reverse_strand(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """On the reverse strand the donor (+1) is the intron's high genomic coordinate."""
        transcript = synthetic_transcripts_by_stable_id["ENST0002"]
        intron = transcript.introns[0]

        low_end = _context_at(transcript, intron.start)
        high_end = _context_at(transcript, intron.end)

        assert get_fact(high_end, "overlap_donor_site") is True
        assert get_fact(high_end, "overlap_acceptor_site") is False
        assert get_fact(low_end, "overlap_acceptor_site") is True
        assert get_fact(low_end, "overlap_donor_site") is False

    @pytest.mark.parametrize("stable_id", ["ENST0001", "ENST0002"])
    @pytest.mark.parametrize(
        ("position", "expected"),
        [
            pytest.param(lambda intron: intron.start, False, id="first_intron_base"),
            pytest.param(lambda intron: intron.start + 1, False, id="second_intron_base"),
            pytest.param(lambda intron: intron.start + 2, True, id="third_intron_base"),
            pytest.param(lambda intron: (intron.start + intron.end) // 2, True, id="middle"),
            pytest.param(lambda intron: intron.end - 2, True, id="third_last_intron_base"),
            pytest.param(lambda intron: intron.end - 1, False, id="second_last_intron_base"),
            pytest.param(lambda intron: intron.end, False, id="last_intron_base"),
        ],
    )
    def test_the_intron_interior_excludes_two_bases_at_each_end(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        stable_id: str,
        position: Callable[[Intron], int],
        expected: bool,
    ) -> None:
        """VEP's ``intronic`` window is ``intron.start + 2 .. intron.end - 2``, whatever the strand."""
        transcript = synthetic_transcripts_by_stable_id[stable_id]
        intron = transcript.introns[0]
        context = _context_at(transcript, position(intron))

        assert get_fact(context, "overlap_intron") is True
        assert get_fact(context, "overlap_intron_interior") is expected

    @pytest.mark.parametrize(
        ("insertion_start", "expected"),
        [
            pytest.param(lambda intron: intron.start + 1, False, id="after_the_second_base"),
            pytest.param(lambda intron: intron.start + 2, True, id="before_the_third_base"),
            pytest.param(lambda intron: intron.end - 1, True, id="after_the_third_last_base"),
            pytest.param(lambda intron: intron.end, False, id="after_the_second_last_base"),
        ],
    )
    def test_an_insertion_at_the_edge_of_the_interior_counts(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        insertion_start: Callable[[Intron], int],
        expected: bool,
    ) -> None:
        """VEP's insertion special case: an insertion flanked by an interior base is intronic."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        intron = transcript.introns[0]
        start = insertion_start(intron)
        overlap = _create_transcript_variation_allele(
            transcript, start=start, end=start - 1, reference="-", alternate="AT"
        )

        assert get_fact(FactContext(overlap=overlap), "overlap_intron_interior") is expected

    def test_splice_facts_are_false_deep_inside_an_intron(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Well away from either edge, only the intron facts remain true."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        intron = transcript.introns[0]
        context = _context_at(transcript, (intron.start + intron.end) // 2)

        assert get_fact(context, "overlap_intron") is True
        assert get_fact(context, "overlap_intron_boundary") is False
        for fact_name in _SPLICE_FACTS:
            assert get_fact(context, fact_name) is False, fact_name


class TestWithinCodingRegionFacts:
    """Facts that need genomic-to-CDS coordinate mapping."""

    @pytest.mark.parametrize("stable_id", ["ENST0001", "ENST0002"])
    def test_overlap_cds_and_utr_are_mutually_exclusive(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript], stable_id: str
    ) -> None:
        """CDS and UTR overlap are mutually exclusive."""
        transcript = synthetic_transcripts_by_stable_id[stable_id]
        five_prime_spans = transcript.five_prime_utr_genomic_spans
        three_prime_spans = transcript.three_prime_utr_genomic_spans
        assert five_prime_spans is not None and three_prime_spans is not None
        cds_position = transcript.cds_offset_to_genomic(30)
        assert cds_position is not None

        five_prime_context = _context_at(transcript, five_prime_spans[0][0])
        three_prime_context = _context_at(transcript, three_prime_spans[0][0])
        cds_context = _context_at(transcript, cds_position)

        assert get_fact(five_prime_context, "overlap_five_prime_utr") is True
        assert get_fact(five_prime_context, "overlap_three_prime_utr") is False
        assert get_fact(five_prime_context, "overlap_cds") is False

        assert get_fact(three_prime_context, "overlap_three_prime_utr") is True
        assert get_fact(three_prime_context, "overlap_five_prime_utr") is False
        assert get_fact(three_prime_context, "overlap_cds") is False

        assert get_fact(cds_context, "overlap_cds") is True
        assert get_fact(cds_context, "overlap_five_prime_utr") is False
        assert get_fact(cds_context, "overlap_three_prime_utr") is False

    def test_the_five_prime_utr_changes_genomic_end_with_the_strand(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """The 5' UTR is at the low genomic end on ``+`` and the high genomic end on ``-``."""
        forward = synthetic_transcripts_by_stable_id["ENST0001"]
        reverse = synthetic_transcripts_by_stable_id["ENST0002"]

        at_forward_start = _context_at(forward, forward.start)
        at_reverse_end = _context_at(reverse, reverse.end)

        assert get_fact(at_forward_start, "overlap_five_prime_utr") is True
        assert get_fact(at_forward_start, "overlap_three_prime_utr") is False
        assert get_fact(at_reverse_end, "overlap_five_prime_utr") is True
        assert get_fact(at_reverse_end, "overlap_three_prime_utr") is False

    def test_utr_facts_are_false_for_a_transcript_with_no_translation(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """A non-coding transcript has no UTR spans, so neither UTR fact can fire in its exons."""
        transcript = synthetic_transcripts_by_stable_id["ENST0003"]
        assert transcript.five_prime_utr_genomic_spans is None
        context = _context_at(transcript, _first_exon(transcript).start + 10)

        assert get_fact(context, "overlap_exon") is True
        assert get_fact(context, "overlap_five_prime_utr") is False
        assert get_fact(context, "overlap_three_prime_utr") is False
        assert get_fact(context, "overlap_cds") is False

    @pytest.mark.parametrize("stable_id", ["ENST0001", "ENST0002"])
    def test_overlaps_start_and_stop_codon(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript], stable_id: str
    ) -> None:
        """Start and stop codons are detected independently, on either strand."""
        transcript = synthetic_transcripts_by_stable_id[stable_id]
        start_spans = transcript.start_codon_genomic_spans()
        stop_spans = transcript.stop_codon_genomic_spans()
        assert start_spans is not None and stop_spans is not None

        start_context = _context_at(transcript, start_spans[0][0])
        stop_context = _context_at(transcript, stop_spans[0][0])

        assert get_fact(start_context, "overlap_start_codon") is True
        assert get_fact(start_context, "overlap_stop_codon") is False
        assert get_fact(stop_context, "overlap_stop_codon") is True
        assert get_fact(stop_context, "overlap_start_codon") is False

    @pytest.mark.parametrize("stable_id", ["ENST0001", "ENST0002"])
    @pytest.mark.parametrize("codon_base", [0, 1, 2])
    def test_every_start_codon_base_is_detected(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        stable_id: str,
        codon_base: int,
    ) -> None:
        """All three start-codon bases count, not just the first."""
        transcript = synthetic_transcripts_by_stable_id[stable_id]
        start_spans = transcript.start_codon_genomic_spans()
        assert start_spans is not None
        span_start, span_end = start_spans[0]
        assert span_end - span_start == 2

        context = _context_at(transcript, span_start + codon_base)
        assert get_fact(context, "overlap_start_codon") is True

    @pytest.mark.parametrize("stable_id", ["ENST0001", "ENST0002"])
    def test_the_base_beyond_the_start_codon_is_not_part_of_it(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript], stable_id: str
    ) -> None:
        """The base immediately outside the start codon is still CDS, but not the start codon."""
        transcript = synthetic_transcripts_by_stable_id[stable_id]
        start_spans = transcript.start_codon_genomic_spans()
        assert start_spans is not None
        span_start, span_end = start_spans[0]
        beyond = span_end + 1 if transcript.strand == Strand.FORWARD else span_start - 1

        context = _context_at(transcript, beyond)
        assert get_fact(context, "overlap_start_codon") is False
        assert get_fact(context, "overlap_cds") is True

    @pytest.mark.parametrize("stable_id", ["ENST0001", "ENST0002"])
    def test_partial_last_codon_is_false_across_a_complete_cds(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript], stable_id: str
    ) -> None:
        """No position in a CDS whose length is a multiple of three sits in a partial codon."""
        transcript = synthetic_transcripts_by_stable_id[stable_id]
        cds_length = transcript.cds_length
        assert cds_length is not None and cds_length % 3 == 0

        for cds_offset in (0, 1, cds_length // 2, cds_length - 2, cds_length - 1):
            genomic_position = transcript.cds_offset_to_genomic(cds_offset)
            assert genomic_position is not None
            context = _context_at(transcript, genomic_position)
            assert get_fact(context, "partial_last_codon") is False, cds_offset
