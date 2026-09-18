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
"""Unit tests for transcript mapping."""

from __future__ import annotations

import pytest

from ensembl.vep.common import CoordSystemLevel, Strand, invalidate, reverse_complement
from ensembl.vep.common.codons import translate_dna
from ensembl.vep.models.genic import Exon, Transcript, Translation
from ensembl.vep.models.sequence import CoordSystem, Slice

CODING_SEQUENCE = "ATG" + "GGG" + "CCC" + "TAA"  # 12 bases, transcript-sense (M G P *)


class DummySlice(Slice):
    """A Slice whose cached sequence can be set/replaced through the unsafe test hook."""

    def __init__(self, length: int, *, strand: Strand = Strand.FORWARD) -> None:
        """Init."""
        super().__init__(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            start=1,
            end=length,
            strand=strand,
        )
        self.set_sequence("N" * length)

    def set_sequence(self, sequence: str) -> None:
        """Replace the cached sequence and invalidate downstream caches would need to see it."""
        with self.unsafe() as unsafe_slice:
            setattr(unsafe_slice, "_seq", sequence)


@pytest.fixture
def forward_transcript() -> Transcript:
    """A forward-strand, two-exon transcript with 3 bases of 5' UTR and CDS split across exons."""
    genomic_slice = DummySlice(40)
    exon_one = Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="E1", slice=genomic_slice)
    exon_two = Exon(start=20, end=30, strand=Strand.FORWARD, stable_id="E2", slice=genomic_slice)
    translation = Translation(start=4, end=24, stable_id="P1", exon_coding_start=4, exon_coding_end=5)
    transcript = Transcript(start=1, end=30, strand=Strand.FORWARD, stable_id="T_FWD", slice=genomic_slice)
    transcript.spliced_exons = [exon_one, exon_two]
    translation.transcript = transcript
    transcript.seed_attrs(translation=translation)

    genomic_chars = list(genomic_slice.seq or "")
    for cds_offset, base in enumerate(CODING_SEQUENCE):
        genomic_position = transcript.cds_offset_to_genomic(cds_offset)
        assert genomic_position is not None
        genomic_chars[genomic_position - 1] = base  # strand is FORWARD: no complement needed
    genomic_slice.set_sequence("".join(genomic_chars))
    invalidate(transcript, "_coding_sequence", "_translated_seq")
    return transcript


@pytest.fixture
def reverse_transcript() -> Transcript:
    """A reverse-strand transcript with the same exon geometry, mirrored, and the same CDS."""
    genomic_slice = DummySlice(40)
    exon_one = Exon(start=11, end=20, strand=Strand.REVERSE, stable_id="E1", slice=genomic_slice)
    exon_two = Exon(start=30, end=40, strand=Strand.REVERSE, stable_id="E2", slice=genomic_slice)
    translation = Translation(start=5, end=25, stable_id="P1", exon_coding_start=5, exon_coding_end=5)
    transcript = Transcript(start=11, end=40, strand=Strand.REVERSE, stable_id="T_REV", slice=genomic_slice)
    # Transcript order: the exon at the *higher* genomic coordinate is 5'-most on the reverse strand.
    transcript.spliced_exons = [exon_two, exon_one]
    translation.transcript = transcript
    transcript.seed_attrs(translation=translation)

    genomic_chars = list(genomic_slice.seq or "")
    for cds_offset, base in enumerate(CODING_SEQUENCE):
        genomic_position = transcript.cds_offset_to_genomic(cds_offset)
        assert genomic_position is not None
        # The stored genomic base is on the forward strand, so store the complement of the
        # transcript-sense base (single-base reverse_complement is just its complement).
        genomic_chars[genomic_position - 1] = reverse_complement(base)
    genomic_slice.set_sequence("".join(genomic_chars))
    invalidate(transcript, "_coding_sequence", "_translated_seq")
    return transcript


class TestForwardCDSMapping:
    """Test cds mapping for forward strand transcript."""

    def test_cds_length(self, forward_transcript: Transcript) -> None:
        """Cds length."""
        assert forward_transcript.cds_length == len(CODING_SEQUENCE)

    def test_cds_genomic_spans_ascending(self, forward_transcript: Transcript) -> None:
        """Cds genomic spans ascending."""
        spans = forward_transcript.cds_genomic_spans
        assert spans is not None
        assert spans == sorted(spans)
        assert spans == [(4, 10), (20, 24)]

    def test_coding_sequence_is_transcript_sense(self, forward_transcript: Transcript) -> None:
        """Coding sequence is transcript sense."""
        assert forward_transcript.coding_sequence == CODING_SEQUENCE

    def test_translated_seq_uses_translate_dna_not_placeholder(self, forward_transcript: Transcript) -> None:
        """Translated seq uses translate dna not placeholder."""
        assert forward_transcript.translated_seq == translate_dna(CODING_SEQUENCE)
        assert forward_transcript.translated_seq == "MGP*"
        assert forward_transcript.translated_seq != "dummy"

    @pytest.mark.parametrize("cds_offset", range(len(CODING_SEQUENCE)))
    def test_offset_genomic_round_trip(self, forward_transcript: Transcript, cds_offset: int) -> None:
        """Offset genomic round trip."""
        genomic_position = forward_transcript.cds_offset_to_genomic(cds_offset)
        assert genomic_position is not None
        assert forward_transcript.genomic_to_cds_offset(genomic_position) == cds_offset

    def test_start_and_stop_codon_genomic_spans(self, forward_transcript: Transcript) -> None:
        """Start and stop codon genomic spans."""
        assert forward_transcript.start_codon_genomic_spans() == [(4, 6)]
        assert forward_transcript.stop_codon_genomic_spans() == [(22, 24)]

    def test_genomic_to_cds_offset_none_outside_the_cds(self, forward_transcript: Transcript) -> None:
        """A position outside every CDS block maps to None (e.g. deep intron, UTR, downstream)."""
        assert forward_transcript.genomic_to_cds_offset(1) is None  # 5' UTR
        assert forward_transcript.genomic_to_cds_offset(15) is None  # intron
        assert forward_transcript.genomic_to_cds_offset(1000) is None  # far outside the transcript

    def test_cds_offset_to_genomic_none_out_of_range(self, forward_transcript: Transcript) -> None:
        """Cds offset to genomic none out of range."""
        assert forward_transcript.cds_offset_to_genomic(-1) is None
        cds_length = forward_transcript.cds_length
        assert cds_length is not None
        assert forward_transcript.cds_offset_to_genomic(cds_length) is None


class TestReverseCDSMapping:
    """Test cds mapping for reverse strand transcript."""

    def test_cds_length_matches_forward(
        self, forward_transcript: Transcript, reverse_transcript: Transcript
    ) -> None:
        """Cds length matches forward."""
        assert reverse_transcript.cds_length == forward_transcript.cds_length == len(CODING_SEQUENCE)

    def test_cds_genomic_spans_transcript_order(self, reverse_transcript: Transcript) -> None:
        """Spans are transcript order (block index 0 is the transcript's 5' exon), which on the
        reverse strand is genomic-*descending*, not ascending -- the mirror image of
        test_cds_genomic_spans_ascending's forward-strand check, which passes precisely because
        transcript order and genomic-ascending order coincide on the forward strand."""
        spans = reverse_transcript.cds_genomic_spans
        assert spans is not None
        assert spans == [(30, 36), (16, 20)]
        assert spans != sorted(spans)

    def test_coding_sequence_is_transcript_sense_not_genomic_forward(
        self, reverse_transcript: Transcript
    ) -> None:
        """The whole point of R1: this must equal CODING_SEQUENCE, not its reverse complement."""
        assert reverse_transcript.coding_sequence == CODING_SEQUENCE
        assert reverse_transcript.coding_sequence != reverse_complement(CODING_SEQUENCE)

    def test_translated_seq_matches_forward_strand_translation(
        self, forward_transcript: Transcript, reverse_transcript: Transcript
    ) -> None:
        """Translated seq matches forward strand translation."""
        assert reverse_transcript.translated_seq == forward_transcript.translated_seq == "MGP*"

    @pytest.mark.parametrize("cds_offset", range(len(CODING_SEQUENCE)))
    def test_offset_genomic_round_trip(self, reverse_transcript: Transcript, cds_offset: int) -> None:
        """Offset genomic round trip."""
        genomic_position = reverse_transcript.cds_offset_to_genomic(cds_offset)
        assert genomic_position is not None
        assert reverse_transcript.genomic_to_cds_offset(genomic_position) == cds_offset

    def test_offset_walks_genomic_descending_across_the_whole_cds(
        self, reverse_transcript: Transcript
    ) -> None:
        """Transcript order (ascending CDS offset) is genomic-descending on the reverse strand."""
        genomic_positions: list[int] = []
        for offset in range(len(CODING_SEQUENCE)):
            genomic_position = reverse_transcript.cds_offset_to_genomic(offset)
            assert genomic_position is not None
            genomic_positions.append(genomic_position)
        assert genomic_positions == sorted(genomic_positions, reverse=True)

    def test_start_codon_is_at_the_high_genomic_coordinate(self, reverse_transcript: Transcript) -> None:
        """5' (start codon) is the high genomic coordinate on the reverse strand."""
        start_spans = reverse_transcript.start_codon_genomic_spans()
        stop_spans = reverse_transcript.stop_codon_genomic_spans()
        assert start_spans is not None and stop_spans is not None
        # Spans are ascending, so span[0][0]/span[-1][1] are the true overall min/max --
        # min()/max() directly over a list of tuples would compare tuples lexicographically.
        assert start_spans[0][0] > stop_spans[-1][1]

    def test_genomic_to_cds_offset_none_outside_the_cds(self, reverse_transcript: Transcript) -> None:
        """A position outside every CDS block maps to None (e.g. deep intron, UTR, downstream)."""
        assert reverse_transcript.genomic_to_cds_offset(11) is None  # 3' UTR
        assert reverse_transcript.genomic_to_cds_offset(25) is None  # intron
        assert reverse_transcript.genomic_to_cds_offset(40) is None  # far outside the transcript

    def test_cds_offset_to_genomic_none_out_of_range(self, reverse_transcript: Transcript) -> None:
        """Cds offset to genomic none out of range."""
        assert reverse_transcript.cds_offset_to_genomic(-1) is None
        cds_length = reverse_transcript.cds_length
        assert cds_length is not None
        assert reverse_transcript.cds_offset_to_genomic(cds_length) is None


class TestNonCodingTranscriptMapping:
    """Test non-coding transcript mapping."""

    def test_non_coding_transcript_mapping_is_all_none(self) -> None:
        """A transcript with no translation has no CDS mapping at all."""
        genomic_slice = DummySlice(20)
        exon = Exon(start=1, end=20, strand=Strand.FORWARD, stable_id="E1", slice=genomic_slice)
        transcript = Transcript(start=1, end=20, strand=Strand.FORWARD, stable_id="T_NC", slice=genomic_slice)
        transcript.spliced_exons = [exon]

        assert transcript.cds_genomic_spans is None
        assert transcript.cds_length is None
        assert transcript.coding_sequence is None
        assert transcript.genomic_to_cds_offset(5) is None
        assert transcript.cds_offset_to_genomic(0) is None
        assert transcript.start_codon_genomic_spans() is None
        assert transcript.stop_codon_genomic_spans() is None


class TestCDNAMapper:
    """Test cdna mapping, which covers the whole spliced transcript, not just the CDS."""

    def test_genomic_spans_cover_every_exon(
        self, forward_transcript: Transcript, reverse_transcript: Transcript
    ) -> None:
        """genomic_spans is built from every spliced exon, transcript order, unlike mapper's
        CDS-clipped blocks."""
        for transcript in (forward_transcript, reverse_transcript):
            expected = [(exon.start, exon.end) for exon in transcript.spliced_exons]
            assert transcript.cdna_mapper is not None
            assert transcript.cdna_mapper.genomic_spans == expected

    def test_offsets_agree_with_cds_mapper_up_to_cdna_coding_start(
        self, forward_transcript: Transcript, reverse_transcript: Transcript
    ) -> None:
        """Every genomic position the CDS mapper resolves, the transcript-wide mapper resolves
        too, offset by exactly (cdna_coding_start - 1). The two mappers must never disagree
        about where a CDS base sits once the UTR-length offset is accounted for."""
        for transcript in (forward_transcript, reverse_transcript):
            mapper = transcript.cds_mapper
            cdna_mapper = transcript.cdna_mapper
            assert mapper is not None and cdna_mapper is not None
            assert transcript.cdna_coding_start is not None
            cdna_coding_start0 = transcript.cdna_coding_start - 1
            for cds_offset in range(mapper.length):
                genomic_position = mapper.offset_to_genomic(cds_offset)
                assert genomic_position is not None
                assert cdna_mapper.genomic_to_offset(genomic_position) == cds_offset + cdna_coding_start0

    def test_resolves_utr_positions_the_cds_mapper_cannot(self, forward_transcript: Transcript) -> None:
        """A 5' UTR genomic position maps to None on the CDS-only mapper but resolves on the
        transcript-wide one -- the whole reason _full_cds_alleles needs a second mapper."""
        assert forward_transcript.cds_mapper is not None
        assert forward_transcript.cds_mapper.genomic_to_offset(1) is None

        assert forward_transcript.cdna_mapper is not None
        assert forward_transcript.cdna_mapper.genomic_to_offset(1) == 0

    def test_none_with_no_exons(self) -> None:
        """A transcript with no exons at all has no transcript-wide mapper either."""
        transcript = Transcript(start=1, end=10, strand=Strand.FORWARD, stable_id="T_EMPTY")
        assert transcript.cdna_mapper is None

    def test_reassigning_spliced_exons_invalidates_the_mapper(self, forward_transcript: Transcript) -> None:
        """Reassigning spliced_exons must not leave a stale cdna_mapper behind."""
        first = forward_transcript.cdna_mapper
        assert first is not None
        forward_transcript.spliced_exons = list(forward_transcript.spliced_exons)
        assert forward_transcript.cdna_mapper is not first
