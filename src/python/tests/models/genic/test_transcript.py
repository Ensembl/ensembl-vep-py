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
"""Unit tests for Transcript class."""

import logging

import pytest

from ensembl.vep.common.enums import CoordSystemLevel, Strand, ExonPhase
from ensembl.vep.common.containers import Attribute
from ensembl.vep.common.methods import invalidate, reverse_complement
from ensembl.vep.models.genic import Transcript, Exon, Intron, Translation, CDS, UTR
from ensembl.vep.models.sequence import CoordSystem, Slice


class DummySlice(Slice):
    """Slice with a sequence supplied through the unsafe test hook."""

    def __init__(self, sequence: str, feature_start: int) -> None:
        padded_sequence = ("N" * (feature_start - 1)) + sequence
        super().__init__(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            start=1,
            end=len(padded_sequence),
            strand=Strand.FORWARD,
        )
        with self.unsafe() as unsafe_slice:
            setattr(unsafe_slice, "_seq", padded_sequence)


@pytest.fixture
def exon() -> Exon:
    """Exon fixture for testing."""
    return Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001")


@pytest.fixture
def intron() -> Intron:
    """Intron fixture for testing."""
    return Intron(start=11, end=19, strand=Strand.FORWARD)


@pytest.fixture
def transcript() -> Transcript:
    """Transcript fixture for testing."""
    return Transcript(start=1, end=50, strand=Strand.FORWARD, stable_id="ENST00000000001")


@pytest.fixture
def translation() -> Translation:
    """Translation fixture for testing."""
    return Translation(start=1, end=50, stable_id="ENSP00000000001")


@pytest.fixture
def spliced_exons_forward() -> list[Exon]:
    """Create exons backed by dummy slices with known sequences."""
    return [
        Exon(
            start=1,
            end=10,
            strand=Strand.FORWARD,
            stable_id="ENSE00000000001",
            slice=DummySlice("A" * 10, feature_start=1),
        ),
        Exon(
            start=20,
            end=30,
            strand=Strand.FORWARD,
            stable_id="ENSE00000000002",
            slice=DummySlice("C" * 11, feature_start=20),
        ),
        Exon(
            start=40,
            end=50,
            strand=Strand.FORWARD,
            stable_id="ENSE00000000003",
            slice=DummySlice("G" * 11, feature_start=40),
        ),
    ]


@pytest.fixture
def reverse_transcript() -> Transcript:
    """Transcript fixture for testing on reverse strand."""
    return Transcript(start=1, end=50, strand=Strand.REVERSE, stable_id="ENST00000000001")


@pytest.fixture
def spliced_exons_reverse() -> list[Exon]:
    """Create exons on the reverse strand"""
    return [
        Exon(
            start=40,
            end=50,
            strand=Strand.REVERSE,
            stable_id="ENSE00000000001",
            slice=DummySlice("G" * 11, feature_start=40),
        ),
        Exon(
            start=20,
            end=30,
            strand=Strand.REVERSE,
            stable_id="ENSE00000000002",
            slice=DummySlice("C" * 11, feature_start=20),
        ),
        Exon(
            start=1,
            end=10,
            strand=Strand.REVERSE,
            stable_id="ENSE00000000003",
            slice=DummySlice("A" * 10, feature_start=1),
        ),
    ]


@pytest.fixture
def spliced_exons_no_sequence() -> list[Exon]:
    """Create exons with no sequences."""
    return [
        Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001"),
        Exon(start=20, end=30, strand=Strand.FORWARD, stable_id="ENSE00000000002"),
        Exon(start=40, end=50, strand=Strand.FORWARD, stable_id="ENSE00000000003"),
    ]


@pytest.fixture
def coding_transcript() -> Transcript:
    """Create a two-exon transcript with a known coding and translated sequence."""
    genomic_sequence = "ATGGGANNNCCCTAA"
    transcript_slice = DummySlice(genomic_sequence, feature_start=1)
    transcript = Transcript(
        start=1,
        end=len(genomic_sequence),
        strand=Strand.FORWARD,
        stable_id="ENST_CODING",
        slice=transcript_slice,
    )
    transcript.spliced_exons = [
        Exon(start=1, end=6, strand=Strand.FORWARD, stable_id="ENSE_CODING_1", slice=transcript_slice),
        Exon(
            start=10,
            end=15,
            strand=Strand.FORWARD,
            stable_id="ENSE_CODING_2",
            slice=transcript_slice,
        ),
    ]
    translation = Translation(start=1, end=15, stable_id="ENSP_CODING")
    translation.transcript = transcript
    transcript.seed_attrs(translation=translation)
    return transcript


class TestTranscriptInitialization:
    """Test Transcript initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = Transcript(
                1, end=100, strand=Strand.FORWARD, stable_id="ENST00000000001"  # type: ignore[misc]
            )


class TestIntrons:
    """Tests for the lazily derived ``introns`` property.

    One Intron per consecutive exon pair, in transcript order, including a pair with no genomic
    gap between them. Parity with the builder path is asserted separately, in
    ``tests/sources/file/gxf/test_transcript.py``.
    """

    def test_no_exons_has_no_introns(self, transcript: Transcript) -> None:
        """A transcript with no exons has no introns."""
        assert transcript.introns == []

    def test_single_exon_has_no_introns(self, transcript: Transcript) -> None:
        """A single-exon transcript has no introns."""
        transcript.spliced_exons = [Exon(start=1, end=50, strand=Strand.FORWARD, stable_id="ENSE1")]

        assert transcript.introns == []

    def test_forward_strand_introns_span_the_gaps(
        self, transcript: Transcript, spliced_exons_forward: list[Exon]
    ) -> None:
        """Forward strand: introns fill the gaps between consecutive exons, in order."""
        transcript.spliced_exons = spliced_exons_forward

        introns = transcript.introns

        assert [(intron.start, intron.end) for intron in introns] == [(11, 19), (31, 39)]
        assert (introns[0].prev_exon, introns[0].next_exon) == (
            spliced_exons_forward[0],
            spliced_exons_forward[1],
        )
        assert (introns[1].prev_exon, introns[1].next_exon) == (
            spliced_exons_forward[1],
            spliced_exons_forward[2],
        )
        assert all(intron.strand == Strand.FORWARD for intron in introns)

    def test_reverse_strand_introns_span_the_gaps(
        self, reverse_transcript: Transcript, spliced_exons_reverse: list[Exon]
    ) -> None:
        """Reverse strand: exons are given in transcript (descending genomic) order, and the
        intron's genomic span is still computed low-to-high."""
        reverse_transcript.spliced_exons = spliced_exons_reverse

        introns = reverse_transcript.introns

        assert [(intron.start, intron.end) for intron in introns] == [(31, 39), (11, 19)]
        # prev_exon/next_exon always follow transcript order, regardless of strand.
        assert (introns[0].prev_exon, introns[0].next_exon) == (
            spliced_exons_reverse[0],
            spliced_exons_reverse[1],
        )
        assert (introns[1].prev_exon, introns[1].next_exon) == (
            spliced_exons_reverse[1],
            spliced_exons_reverse[2],
        )

    def test_adjacent_exons_produce_a_zero_length_intron(self) -> None:
        """Exons with no genomic gap between them still get an intron, with start > end."""
        exon_1 = Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE1")
        exon_2 = Exon(start=11, end=20, strand=Strand.FORWARD, stable_id="ENSE2")
        transcript = Transcript(start=1, end=20, strand=Strand.FORWARD, stable_id="ENST_ADJACENT")
        transcript.spliced_exons = [exon_1, exon_2]

        (intron,) = transcript.introns

        assert (intron.start, intron.end) == (11, 10)

    def test_introns_invalidated_when_spliced_exons_reassigned(self) -> None:
        """Reassigning spliced_exons clears the cached introns."""
        transcript = Transcript(start=1, end=100, strand=Strand.FORWARD, stable_id="ENST_INVALIDATE")
        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE1"),
            Exon(start=20, end=30, strand=Strand.FORWARD, stable_id="ENSE2"),
        ]
        assert len(transcript.introns) == 1

        transcript.spliced_exons = [Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE1")]

        assert transcript.introns == []

    def test_seeding_overrides_the_derived_value(self) -> None:
        """Seeding introns bypasses derivation entirely, for callers such as the synthetic corpus
        that build their own (e.g. skipping zero-gap introns). Seeding after spliced_exons is set
        matters: spliced_exons' own on_setattr hook invalidates _introns, so seeding first would
        be discarded."""
        transcript = Transcript(start=1, end=100, strand=Strand.FORWARD, stable_id="ENST_OVERRIDE")
        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE1"),
            Exon(start=11, end=20, strand=Strand.FORWARD, stable_id="ENSE2"),
        ]
        custom_intron = Intron(start=100, end=200, strand=Strand.FORWARD)

        transcript.seed_attrs(introns=[custom_intron])

        assert transcript.introns == [custom_intron]


class TestTranscriptProperties:  # pylint: disable=too-many-public-methods
    """Test Transcript properties."""

    def test_core_properties(self, transcript: Transcript) -> None:
        """Test that Transcript core properties are set correctly."""
        assert transcript.start == 1
        assert transcript.end == 50
        assert transcript.strand == Strand.FORWARD
        assert transcript.stable_id == "ENST00000000001"

    def test_externally_settable_properties(
        self, transcript: Transcript, exon: Exon, intron: Intron, translation: Translation
    ) -> None:
        """Test that externally settable properties can be set and accessed."""
        assert transcript.spliced_exons == []
        transcript.spliced_exons = [exon]
        assert transcript.spliced_exons == [exon]

        assert transcript.introns == []
        transcript.seed_attrs(introns=[intron])
        assert transcript.introns == [intron]

        assert transcript.translation is None
        transcript.seed_attrs(translation=translation)
        assert transcript.translation == translation

        assert transcript.alternative_translations == []
        transcript.alternative_translations = [translation]
        assert transcript.alternative_translations == [translation]

        assert transcript.attributes == []
        transcript.attributes = [Attribute(name="test_attribute", value="test_value")]
        assert transcript.attributes == [Attribute(name="test_attribute", value="test_value")]

    def test_first_coding_exon_index(
        self, transcript: Transcript, spliced_exons_forward: list[Exon], translation: Translation
    ) -> None:
        """Assert first_coding_exon_index is calculated correctly."""

        # With no seeded UTR/CDS and translation, first_coding_exon_index is None.
        assert transcript.first_coding_exon_index is None

        transcript.spliced_exons = spliced_exons_forward
        assert transcript.first_coding_exon_index is None

        transcript.seed_attrs(translation=translation)
        assert transcript.first_coding_exon_index is None

        transcript.seed_and_invalidate_attrs(translation=translation)
        assert transcript.first_coding_exon_index == 0

        # set coding start to somewhere outside the exons, so that _first_coding_exon_index
        # can't find the index for the first coding exon
        transcript.seed_and_invalidate_attrs(cds=CDS(start=55, end=80, strand=Strand.FORWARD))
        assert transcript.first_coding_exon_index is None

    def test_last_coding_exon_index(
        self, transcript: Transcript, spliced_exons_forward: list[Exon], translation: Translation
    ) -> None:
        """Assert last_coding_exon_index is calculated correctly."""

        # With no seeded UTR/CDS and translation, last_coding_exon_index is None.
        assert transcript.last_coding_exon_index is None

        transcript.spliced_exons = spliced_exons_forward
        assert transcript.last_coding_exon_index is None

        transcript.seed_attrs(translation=translation)
        assert transcript.last_coding_exon_index is None

        transcript.seed_and_invalidate_attrs(translation=translation)
        assert transcript.last_coding_exon_index == 2

        # set coding end to somewhere outside the exons, so that _last_coding_exon_index
        # can't find the index for the first coding exon
        transcript.seed_and_invalidate_attrs(cds=CDS(start=55, end=80, strand=Strand.FORWARD))
        assert transcript.last_coding_exon_index is None

    def test_first_and_last_coding_exon_index_are_directly_seedable(
        self, transcript: Transcript, spliced_exons_forward: list[Exon]
    ) -> None:
        """first_coding_exon_index/last_coding_exon_index (and exon_coding_start/exon_coding_end)
        can be seeded directly."""
        transcript.spliced_exons = spliced_exons_forward
        transcript.seed_attrs(first_coding_exon_index=1, last_coding_exon_index=2)
        assert transcript.first_coding_exon_index == 1
        assert transcript.last_coding_exon_index == 2

    def test_exon_coding_start(
        self, transcript: Transcript, spliced_exons_forward: list[Exon], translation: Translation
    ) -> None:
        """Test that exon_coding_start is calculated correctly."""
        assert transcript.exon_coding_start is None

        transcript.seed_attrs(translation=translation)
        transcript.spliced_exons = spliced_exons_forward
        assert transcript.exon_coding_start == 1

    def test_exon_coding_end(
        self, transcript: Transcript, spliced_exons_forward: list[Exon], translation: Translation
    ) -> None:
        """Test that exon_coding_end is calculated correctly."""
        assert transcript.exon_coding_end is None

        transcript.seed_attrs(translation=translation)
        transcript.spliced_exons = spliced_exons_forward
        assert transcript.exon_coding_end == 11

    def test_exon_coding_start_and_end_are_directly_seedable(
        self, transcript: Transcript, spliced_exons_forward: list[Exon]
    ) -> None:
        """exon_coding_start/exon_coding_end can be seeded directly."""
        transcript.spliced_exons = spliced_exons_forward
        transcript.seed_attrs(exon_coding_start=3, exon_coding_end=6)
        assert transcript.exon_coding_start == 3
        assert transcript.exon_coding_end == 6

    def test_lengths_with_cds_phase(self, transcript: Transcript) -> None:
        """Test cds lengths are calculated properly when a CDS has phase other than ZERO."""
        # Create a two exon transcript where the CDS starts in ExonPhase.ONE
        genomic_sequence = "GGGTGGGAGANNNNNNNNNCAAGTAAGGGG"
        transcript_slice = DummySlice(genomic_sequence, feature_start=1)
        exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001", slice=transcript_slice),
            Exon(
                start=20, end=30, strand=Strand.FORWARD, stable_id="ENSE00000000002", slice=transcript_slice
            ),
        ]
        exons[0].phase = ExonPhase.ONE
        transcript.spliced_exons = exons
        transcript.seed_attrs(
            first_coding_exon_index=0, last_coding_exon_index=1, exon_coding_start=4, exon_coding_end=7
        )
        # compare the different ways of calculating the CDS length
        cds_length = transcript.cds_length
        cds = transcript.cds
        assert cds is not None
        assert cds.seq is not None
        cds_seq_len = len(cds.seq)
        cds_nucleotide_length = cds.nucleotide_length
        assert transcript.cds_mapper is not None
        mapper_length = transcript.cds_mapper.length
        assert mapper_length is not None
        # mapper length is one larger than all the other lengths
        assert cds_length == cds_seq_len == cds_nucleotide_length == mapper_length - 1

        cds_protein_length = cds.protein_length if cds else None
        assert cds_protein_length == cds_length // 3

    def test_cds_cdna_coding_start(
        self, transcript: Transcript, translation: Translation, spliced_exons_forward: list[Exon]
    ) -> None:
        """Test that cds_cdna_coding_start is calculated correctly."""
        assert transcript.cdna_coding_start is None

        transcript.seed_attrs(translation=translation)
        transcript.spliced_exons = spliced_exons_forward
        assert transcript.cdna_coding_start == 1

    def test_cds_cdna_coding_start_not_in_first_exon(
        self, transcript: Transcript, spliced_exons_forward: list[Exon]
    ) -> None:
        """Test cdna_coding_start when cds doesn't start in the first exon"""
        transcript.spliced_exons = spliced_exons_forward
        transcript.seed_and_invalidate_attrs(cds=CDS(start=25, end=50, strand=Strand.FORWARD))
        assert transcript.cdna_coding_start == 16

    def test_cds_cdna_coding_end(
        self, transcript: Transcript, translation: Translation, spliced_exons_forward: list[Exon]
    ) -> None:
        """Test that cds_cdna_coding_end is calculated correctly."""
        assert transcript.cdna_coding_end is None

        transcript.seed_attrs(translation=translation)
        transcript.spliced_exons = spliced_exons_forward
        assert transcript.cdna_coding_end == 32

    def test_cds_cdna_coding_end_not_in_first_exon(
        self, transcript: Transcript, spliced_exons_forward: list[Exon]
    ) -> None:
        """Test cdna_coding_end when cds doesn't end in the first exon"""
        transcript.spliced_exons = spliced_exons_forward
        transcript.seed_and_invalidate_attrs(cds=CDS(start=1, end=25, strand=Strand.FORWARD))
        assert transcript.cdna_coding_end == 16

    def test_coding_start_using_cds(self, transcript: Transcript) -> None:
        """Test that coding_start is calculated correctly from a CDS."""
        assert transcript.coding_start is None

        transcript.seed_and_invalidate_attrs(cds=CDS(start=6, end=50, strand=Strand.FORWARD))
        assert transcript.coding_start == 6

    def test_coding_start_using_five_prime_utr(self, transcript: Transcript) -> None:
        """Test that coding_start is calculated correctly from a 5' UTR."""
        assert transcript.coding_start is None

        transcript.seed_and_invalidate_attrs(five_prime_utr=UTR(start=1, end=5, strand=Strand.FORWARD))
        assert transcript.coding_start == 6

    def test_coding_start_using_translation(
        self, transcript: Transcript, translation: Translation, spliced_exons_forward: list[Exon]
    ) -> None:
        """Test that coding_start is calculated correctly from a translation."""
        assert transcript.coding_start is None

        transcript.seed_and_invalidate_attrs(translation=translation)
        assert transcript.coding_start == 1

        transcript.spliced_exons = spliced_exons_forward
        assert transcript.coding_start == 1

    def test_coding_start_using_first_coding_exon(
        self, transcript: Transcript, spliced_exons_forward: list[Exon]
    ) -> None:
        """Test that coding_start is calculated correctly from the first coding exon index
        and exon_coding_start."""
        assert transcript.coding_start is None

        transcript.seed_and_invalidate_attrs(first_coding_exon_index=0, exon_coding_start=6)
        assert transcript.coding_start is None

        transcript.spliced_exons = spliced_exons_forward
        transcript.seed_and_invalidate_attrs(first_coding_exon_index=0, exon_coding_start=6)
        assert transcript.coding_start == 6

    def test_coding_start_is_none_for_incomplete_seeded_exon_boundary(self) -> None:
        """A seeded exon-relative coding start still needs exons and a complete boundary."""
        transcript = Transcript(start=1, end=50, strand=Strand.FORWARD, stable_id="T")
        transcript.seed_attrs(first_coding_exon_index=0, exon_coding_start=6)
        assert transcript.coding_start is None

        transcript = Transcript(start=1, end=50, strand=Strand.FORWARD, stable_id="T")
        transcript.spliced_exons = [Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="E0")]
        transcript.seed_attrs(first_coding_exon_index=None, exon_coding_start=6)
        assert transcript.coding_start is None

    def test_coding_end_using_cds(self, transcript: Transcript) -> None:
        """Test that coding_end is calculated correctly using the CDS."""
        assert transcript.coding_end is None

        transcript.seed_and_invalidate_attrs(cds=CDS(start=6, end=50, strand=Strand.FORWARD))
        assert transcript.coding_end == 50

    def test_coding_end_using_five_prime_utr(self, transcript: Transcript) -> None:
        """Test that coding_end is calculated correctly from a 3' UTR."""
        assert transcript.coding_end is None

        transcript.seed_and_invalidate_attrs(three_prime_utr=UTR(start=51, end=100, strand=Strand.FORWARD))
        assert transcript.coding_end == 50

    def test_coding_end_using_translation(
        self, transcript: Transcript, translation: Translation, spliced_exons_forward: list[Exon]
    ) -> None:
        """Test that coding_end is calculated correctly."""
        assert transcript.coding_end is None

        transcript.seed_and_invalidate_attrs(translation=translation)
        assert transcript.coding_end == 50

        transcript.spliced_exons = spliced_exons_forward
        assert transcript.coding_end == 50

    def test_coding_end_using_last_coding_exon(
        self, transcript: Transcript, spliced_exons_forward: list[Exon]
    ) -> None:
        """Test that coding_end is calculated correctly from the last coding exon index
        and exon_coding_end."""
        assert transcript.coding_end is None

        transcript.seed_and_invalidate_attrs(last_coding_exon_index=0, exon_coding_end=45)
        assert transcript.coding_end is None

        transcript.spliced_exons = spliced_exons_forward
        transcript.seed_and_invalidate_attrs(last_coding_exon_index=0, exon_coding_end=45)
        assert transcript.coding_end == 45

    def test_coding_end_is_none_for_incomplete_seeded_exon_boundary(self) -> None:
        """A seeded exon-relative coding end still needs exons and a complete boundary."""
        transcript = Transcript(start=1, end=50, strand=Strand.FORWARD, stable_id="T")
        transcript.seed_attrs(last_coding_exon_index=0, exon_coding_end=6)
        assert transcript.coding_end is None

        transcript = Transcript(start=1, end=50, strand=Strand.FORWARD, stable_id="T")
        transcript.spliced_exons = [Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="E0")]
        transcript.seed_attrs(last_coding_exon_index=None, exon_coding_end=6)
        assert transcript.coding_end is None

    def test_cds(
        self, transcript: Transcript, translation: Translation, spliced_exons_forward: list[Exon]
    ) -> None:
        """Test that cds is calculated correctly."""
        assert transcript.cds is None

        transcript.seed_and_invalidate_attrs(translation=translation)
        translation.transcript = transcript

        # Mock spliced_seq for testing
        transcript.spliced_exons = spliced_exons_forward

        cds = CDS(start=1, end=50, strand=Strand.FORWARD)
        cds.seq = "A" * 10 + "C" * 11 + "G" * 11
        cds.transcript = transcript
        cds.translation = translation

        assert transcript.cds == cds

    def test_exon_with_no_sequence(
        self,
        transcript: Transcript,
        spliced_exons_no_sequence: list[Exon],
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Test that cds is calculated correctly."""
        assert transcript.cds is None
        transcript.spliced_exons = spliced_exons_no_sequence
        transcript.seed_and_invalidate_attrs(
            first_coding_exon_index=0,
            last_coding_exon_index=2,
            exon_coding_start=6,
            exon_coding_end=6,
        )

        with caplog.at_level(logging.WARNING):
            cds = transcript.cds
            assert transcript.cds is not None
            assert (
                "Exon ENSE00000000001 has no sequence information" in caplog.text
            ), "Incorrect warning in cds when exons have no sequence"
            assert cds.seq == "NNNNNNNNNNNNNNNNNNNNNN"

        with caplog.at_level(logging.WARNING):
            spliced_seq = transcript.spliced_seq
            assert (
                "Exon ENSE00000000001 has no sequence information" in caplog.text
            ), "Incorrect warning in spliced_seq when exons have no sequence"
            assert spliced_seq == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"

    @pytest.mark.parametrize("strand", [Strand.FORWARD, Strand.REVERSE])
    def test_cds_when_coding_starts_in_a_later_exon(self, strand: Strand) -> None:
        """CDS bounds are correct when the first coding exon is not exon 0, on both strands."""
        dummy_slice = DummySlice("A" * 60, feature_start=1)
        exons = [
            Exon(start=1, end=10, strand=strand, stable_id="E0", slice=dummy_slice),
            Exon(start=20, end=31, strand=strand, stable_id="E1", slice=dummy_slice),
            Exon(start=40, end=49, strand=strand, stable_id="E2", slice=dummy_slice),
        ]
        transcript = Transcript(start=1, end=60, strand=strand, stable_id="T", slice=dummy_slice)
        # spliced_exons is transcript order: 5'-most first, so reversed on the reverse strand.
        transcript.spliced_exons = exons if strand == Strand.FORWARD else list(reversed(exons))
        transcript.seed_attrs(
            first_coding_exon_index=1,
            last_coding_exon_index=2,
            exon_coding_start=3,
            exon_coding_end=6,
        )

        cds = transcript.cds
        assert cds is not None
        if strand == Strand.FORWARD:
            # first coding exon E1 (20-31), 3rd base in -> 22; last E2 (40-49), 6th base -> 45.
            assert transcript.cds_genomic_spans == [(22, 31), (40, 45)]
            assert (cds.start, cds.end) == (22, 45)
        else:
            # transcript order is [E2, E1, E0]; first coding exon E1 (20-31) counted from its own
            # 5' end (31) -> 29; last coding exon E0 (1-10) counted from 10 -> 5. Ascending: 5-29.
            assert transcript.cds_genomic_spans == [(20, 29), (5, 10)]
            assert (cds.start, cds.end) == (5, 29)

    @pytest.mark.parametrize("strand", [Strand.FORWARD, Strand.REVERSE])
    def test_cds_when_coding_ends_in_an_earlier_exon(self, strand: Strand) -> None:
        """CDS bounds are correct when the last coding exon is not the final exon, both strands."""
        dummy_slice = DummySlice("A" * 60, feature_start=1)
        exons = [
            Exon(start=1, end=10, strand=strand, stable_id="E0", slice=dummy_slice),
            Exon(start=20, end=31, strand=strand, stable_id="E1", slice=dummy_slice),
            Exon(start=40, end=49, strand=strand, stable_id="E2", slice=dummy_slice),
        ]
        transcript = Transcript(start=1, end=60, strand=strand, stable_id="T", slice=dummy_slice)
        # spliced_exons is transcript order: 5'-most first, so reversed on the reverse strand.
        transcript.spliced_exons = exons if strand == Strand.FORWARD else list(reversed(exons))
        transcript.seed_attrs(
            first_coding_exon_index=0,
            last_coding_exon_index=1,
            exon_coding_start=3,
            exon_coding_end=6,
        )

        cds = transcript.cds
        assert cds is not None
        if strand == Strand.FORWARD:
            # first coding exon E0 (1-10), 3rd base in -> 3; last coding exon E1 (20-31), 6th
            # base -> 25. E2 is past the last coding exon and contributes nothing.
            assert (cds.start, cds.end) == (3, 25)
        else:
            # transcript order is [E2, E1, E0]; first coding exon E2 (40-49) counted from its own
            # 5' end (49) -> 47; last coding exon E1 (20-31) counted from 31 -> 26. E0 is unused.
            # Ascending: 26-47.
            assert (cds.start, cds.end) == (26, 47)
        # 8 bases from the first coding exon (offset 3 to its end) + 6 from the last.
        assert cds.seq is not None
        assert len(cds.seq) == 14

    @staticmethod
    def _transcript_with_utrs(
        strand: Strand, *, translation_start: int = 3, translation_end: int = 45
    ) -> Transcript:
        """Build a coding transcript whose exons are ordered in transcript orientation."""
        exons = [
            Exon(start=1, end=10, strand=strand, stable_id="E0"),
            Exon(start=20, end=31, strand=strand, stable_id="E1"),
            Exon(start=40, end=49, strand=strand, stable_id="E2"),
        ]
        transcript = Transcript(start=1, end=49, strand=strand, stable_id="T")
        transcript.spliced_exons = exons if strand == Strand.FORWARD else list(reversed(exons))
        transcript.seed_attrs(
            translation=Translation(
                start=translation_start,
                end=translation_end,
                stable_id="P",
            )
        )
        return transcript

    @pytest.mark.parametrize(
        ("strand", "expected_start", "expected_end"),
        [
            (Strand.FORWARD, 1, 2),
            (Strand.REVERSE, 48, 49),
        ],
    )
    def test_five_prime_utr(self, strand: Strand, expected_start: int, expected_end: int) -> None:
        """The 5' UTR is the envelope before the coding start in transcript orientation."""
        transcript = self._transcript_with_utrs(strand)

        five_prime_utr = transcript.five_prime_utr

        assert five_prime_utr is not None
        assert (five_prime_utr.start, five_prime_utr.end) == (expected_start, expected_end)
        assert five_prime_utr.strand == strand
        assert five_prime_utr.transcript is transcript

    @pytest.mark.parametrize(
        ("strand", "expected_start", "expected_end"),
        [
            (Strand.FORWARD, 46, 49),
            (Strand.REVERSE, 1, 4),
        ],
    )
    def test_three_prime_utr(self, strand: Strand, expected_start: int, expected_end: int) -> None:
        """The 3' UTR is the envelope after the coding end in transcript orientation."""
        transcript = self._transcript_with_utrs(strand)

        three_prime_utr = transcript.three_prime_utr

        assert three_prime_utr is not None
        assert (three_prime_utr.start, three_prime_utr.end) == (expected_start, expected_end)
        assert three_prime_utr.strand == strand
        assert three_prime_utr.transcript is transcript

    @pytest.mark.parametrize("strand", [Strand.FORWARD, Strand.REVERSE])
    def test_utrs_are_none_when_coding_region_reaches_transcript_ends(self, strand: Strand) -> None:
        """A coding region spanning the transcript leaves no UTR on either side."""
        transcript = self._transcript_with_utrs(strand, translation_start=1, translation_end=49)

        assert transcript.five_prime_utr is None
        assert transcript.three_prime_utr is None

    def test_reverse_strand_cds_read_first_keeps_coding_boundaries(self) -> None:
        """Reading the derived CDS before coding_start/end does not flip them on the reverse strand."""
        dummy_slice = DummySlice("ACGT" * 25, feature_start=1)
        transcript = Transcript(start=1, end=100, strand=Strand.REVERSE, stable_id="T", slice=dummy_slice)
        transcript.spliced_exons = [
            Exon(start=60, end=100, strand=Strand.REVERSE, stable_id="E0", slice=dummy_slice),
            Exon(start=1, end=40, strand=Strand.REVERSE, stable_id="E1", slice=dummy_slice),
        ]
        transcript.seed_attrs(
            first_coding_exon_index=0, last_coding_exon_index=1, exon_coding_start=11, exon_coding_end=31
        )

        cds = transcript.cds

        assert cds is not None
        assert (cds.start, cds.end) == (10, 90)
        assert (transcript.coding_start, transcript.coding_end) == (90, 10)
        assert transcript.five_prime_utr is not None
        assert (transcript.five_prime_utr.start, transcript.five_prime_utr.end) == (91, 100)
        assert transcript.three_prime_utr is not None
        assert (transcript.three_prime_utr.start, transcript.three_prime_utr.end) == (1, 9)

    @pytest.mark.parametrize(
        ("strand", "expected_seq"),
        [(Strand.FORWARD, "ATGAAATAA"), (Strand.REVERSE, "ATGAAATAA")],
    )
    def test_single_exon_cds_seq_excludes_both_utrs(self, strand: Strand, expected_seq: str) -> None:
        """A CDS inside one exon with UTR on both sides yields exactly the coding bases."""
        transcript_sense = "CC" + expected_seq + "GGG"
        genomic = transcript_sense if strand == Strand.FORWARD else reverse_complement(transcript_sense)
        dummy_slice = DummySlice(genomic, feature_start=1)
        transcript = Transcript(start=1, end=14, strand=strand, stable_id="T", slice=dummy_slice)
        transcript.spliced_exons = [Exon(start=1, end=14, strand=strand, stable_id="E0", slice=dummy_slice)]
        transcript.seed_attrs(
            first_coding_exon_index=0, last_coding_exon_index=0, exon_coding_start=3, exon_coding_end=11
        )

        cds = transcript.cds

        assert cds is not None
        assert (cds.start, cds.end) == ((3, 11) if strand == Strand.FORWARD else (4, 12))
        assert cds.seq == expected_seq

    def test_cds_is_none_when_coding_offsets_do_not_fit_exons(self) -> None:
        """An exon_coding_start past the first coding exon's end yields no CDS."""
        dummy_slice = DummySlice("A" * 40, feature_start=1)
        transcript = Transcript(start=1, end=40, strand=Strand.FORWARD, stable_id="T", slice=dummy_slice)
        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="E0", slice=dummy_slice),
            Exon(start=20, end=40, strand=Strand.FORWARD, stable_id="E1", slice=dummy_slice),
        ]
        transcript.seed_attrs(
            first_coding_exon_index=0, last_coding_exon_index=1, exon_coding_start=12, exon_coding_end=5
        )

        assert transcript.cds is None

    def test_cds_genomic_spans_is_none_when_offsets_clip_away_all_blocks(self) -> None:
        """Inconsistent coding offsets can leave no valid genomic CDS block."""
        transcript = Transcript(start=1, end=10, strand=Strand.FORWARD, stable_id="T")
        transcript.spliced_exons = [Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="E0")]
        transcript.seed_attrs(
            first_coding_exon_index=0,
            last_coding_exon_index=0,
            exon_coding_start=12,
            exon_coding_end=11,
        )

        assert transcript.cds_genomic_spans is None

    def test_five_prime_utr_is_None(self, transcript: Transcript, spliced_exons_forward: list[Exon]) -> None:
        """Test five_prime_utr is None if there are no exons or no coding_start"""
        assert transcript.five_prime_utr is None
        transcript.spliced_exons = spliced_exons_forward
        assert transcript.five_prime_utr is None

    def test_three_prime_utr_is_None(self, transcript: Transcript, spliced_exons_forward: list[Exon]) -> None:
        """Test three_prime_utr is None if there are no exons or no coding_end"""
        assert transcript.three_prime_utr is None
        transcript.spliced_exons = spliced_exons_forward
        assert transcript.three_prime_utr is None

    @pytest.mark.parametrize(
        ("strand", "expected_five_prime", "expected_three_prime"),
        [
            (Strand.FORWARD, [(1, 2)], [(46, 49)]),
            (Strand.REVERSE, [(48, 49)], [(1, 4)]),
        ],
    )
    def test_utr_genomic_spans(
        self,
        strand: Strand,
        expected_five_prime: list[tuple[int, int]],
        expected_three_prime: list[tuple[int, int]],
    ) -> None:
        """UTR genomic spans include transcript-order exon blocks clipped to the UTR."""
        transcript = self._transcript_with_utrs(strand)

        assert transcript.five_prime_utr_genomic_spans == expected_five_prime
        assert transcript.three_prime_utr_genomic_spans == expected_three_prime

    def test_utr_genomic_spans_include_whole_non_coding_exons(self) -> None:
        """UTR genomic spans include whole exons before/after the CDS plus clipped coding exons."""
        transcript = Transcript(start=1, end=49, strand=Strand.FORWARD, stable_id="T")
        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="E0"),
            Exon(start=20, end=31, strand=Strand.FORWARD, stable_id="E1"),
            Exon(start=40, end=49, strand=Strand.FORWARD, stable_id="E2"),
        ]
        transcript.seed_attrs(cds=CDS(start=25, end=45, strand=Strand.FORWARD))

        assert transcript.five_prime_utr_genomic_spans == [(1, 10), (20, 24)]
        assert transcript.three_prime_utr_genomic_spans == [(46, 49)]

    def test_utr_genomic_spans_are_none_without_genomic_spans(self) -> None:
        """A transcript with no cds_genomic_spans has no UTR genomic spans."""
        transcript = Transcript(start=1, end=10, strand=Strand.FORWARD, stable_id="T")
        transcript.spliced_exons = [Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="E0")]

        assert transcript.five_prime_utr_genomic_spans is None
        assert transcript.three_prime_utr_genomic_spans is None

    def test_utr_genomic_spans_are_none_without_utr_bases(self) -> None:
        """A CDS spanning the full spliced transcript leaves no UTR genomic spans."""
        transcript = Transcript(start=1, end=10, strand=Strand.FORWARD, stable_id="T")
        transcript.spliced_exons = [Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="E0")]
        transcript.seed_attrs(
            first_coding_exon_index=0,
            last_coding_exon_index=0,
            exon_coding_start=1,
            exon_coding_end=10,
        )

        assert transcript.five_prime_utr_genomic_spans is None
        assert transcript.three_prime_utr_genomic_spans is None

    def test_translation(self, coding_transcript: Transcript) -> None:
        """Test that a translation is derived from resolved coding boundaries."""
        invalidate(coding_transcript, "_translation")
        coding_transcript.seed_attrs(
            first_coding_exon_index=0,
            last_coding_exon_index=1,
            exon_coding_start=1,
            exon_coding_end=6,
        )

        translation = coding_transcript.translation

        assert translation is not None
        assert (translation.start, translation.end) == (1, 15)
        assert translation.stable_id == "ENST_CODING_derived_translation"
        assert (translation.exon_coding_start, translation.exon_coding_end) == (1, 6)
        assert translation.seq == "MGP*"

    def test_translation_reverse_strand(self, reverse_transcript: Transcript) -> None:
        """Test the translation on a transcript on the reverse strand"""
        reverse_transcript.seed_and_invalidate_attrs(cds=CDS(start=25, end=50, strand=Strand.REVERSE))
        translation = reverse_transcript.translation
        assert translation is not None
        assert (translation.start, translation.end) == (1, 26)

    def test_spliced_seq(self, transcript: Transcript, spliced_exons_forward: list[Exon]) -> None:
        """Test that spliced_seq is calculated correctly."""
        assert transcript.spliced_seq is None

        transcript.spliced_exons = spliced_exons_forward

        assert transcript.spliced_seq == "A" * 10 + "C" * 11 + "G" * 11

    def test_coding_seq(self, coding_transcript: Transcript) -> None:
        """Test that coding_seq is calculated correctly."""
        assert coding_transcript.coding_sequence == "ATGGGACCCTAA"

    def test_coding_seq_returns_none_if_slice_seq_is_none(self, transcript: Transcript) -> None:
        """Test that coding_seq is None if the transcript slice has no sequence."""
        slice_no_seq = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            start=1,
            end=100,
            strand=Strand.FORWARD,
            sequence_source=None,
        )
        transcript.slice = slice_no_seq
        assert transcript.coding_sequence is None

    @pytest.mark.parametrize("strand", [Strand.FORWARD, Strand.REVERSE])
    def test_coding_seq_no_cds_seq(self, strand: Strand) -> None:
        """Test that coding_seq is calculated without a CDS."""
        expected_seq = "ATGAAATAA"
        transcript_sense = "CC" + expected_seq + "GGG"
        genomic = transcript_sense if strand == Strand.FORWARD else reverse_complement(transcript_sense)
        dummy_slice = DummySlice(genomic, feature_start=1)
        transcript = Transcript(start=1, end=14, strand=strand, stable_id="T", slice=dummy_slice)
        transcript.spliced_exons = [Exon(start=1, end=14, strand=strand, stable_id="E0", slice=dummy_slice)]
        transcript.seed_attrs(
            first_coding_exon_index=0, last_coding_exon_index=0, exon_coding_start=3, exon_coding_end=11
        )

        assert transcript.coding_sequence == expected_seq

    @pytest.mark.parametrize("strand", [Strand.FORWARD, Strand.REVERSE])
    def test_coding_seq_uses_slice_sequence_when_seeded_cds_has_no_sequence(self, strand: Strand) -> None:
        """A seeded CDS with no seq can still derive coding sequence from transcript slice spans."""
        expected_seq = "ATGAAATAA"
        transcript_sense = "CC" + expected_seq + "GGG"
        genomic = transcript_sense if strand == Strand.FORWARD else reverse_complement(transcript_sense)
        dummy_slice = DummySlice(genomic, feature_start=1)
        transcript = Transcript(start=1, end=14, strand=strand, stable_id="T", slice=dummy_slice)
        transcript.spliced_exons = [Exon(start=1, end=14, strand=strand, stable_id="E0", slice=dummy_slice)]
        cds = CDS(
            start=3 if strand == Strand.FORWARD else 4,
            end=11 if strand == Strand.FORWARD else 12,
            strand=strand,
        )
        transcript.seed_attrs(cds=cds)

        assert transcript.coding_sequence == expected_seq

    def test_translated_seq(self, coding_transcript: Transcript) -> None:
        """Test that translated_seq is calculated correctly."""
        translation = coding_transcript.translation
        assert translation is not None
        assert translation.seq is None
        assert coding_transcript.translated_seq == "MGP*"

    def test_translated_seq_prefers_seeded_translation_sequence(self, transcript: Transcript) -> None:
        """translated_seq returns a seeded protein sequence without deriving from CDS."""
        translation = Translation(start=1, end=9, stable_id="P")
        translation.seq = "MKT"
        transcript.seed_attrs(translation=translation)

        assert transcript.translated_seq == "MKT"

    def test_translated_seq_is_none_when_no_translation_or_cds_sequence(
        self, transcript: Transcript, caplog: pytest.LogCaptureFixture
    ) -> None:
        """translated_seq warns and returns None when neither translation nor CDS sequence exists."""
        with caplog.at_level(logging.WARNING):
            assert transcript.translated_seq is None

        assert "translation and CDS sequence are both not available" in caplog.text

    def test_cdna_coding_end_uses_whole_last_coding_exon_when_end_offset_is_unknown(
        self, spliced_exons_forward: list[Exon]
    ) -> None:
        """If exon_coding_end is unknown, cdna_coding_end falls back to the full exon length."""
        transcript = Transcript(start=1, end=50, strand=Strand.FORWARD, stable_id="T")
        transcript.spliced_exons = spliced_exons_forward
        transcript.seed_attrs(last_coding_exon_index=1, exon_coding_end=None)

        assert transcript.cdna_coding_end == 21

    def test_cdna_mapper_no_cds(self, spliced_exons_forward: list[Exon]) -> None:
        """A non-coding transcript maps cDNA bases with an empty coding range."""
        transcript = Transcript(start=1, end=50, strand=Strand.FORWARD, stable_id="T")
        transcript.spliced_exons = spliced_exons_forward
        assert transcript.cdna_mapper is not None
        assert transcript.cdna_mapper.coding_start_offset == 0
        assert transcript.cdna_mapper.coding_end_offset == -1

    def test_cdna_mapper(self, spliced_exons_forward: list[Exon]) -> None:
        """A coding transcript's cDNA mapper records the coding interval."""
        transcript = Transcript(start=1, end=50, strand=Strand.FORWARD, stable_id="T")
        transcript.spliced_exons = spliced_exons_forward
        transcript.seed_and_invalidate_attrs(
            first_coding_exon_index=0,
            last_coding_exon_index=2,
            exon_coding_start=5,
            exon_coding_end=5,
        )
        assert transcript.cdna_mapper is not None
        assert transcript.cdna_mapper.coding_start_offset == 5
        assert transcript.cdna_mapper.coding_end_offset == 26


class TestTranscriptEquality:
    """Test equality comparison."""

    def test_equal_transcripts(self, transcript: Transcript) -> None:
        """Transcripts with the same stable ID and version are equal."""
        other = Transcript(
            start=transcript.start,
            end=transcript.end,
            strand=transcript.strand,
            stable_id=transcript.stable_id,
        )

        assert transcript == other

    def test_different_stable_id(self, transcript: Transcript) -> None:
        """Transcripts with different stable IDs are not equal."""
        other = Transcript(
            start=transcript.start,
            end=transcript.end,
            strand=transcript.strand,
            stable_id="ENST00000000002",
        )

        assert transcript != other

    def test_different_version(self, transcript: Transcript) -> None:
        """Transcripts with different versions are not equal."""
        other = Transcript(
            start=transcript.start,
            end=transcript.end,
            strand=transcript.strand,
            stable_id=transcript.stable_id,
        )
        other.version = 2

        assert transcript != other

    def test_different_type(self, transcript: Transcript, exon: Exon) -> None:
        """A transcript is not equal to another feature type."""
        assert transcript != exon
