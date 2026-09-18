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

import pytest

from ensembl.vep.common.enums import CoordSystemLevel, Strand
from ensembl.vep.common.containers import Attribute
from ensembl.vep.common.methods import invalidate, reverse_complement
from ensembl.vep.models.genic import Transcript, Exon, Intron, Translation, CDS
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
    return Transcript(start=1, end=100, strand=Strand.FORWARD, stable_id="ENST00000000001")


@pytest.fixture
def translation() -> Translation:
    """Translation fixture for testing."""
    return Translation(start=1, end=50, stable_id="ENSP00000000001")


@pytest.fixture
def spliced_exons() -> list[Exon]:
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


class TestTranscriptProperties:
    """Test Transcript properties."""

    def test_core_properties(self, transcript: Transcript) -> None:
        """Test that Transcript core properties are set correctly."""
        assert transcript.start == 1
        assert transcript.end == 100
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
        transcript.introns = [intron]
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
        self, transcript: Transcript, spliced_exons: list[Exon], translation: Translation
    ) -> None:
        """Assert first_coding_exon_index is calculated correctly."""

        # With no seeded UTR/CDS and translation, first_coding_exon_index is None.
        assert transcript.first_coding_exon_index is None

        transcript.spliced_exons = spliced_exons
        assert transcript.first_coding_exon_index is None

        transcript.seed_attrs(translation=translation)
        assert transcript.first_coding_exon_index is None

        transcript.seed_and_invalidate_attrs(translation=translation)
        assert transcript.first_coding_exon_index == 0

    def test_last_coding_exon_index(
        self, transcript: Transcript, spliced_exons: list[Exon], translation: Translation
    ) -> None:
        """Assert last_coding_exon_index is calculated correctly."""

        # With no seeded UTR/CDS and translation, last_coding_exon_index is None.
        assert transcript.last_coding_exon_index is None

        transcript.spliced_exons = spliced_exons
        assert transcript.last_coding_exon_index is None

        transcript.seed_attrs(translation=translation)
        assert transcript.last_coding_exon_index is None

        transcript.seed_and_invalidate_attrs(translation=translation)
        assert transcript.last_coding_exon_index == 2

    def test_first_and_last_coding_exon_index_are_directly_seedable(self, transcript: Transcript) -> None:
        """first_coding_exon_index/last_coding_exon_index (and exon_coding_start/exon_coding_end)
        can be seeded directly."""
        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001"),
            Exon(start=20, end=30, strand=Strand.FORWARD, stable_id="ENSE00000000002"),
            Exon(start=40, end=50, strand=Strand.FORWARD, stable_id="ENSE00000000003"),
        ]
        transcript.seed_attrs(first_coding_exon_index=1, last_coding_exon_index=2)
        assert transcript.first_coding_exon_index == 1
        assert transcript.last_coding_exon_index == 2

    def test_exon_coding_start(
        self, transcript: Transcript, spliced_exons: list[Exon], translation: Translation
    ) -> None:
        """Test that exon_coding_start is calculated correctly."""
        assert transcript.exon_coding_start is None

        transcript.seed_attrs(translation=translation)
        transcript.spliced_exons = spliced_exons
        assert transcript.exon_coding_start == 1

    def test_exon_coding_end(
        self, transcript: Transcript, spliced_exons: list[Exon], translation: Translation
    ) -> None:
        """Test that exon_coding_end is calculated correctly."""
        assert transcript.exon_coding_end is None

        transcript.seed_attrs(translation=translation)
        transcript.spliced_exons = spliced_exons
        assert transcript.exon_coding_end == 11

    def test_exon_coding_start_and_end_are_directly_seedable(self, transcript: Transcript) -> None:
        """exon_coding_start/exon_coding_end can be seeded directly."""
        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001"),
            Exon(start=20, end=30, strand=Strand.FORWARD, stable_id="ENSE00000000002"),
            Exon(start=40, end=50, strand=Strand.FORWARD, stable_id="ENSE00000000003"),
        ]
        transcript.seed_attrs(exon_coding_start=3, exon_coding_end=6)
        assert transcript.exon_coding_start == 3
        assert transcript.exon_coding_end == 6

    def test_cds_cdna_coding_start(self, transcript: Transcript, translation: Translation) -> None:
        """Test that cds_cdna_coding_start is calculated correctly."""
        assert transcript.cdna_coding_start is None

        transcript.seed_attrs(translation=translation)
        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001"),
            Exon(start=20, end=30, strand=Strand.FORWARD, stable_id="ENSE00000000002"),
            Exon(start=40, end=50, strand=Strand.FORWARD, stable_id="ENSE00000000003"),
        ]
        assert transcript.cdna_coding_start == 1

    def test_cds_cdna_coding_end(self, transcript: Transcript, translation: Translation) -> None:
        """Test that cds_cdna_coding_end is calculated correctly."""
        assert transcript.cdna_coding_end is None

        transcript.seed_attrs(translation=translation)
        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001"),
            Exon(start=20, end=30, strand=Strand.FORWARD, stable_id="ENSE00000000002"),
            Exon(start=40, end=50, strand=Strand.FORWARD, stable_id="ENSE00000000003"),
        ]
        assert transcript.cdna_coding_end == 32

    def test_coding_start(self, transcript: Transcript, translation: Translation) -> None:
        """Test that coding_start is calculated correctly."""
        assert transcript.coding_start is None

        # Seeding invalidates the cached None above; the translation alone fixes coding_start,
        # before any exon is known.
        transcript.seed_and_invalidate_attrs(translation=translation)
        assert transcript.coding_start == 1

        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001"),
            Exon(start=20, end=30, strand=Strand.FORWARD, stable_id="ENSE00000000002"),
            Exon(start=40, end=50, strand=Strand.FORWARD, stable_id="ENSE00000000003"),
        ]
        assert transcript.coding_start == 1

    def test_coding_end(self, transcript: Transcript, translation: Translation) -> None:
        """Test that coding_end is calculated correctly."""
        assert transcript.coding_end is None

        # Seeding invalidates the cached None above; the translation alone fixes coding_end.
        transcript.seed_and_invalidate_attrs(translation=translation)
        assert transcript.coding_end == 50

        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001"),
            Exon(start=20, end=30, strand=Strand.FORWARD, stable_id="ENSE00000000002"),
            Exon(start=40, end=50, strand=Strand.FORWARD, stable_id="ENSE00000000003"),
        ]
        assert transcript.coding_end == 50

    def test_cds(self, transcript: Transcript, translation: Translation, spliced_exons: list[Exon]) -> None:
        """Test that cds is calculated correctly."""
        assert transcript.cds is None

        transcript.seed_and_invalidate_attrs(translation=translation)
        translation.transcript = transcript

        # Mock spliced_seq for testing
        transcript.spliced_exons = spliced_exons

        cds = CDS(start=1, end=50, strand=Strand.FORWARD)
        cds.seq = "A" * 10 + "C" * 11 + "G" * 11
        cds.transcript = transcript
        cds.translation = translation

        assert transcript.cds == cds

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
            assert (cds.start, cds.end) == (22, 45)
        else:
            # transcript order is [E2, E1, E0]; first coding exon E1 (20-31) counted from its own
            # 5' end (31) -> 29; last coding exon E0 (1-10) counted from 10 -> 5. Ascending: 5-29.
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
        assert cds.seq == expected_seq
        assert (cds.start, cds.end) == ((3, 11) if strand == Strand.FORWARD else (4, 12))

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

    def test_spliced_seq(self, transcript: Transcript, spliced_exons: list[Exon]) -> None:
        """Test that spliced_seq is calculated correctly."""
        assert transcript.spliced_seq is None

        transcript.spliced_exons = spliced_exons

        assert transcript.spliced_seq == "A" * 10 + "C" * 11 + "G" * 11

    def test_coding_seq(self, coding_transcript: Transcript) -> None:
        """Test that coding_seq is calculated correctly."""
        assert coding_transcript.coding_sequence == "ATGGGACCCTAA"

    def test_translated_seq(self, coding_transcript: Transcript) -> None:
        """Test that translated_seq is calculated correctly."""
        translation = coding_transcript.translation
        assert translation is not None
        assert translation.seq is None
        assert coding_transcript.translated_seq == "MGP*"


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
