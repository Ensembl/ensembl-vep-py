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

from ensembl.vep.common import Attribute, Strand
from ensembl.vep.models.genic import Transcript, Exon, Intron, Translation, CDS


@pytest.fixture
def exon() -> Exon:
    return Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001")


@pytest.fixture
def intron() -> Intron:
    return Intron(start=11, end=19, strand=Strand.FORWARD)


@pytest.fixture
def transcript() -> Transcript:
    return Transcript(start=1, end=100, strand=Strand.FORWARD, stable_id="ENST00000000001")


@pytest.fixture
def translation() -> Translation:
    return Translation(start=10, end=90, stable_id="ENSP00000000001")


class TestTranscriptInitialization:
    """Test Transcript initialization."""

    def test_initialization_without_keywords(self) -> None:
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = Transcript(start=1, end=100, strand=Strand.FORWARD, stable_id="ENST00000000001")


class TestTranscriptProperties:
    """Test Transcript properties."""

    def test_core_properties(self, transcript: Transcript) -> None:
        """Test that Transcript core properties are set correctly."""
        assert transcript.start == 1
        assert transcript.end == 100
        assert transcript.strand == 1
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
        transcript.translation = translation
        assert transcript.translation == translation

        assert transcript.alternative_translations == []
        transcript.alternative_translations = [translation]
        assert transcript.alternative_translations == [translation]

        assert transcript.attributes == []
        transcript.attributes = [Attribute(name="test_attribute", value="test_value")]
        assert transcript.attributes == [Attribute(name="test_attribute", value="test_value")]

    def test_start_exon_index(self, transcript: Transcript) -> None:
        """Test that start_exon_index is calculated correctly."""
        assert transcript.start_exon_index == 0

    def test_end_exon_index(self, transcript: Transcript) -> None:
        """Test that end_exon_index is calculated correctly."""
        assert transcript.end_exon_index == -1

        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001"),
            Exon(start=20, end=30, strand=Strand.FORWARD, stable_id="ENSE00000000002"),
            Exon(start=40, end=50, strand=Strand.FORWARD, stable_id="ENSE00000000003"),
        ]

        assert transcript.end_exon_index == 2

    def test_spliced_seq(self, transcript: Transcript) -> None:
        """Test that spliced_seq is calculated correctly."""
        assert transcript.spliced_seq is None

        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001"),
            Exon(start=20, end=30, strand=Strand.FORWARD, stable_id="ENSE00000000002"),
            Exon(start=40, end=50, strand=Strand.FORWARD, stable_id="ENSE00000000003"),
        ]

        # Mock exon sequences for testing
        transcript.spliced_exons[0].seq = "A" * 10
        transcript.spliced_exons[1].seq = "C" * 11
        transcript.spliced_exons[2].seq = "G" * 11

        assert transcript.spliced_seq == "A" * 10 + "C" * 11 + "G" * 11

    def test_cds(self, transcript: Transcript, translation: Translation) -> None:
        """Test that cds is calculated correctly."""
        assert transcript.cds is None

        transcript.translation = translation
        translation.transcript = transcript

        # Mock spliced_seq for testing
        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001"),
            Exon(start=20, end=30, strand=Strand.FORWARD, stable_id="ENSE00000000002"),
            Exon(start=40, end=50, strand=Strand.FORWARD, stable_id="ENSE00000000003"),
        ]

        transcript.spliced_exons[0].seq = "A" * 10
        transcript.spliced_exons[1].seq = "C" * 11
        transcript.spliced_exons[2].seq = "G" * 11

        cds = CDS(start=1, end=50, strand=Strand.FORWARD)
        cds.seq = "A" * 10 + "C" * 11 + "G" * 11
        cds.transcript = transcript
        cds.translation = translation

        assert transcript.cds == cds

        # should be way to assign exon_coding_start for translation?

    def test_cds_cdna_coding_start(self, transcript: Transcript, translation: Translation) -> None:
        """Test that cds_cdna_coding_start is calculated correctly."""
        assert transcript.cdna_coding_start is None

        transcript.translation = translation
        assert transcript.cdna_coding_start == 1

        # should be way to assign start_exon_index?

    def test_cds_cdna_coding_end(self, transcript: Transcript, translation: Translation) -> None:
        """Test that cds_cdna_coding_end is calculated correctly."""
        assert transcript.cdna_coding_end is None

        transcript.translation = translation
        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001"),
            Exon(start=20, end=30, strand=Strand.FORWARD, stable_id="ENSE00000000002"),
            Exon(start=40, end=50, strand=Strand.FORWARD, stable_id="ENSE00000000003"),
        ]
        assert transcript.cdna_coding_end == 32

        # should be way to assign end_exon_index?

    def test_cds_start(self, transcript: Transcript, translation: Translation) -> None:
        """Test that cds_start is calculated correctly."""
        assert transcript.cds_start is None

        transcript.translation = translation
        assert transcript.cds_start is None

        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001"),
            Exon(start=20, end=30, strand=Strand.FORWARD, stable_id="ENSE00000000002"),
            Exon(start=40, end=50, strand=Strand.FORWARD, stable_id="ENSE00000000003"),
        ]
        assert transcript.cds_start == 1

    def test_cds_end(self, transcript: Transcript, translation: Translation) -> None:
        """Test that cds_end is calculated correctly."""
        assert transcript.cds_end is None

        transcript.translation = translation
        assert transcript.cds_end is None

        transcript.spliced_exons = [
            Exon(start=1, end=10, strand=Strand.FORWARD, stable_id="ENSE00000000001"),
            Exon(start=20, end=30, strand=Strand.FORWARD, stable_id="ENSE00000000002"),
            Exon(start=40, end=50, strand=Strand.FORWARD, stable_id="ENSE00000000003"),
        ]
        assert transcript.cds_end == 50

    def test_5_prime_utr(self, transcript: Transcript, translation: Translation) -> None:
        """Test that 5_prime_utr is calculated correctly."""
        # assert transcript.five_prime_utr is None

        # 5_prime_utr accesses spliced_seq and gets into deadlock, workaround?

    def test_3_prime_utr(self, transcript: Transcript, translation: Translation) -> None:
        """Test that 3_prime_utr is calculated correctly."""
        # assert transcript.three_prime_utr is None

        # 3_prime_utr accesses spliced_seq and gets into deadlock, workaround?

    def test_translated_seq(self, _transcript: Transcript, _translation: Translation) -> None:
        """Test that translated_seq is calculated correctly."""
        pass


class TestTranscriptEquality:
    """Test equality comparison."""
    pass
