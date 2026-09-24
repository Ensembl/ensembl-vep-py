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
"""Unit tests for TranscriptVariationAllele class."""

import pytest

from ensembl.vep.common import CoordSystemLevel, Readthrough, Strand, is_resolved, reverse_complement
from ensembl.vep.common.codons import translate_dna
from ensembl.vep.models.genic import Exon, Transcript, Translation
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.overlap.transcript_variation_allele import MAX_READTHROUGH_BASES, _edited_boundary
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.models.variation import VariationAllele

CODING_SEQUENCE = "ATG" + "GGG" * 10 + "TAA"  # (M + 10*G + stop), transcript-sense
REFERENCE_PEPTIDE = "M" + "G" * 10 + "*"


class DummySlice(Slice):
    """Slice with a sequence supplied through the unsafe test hook."""

    def __init__(self, sequence: str) -> None:
        super().__init__(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="chrPin",
            start=1,
            end=len(sequence),
            strand=Strand.FORWARD,
        )
        with self.unsafe() as unsafe_slice:
            setattr(unsafe_slice, "_seq", sequence)


def _write_sequence_at_position(chars: list[str], position: int, sequence: str) -> None:
    """Write a sequence into a character list at a 1-based position."""
    for offset, base in enumerate(sequence):
        chars[position - 1 + offset] = base


def _forward_exons(genomic_slice: Slice | None) -> list[Exon]:
    """Create the two forward-strand exons shared by every forward transcript in this module."""
    return [
        Exon(start=1, end=20, strand=Strand.FORWARD, stable_id="ENST_FWD_e0", slice=genomic_slice),
        Exon(start=30, end=74, strand=Strand.FORWARD, stable_id="ENST_FWD_e1", slice=genomic_slice),
    ]


def _seed_coding_boundaries(transcript: Transcript) -> None:
    """Seed the coding-exon boundaries shared by every transcript in this module."""
    transcript.seed_attrs(
        first_coding_exon_index=0, last_coding_exon_index=1, exon_coding_start=10, exon_coding_end=25
    )


def _forward_transcript(
    *,
    five_prime_utr: str = "N" * 9,
    three_prime_utr: str = "",
    downstream_extra: str = "",
    genome_length: int = 74,
) -> Transcript:
    """Create a forward-strand transcript with two exons and configurable flanking sequence.

    Exon one (genomic 1-20) holds the 5' UTR and the first 11 CDS bases from genomic 10. Exon two
    (genomic 30-74) holds the remaining CDS and the 3' UTR from genomic 55. Downstream sequence is
    placed beyond the transcript's genomic end.
    """
    chars = list("N" * max(genome_length, 74 + len(downstream_extra)))
    _write_sequence_at_position(chars, 1, five_prime_utr)
    _write_sequence_at_position(chars, 10, CODING_SEQUENCE[:11])
    _write_sequence_at_position(chars, 30, CODING_SEQUENCE[11:])
    _write_sequence_at_position(chars, 55, three_prime_utr)
    _write_sequence_at_position(chars, 75, downstream_extra)
    genomic_slice = DummySlice("".join(chars))

    transcript = Transcript(start=1, end=74, strand=Strand.FORWARD, stable_id="ENST_FWD", slice=genomic_slice)
    transcript.spliced_exons = _forward_exons(genomic_slice)
    _seed_coding_boundaries(transcript)
    return transcript


def _reverse_transcript() -> Transcript:
    """Create a reverse-strand transcript mirroring the forward layout.

    Exon one (genomic 55-74) holds the 5' UTR and the first 11 CDS bases, so the CDS starts at genomic
    65. Exon two (genomic 1-45) holds the remaining CDS.
    """
    chars = list("N" * 74)
    _write_sequence_at_position(chars, 55, reverse_complement("N" * 9 + CODING_SEQUENCE[:11]))
    exon_two_sequence = CODING_SEQUENCE[11:]
    _write_sequence_at_position(chars, 45 - len(exon_two_sequence) + 1, reverse_complement(exon_two_sequence))
    genomic_slice = DummySlice("".join(chars))

    transcript = Transcript(start=1, end=74, strand=Strand.REVERSE, stable_id="ENST_REV", slice=genomic_slice)
    transcript.spliced_exons = [
        Exon(start=55, end=74, strand=Strand.REVERSE, stable_id="ENST_REV_e0", slice=genomic_slice),
        Exon(start=1, end=45, strand=Strand.REVERSE, stable_id="ENST_REV_e1", slice=genomic_slice),
    ]
    _seed_coding_boundaries(transcript)
    return transcript


def _reassign_forward_coding_exons(transcript: Transcript, coding_sequence: str) -> None:
    """Replace a forward transcript's exons with fresh ones carrying ``coding_sequence``.

    Fresh Exon and Slice objects sidestep each Exon's own ``_seq`` cache, which mutating the existing
    slice would not invalidate.
    """
    chars = list("N" * 74)
    _write_sequence_at_position(chars, 10, coding_sequence[:11])
    _write_sequence_at_position(chars, 30, coding_sequence[11:])
    transcript.spliced_exons = _forward_exons(DummySlice("".join(chars)))
    # Reassigning spliced_exons invalidates the seeded coding boundaries too.
    _seed_coding_boundaries(transcript)


def _variation_allele(
    transcript: Transcript, start: int, end: int, reference: str, alternate: str
) -> VariationAllele:
    """Create a variation allele on the transcript's strand and slice."""
    return VariationAllele(
        start=start,
        end=end,
        strand=transcript.strand,
        slice=transcript.slice,
        reference_allele=reference,
        alternate_allele=alternate,
    )


def _transcript_variation_allele(
    transcript: Transcript,
    variation_allele: VariationAllele,
    *,
    readthrough: Readthrough = Readthrough.NONE,
) -> TranscriptVariationAllele:
    """Create a transcript overlap for an existing variation allele."""
    return TranscriptVariationAllele(
        feature=transcript,
        overlapping_feature=variation_allele,
        readthrough=readthrough,
    )


@pytest.fixture
def transcript() -> Transcript:
    """Forward-strand coding transcript fixture for testing."""
    return _forward_transcript()


@pytest.fixture
def reverse_transcript() -> Transcript:
    """Reverse-strand coding transcript fixture for testing."""
    return _reverse_transcript()


@pytest.fixture
def transcript_with_utr() -> Transcript:
    """Forward-strand transcript fixture with a 20-base 3' UTR."""
    return _forward_transcript(three_prime_utr="N" * 20)


@pytest.fixture
def readthrough_transcript() -> Transcript:
    """Forward-strand transcript fixture with an in-frame stop codon in its 3' UTR."""
    return _forward_transcript(three_prime_utr="TAA" + "N" * 17)


@pytest.fixture
def long_readthrough_transcript() -> Transcript:
    """Forward-strand transcript fixture on a long slice without a downstream stop codon."""
    return _forward_transcript(three_prime_utr="N" * 20, genome_length=4000)


@pytest.fixture
def transcript_without_exons() -> Transcript:
    """Transcript fixture with no spliced exons."""
    return Transcript(start=1, end=10, strand=Strand.FORWARD, stable_id="ENST_EMPTY")


@pytest.fixture
def variation_allele(transcript: Transcript) -> VariationAllele:
    """Variation allele fixture on the coding transcript."""
    return _variation_allele(transcript, 50, 150, "A", "T")


@pytest.fixture
def transcript_variation_allele(
    transcript: Transcript, variation_allele: VariationAllele
) -> TranscriptVariationAllele:
    """TranscriptVariationAllele fixture for testing."""
    return _transcript_variation_allele(transcript, variation_allele)


class TestTranscriptVariationAlleleInitialization:
    """Test TranscriptVariationAllele initialization."""

    def test_initialization_without_keywords(
        self, transcript: Transcript, variation_allele: VariationAllele
    ) -> None:
        """Test that initializing without keywords raises error."""
        with pytest.raises(TypeError, match="takes 1 positional argument"):
            _ = TranscriptVariationAllele(transcript, variation_allele)  # type: ignore[misc]


class TestTranscriptVariationAlleleProperties:
    """Test TranscriptVariationAllele properties."""

    def test_core_properties(
        self,
        transcript_variation_allele: TranscriptVariationAllele,
        transcript: Transcript,
        variation_allele: VariationAllele,
    ) -> None:
        """Test that feature and overlapping_feature are set correctly."""
        assert transcript_variation_allele.feature == transcript
        assert transcript_variation_allele.overlapping_feature == variation_allele

    def test_default_properties(self, transcript_variation_allele: TranscriptVariationAllele) -> None:
        """Test that default properties are set correctly."""
        assert transcript_variation_allele.is_reference is False
        assert transcript_variation_allele.readthrough == Readthrough.NONE


class TestCodonAlleles:
    """Test codon and peptide window alleles."""

    @pytest.mark.parametrize(
        (
            "start",
            "end",
            "reference",
            "alternate",
            "reference_codon",
            "alternate_codon",
            "reference_peptide",
            "alternate_peptide",
        ),
        [
            pytest.param(14, 14, "G", "A", "GGG", "GAG", "G", "E", id="missense_snv"),
            pytest.param(11, 11, "T", "A", "ATG", "AAG", "M", "K", id="start_codon_snv"),
            pytest.param(14, 15, "GG", "AC", "GGG", "GAC", "G", "D", id="mnv"),
            pytest.param(30, 30, "G", "A", "GGG", "GGA", "G", "G", id="snv_in_codon_split_across_exons"),
            pytest.param(15, 17, "GGG", "-", "GGGGGG", "GGG", "GG", "G", id="inframe_deletion"),
            pytest.param(15, 14, "-", "GGG", "GGG", "GGGGGG", "G", "GG", id="inframe_insertion"),
            pytest.param(15, 15, "G", "-", "GGG", "GG", "G", "X", id="frameshift_deletion"),
            pytest.param(15, 14, "-", "AA", "GGG", "GGAAG", "G", "GX", id="frameshift_insertion"),
            pytest.param(15, 17, "GGG", "A", "GGGGGG", "GGAG", "GG", "GX", id="delins"),
            pytest.param(8, 11, "NNAT", "-", "ATG", "G", "M", "X", id="deletion_from_utr_into_start_codon"),
            pytest.param(19, 22, "GGNN", "-", "GGG", "G", "G", "X", id="deletion_from_exon_into_intron"),
        ],
    )
    def test_codon_alleles(
        self,
        transcript: Transcript,
        start: int,
        end: int,
        reference: str,
        alternate: str,
        reference_codon: str,
        alternate_codon: str,
        reference_peptide: str,
        alternate_peptide: str,
    ) -> None:
        """Test that codon and peptide alleles cover every codon the variant touches."""
        variation_allele = _variation_allele(transcript, start, end, reference, alternate)
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert transcript_variation_allele.reference_codon == reference_codon
        assert transcript_variation_allele.alternate_codon == alternate_codon
        assert transcript_variation_allele.reference_peptide == reference_peptide
        assert transcript_variation_allele.alternate_peptide == alternate_peptide

    @pytest.mark.parametrize(
        (
            "start",
            "end",
            "reference",
            "alternate",
            "reference_codon",
            "alternate_codon",
            "reference_peptide",
            "alternate_peptide",
        ),
        [
            pytest.param(65, 65, "T", "A", "ATG", "TTG", "M", "L", id="start_codon_first_base_snv"),
            pytest.param(63, 63, "C", "T", "ATG", "ATA", "M", "I", id="start_codon_last_base_snv"),
            pytest.param(60, 62, "CCC", "-", "GGG", "", "G", "", id="inframe_deletion"),
            pytest.param(62, 61, "-", "AAA", "GGG", "GTTTGG", "G", "VW", id="inframe_insertion"),
        ],
    )
    def test_codon_alleles_reverse_strand(
        self,
        reverse_transcript: Transcript,
        start: int,
        end: int,
        reference: str,
        alternate: str,
        reference_codon: str,
        alternate_codon: str,
        reference_peptide: str,
        alternate_peptide: str,
    ) -> None:
        """Test that reverse-strand codon alleles are transcript-sense."""
        variation_allele = _variation_allele(reverse_transcript, start, end, reference, alternate)
        transcript_variation_allele = _transcript_variation_allele(reverse_transcript, variation_allele)

        assert transcript_variation_allele.reference_codon == reference_codon
        assert transcript_variation_allele.alternate_codon == alternate_codon
        assert transcript_variation_allele.reference_peptide == reference_peptide
        assert transcript_variation_allele.alternate_peptide == alternate_peptide

    def test_codon_alleles_include_the_whole_stop_codon(self, transcript: Transcript) -> None:
        """Test that a stop-codon SNV's codon window covers all three stop-codon bases."""
        variation_allele = _variation_allele(transcript, 52, 52, "T", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert transcript_variation_allele.reference_codon == "TAA"
        assert transcript_variation_allele.alternate_codon == "AAA"
        assert transcript_variation_allele.reference_peptide == "*"
        assert transcript_variation_allele.alternate_peptide == "K"

    @pytest.mark.parametrize(
        ("start", "end", "reference", "alternate"),
        [
            pytest.param(5, 5, "N", "A", id="five_prime_utr_snv"),
            pytest.param(24, 24, "N", "A", id="intronic_snv"),
            pytest.param(60, 60, "N", "A", id="three_prime_utr_snv"),
        ],
    )
    def test_codon_alleles_are_none_outside_cds(
        self, transcript: Transcript, start: int, end: int, reference: str, alternate: str
    ) -> None:
        """Test that a variant outside the CDS has no codon or peptide alleles."""
        variation_allele = _variation_allele(transcript, start, end, reference, alternate)
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert transcript_variation_allele.reference_codon is None
        assert transcript_variation_allele.alternate_codon is None
        assert transcript_variation_allele.reference_peptide is None
        assert transcript_variation_allele.alternate_peptide is None

    def test_codon_alleles_are_none_without_spliced_exons(self, transcript_without_exons: Transcript) -> None:
        """Test that a transcript without exons has no codon or peptide alleles."""
        variation_allele = _variation_allele(transcript_without_exons, 5, 5, "N", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript_without_exons, variation_allele)

        assert transcript_variation_allele.reference_codon is None
        assert transcript_variation_allele.alternate_peptide is None


class TestCdsLengthChange:
    """Test the net coding length change of a variant."""

    @pytest.mark.parametrize(
        ("start", "end", "reference", "alternate", "cds_length_change"),
        [
            pytest.param(15, 15, "G", "A", 0, id="snv"),
            pytest.param(14, 15, "GG", "AC", 0, id="mnv"),
            pytest.param(15, 17, "GGG", "-", -3, id="inframe_deletion"),
            pytest.param(15, 14, "-", "GGG", 3, id="inframe_insertion"),
            pytest.param(15, 15, "G", "-", -1, id="frameshift_deletion"),
            pytest.param(15, 14, "-", "AA", 2, id="frameshift_insertion"),
            pytest.param(15, 17, "GGG", "A", -2, id="delins"),
            pytest.param(8, 12, "NNATG", "-", -3, id="deletion_from_utr_into_start_codon"),
            pytest.param(19, 22, "GGNN", "-", -2, id="deletion_from_exon_into_intron"),
            pytest.param(20, 30, "G" + "N" * 9 + "G", "-", -2, id="deletion_across_whole_intron"),
            pytest.param(53, 56, "AANN", "-", -2, id="deletion_from_stop_codon_into_utr"),
        ],
    )
    def test_cds_length_change(
        self,
        transcript: Transcript,
        start: int,
        end: int,
        reference: str,
        alternate: str,
        cds_length_change: int,
    ) -> None:
        """Test that the length change counts only the coding bases the variant edits."""
        variation_allele = _variation_allele(transcript, start, end, reference, alternate)
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert transcript_variation_allele.cds_length_change == cds_length_change

    @pytest.mark.parametrize(
        ("start", "end", "reference", "alternate", "cds_length_change"),
        [
            pytest.param(63, 63, "C", "T", 0, id="snv"),
            pytest.param(60, 62, "CCC", "-", -3, id="inframe_deletion"),
            pytest.param(62, 61, "-", "AAA", 3, id="inframe_insertion"),
            pytest.param(60, 60, "C", "-", -1, id="frameshift_deletion"),
        ],
    )
    def test_cds_length_change_reverse_strand(
        self,
        reverse_transcript: Transcript,
        start: int,
        end: int,
        reference: str,
        alternate: str,
        cds_length_change: int,
    ) -> None:
        """Test that the length change is strand-independent."""
        variation_allele = _variation_allele(reverse_transcript, start, end, reference, alternate)
        transcript_variation_allele = _transcript_variation_allele(reverse_transcript, variation_allele)

        assert transcript_variation_allele.cds_length_change == cds_length_change

    @pytest.mark.parametrize(
        ("start", "end", "reference", "alternate"),
        [
            pytest.param(3, 3, "N", "A", id="five_prime_utr_snv"),
            pytest.param(24, 24, "N", "A", id="intronic_snv"),
            pytest.param(60, 62, "NNN", "-", id="three_prime_utr_deletion"),
        ],
    )
    def test_cds_length_change_is_none_outside_cds(
        self, transcript: Transcript, start: int, end: int, reference: str, alternate: str
    ) -> None:
        """Test that a variant outside the CDS has no coding length change."""
        variation_allele = _variation_allele(transcript, start, end, reference, alternate)
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert transcript_variation_allele.cds_length_change is None

    def test_cds_length_change_is_none_without_spliced_exons(
        self, transcript_without_exons: Transcript
    ) -> None:
        """Test that a transcript without exons has no coding length change."""
        variation_allele = _variation_allele(transcript_without_exons, 5, 5, "N", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript_without_exons, variation_allele)

        assert transcript_variation_allele.cds_length_change is None

    def test_cds_length_change_excludes_cdna_only_bases(self, transcript: Transcript) -> None:
        """Test that a deletion across the CDS start counts coding bases, not cDNA bases."""
        variation_allele = _variation_allele(transcript, 8, 12, "NNATG", "-")
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)
        cdna_edit = transcript_variation_allele._cdna_edit  # pylint: disable=protected-access

        assert transcript_variation_allele.cds_length_change == -3
        assert cdna_edit is not None
        assert cdna_edit.length_change == -5


class TestSequenceEditGeometry:
    """Test CDS and cDNA edit geometry."""

    def test_cds_and_cdna_edit_footprints_differ_at_cds_start(self, transcript: Transcript) -> None:
        """Test that a variant across the CDS start has distinct CDS and cDNA footprints."""
        # The 5' UTR is genomic 1-9 and the CDS starts at 10, so this removes 2 UTR and 3 CDS bases.
        variation_allele = _variation_allele(transcript, 8, 12, "NNATG", "-")
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)
        cds_edit = transcript_variation_allele._cds_edit  # pylint: disable=protected-access
        cdna_edit = transcript_variation_allele._cdna_edit  # pylint: disable=protected-access

        assert cds_edit is not None
        assert cds_edit.footprint == (0, 2)
        assert cdna_edit is not None
        assert cdna_edit.footprint == (7, 11)

    @pytest.mark.parametrize(
        ("original_boundary", "edit_start", "edit_end", "length_change", "expected"),
        [
            pytest.param(5, 10, 12, -2, 5, id="edit_after_boundary"),
            pytest.param(15, 10, 12, -2, 13, id="deletion_before_boundary"),
            pytest.param(15, 10, 10, 3, 18, id="insertion_before_boundary"),
            pytest.param(10, 10, 10, 3, 13, id="insertion_at_boundary"),
            pytest.param(11, 10, 12, -2, 10, id="edit_over_boundary"),
        ],
    )
    def test_edited_boundary(
        self, original_boundary: int, edit_start: int, edit_end: int, length_change: int, expected: int
    ) -> None:
        """Test that a boundary is kept, shifted, or anchored to the edit start."""
        assert _edited_boundary(original_boundary, edit_start, edit_end, length_change) == expected


class TestFullCdsAlleles:
    """Test reference and alternate CDS and peptide properties."""

    def test_reference_cds_and_peptide(self, transcript: Transcript) -> None:
        """Test that a forward transcript exposes its unedited CDS and peptide."""
        variation_allele = _variation_allele(transcript, 5, 5, "N", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert transcript_variation_allele.reference_cds == CODING_SEQUENCE
        assert transcript_variation_allele.reference_cds_peptide == REFERENCE_PEPTIDE

    @pytest.mark.parametrize(
        ("start", "end", "reference", "alternate"),
        [
            pytest.param(5, 5, "N", "A", id="five_prime_utr_snv"),
            pytest.param(3, 6, "NNNN", "-", id="five_prime_utr_deletion"),
            pytest.param(10, 9, "-", "GGG", id="insertion_before_start_codon"),
            pytest.param(24, 24, "N", "A", id="intronic_snv"),
            pytest.param(23, 26, "NNNN", "-", id="intronic_deletion"),
        ],
    )
    def test_alternate_cds_is_unchanged_by_non_coding_variant(
        self, transcript: Transcript, start: int, end: int, reference: str, alternate: str
    ) -> None:
        """Test that a variant outside the CDS leaves the alternate CDS and peptide unchanged."""
        variation_allele = _variation_allele(transcript, start, end, reference, alternate)
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert transcript_variation_allele.alternate_cds == CODING_SEQUENCE
        assert transcript_variation_allele.alternate_cds_peptide == REFERENCE_PEPTIDE

    def test_alternate_cds_is_unchanged_by_variant_outside_transcript(self) -> None:
        """Test that a variant beyond the transcript leaves the alternate CDS unchanged."""
        transcript = _forward_transcript(genome_length=120)
        variation_allele = _variation_allele(transcript, 100, 100, "N", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert transcript_variation_allele.alternate_cds == CODING_SEQUENCE

    @pytest.mark.parametrize(
        ("start", "end", "reference", "alternate"),
        [
            pytest.param(8, 11, "NNAT", "-", id="deletion_from_utr_into_start_codon"),
            pytest.param(10, 12, "ATG", "-", id="whole_start_codon_deletion"),
        ],
    )
    def test_alternate_cds_peptide_loses_start_methionine(
        self, transcript: Transcript, start: int, end: int, reference: str, alternate: str
    ) -> None:
        """Test that deleting start-codon bases removes the start methionine."""
        variation_allele = _variation_allele(transcript, start, end, reference, alternate)
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert transcript_variation_allele.alternate_cds_peptide is not None
        assert not transcript_variation_allele.alternate_cds_peptide.startswith("M")

    def test_alternate_cds_peptide_keeps_start_methionine(self, transcript: Transcript) -> None:
        """Test that an edit preserving the start codon retains the start methionine."""
        variation_allele = _variation_allele(transcript, 11, 10, "A", "AT")
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert transcript_variation_allele.alternate_cds_peptide is not None
        assert transcript_variation_allele.alternate_cds_peptide.startswith("M")

    def test_alternate_cds_anchors_at_start_codon_delins(self, transcript: Transcript) -> None:
        """Test that a delins across the start codon anchors the alternate CDS at the edit start."""
        assert transcript.spliced_seq is not None
        # Replace 2 bases of 5' UTR and the first ATG base; keep the untouched CDS remainder.
        expected_alternate_cds = "G" * 10 + transcript.spliced_seq[10:45]
        variation_allele = _variation_allele(transcript, 8, 10, "NNA", "G" * 10)
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert transcript_variation_allele.alternate_cds == expected_alternate_cds
        assert transcript_variation_allele.alternate_cds_peptide is not None
        assert not transcript_variation_allele.alternate_cds_peptide.startswith("M")

    def test_alternate_cds_peptide_loses_stop_codon(self, transcript_with_utr: Transcript) -> None:
        """Test that a stop-codon SNV stays bounded by the annotated CDS end."""
        variation_allele = _variation_allele(transcript_with_utr, 52, 52, "T", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript_with_utr, variation_allele)

        assert transcript_variation_allele.reference_cds_peptide == REFERENCE_PEPTIDE
        assert transcript_variation_allele.alternate_cds_peptide == "M" + "G" * 10 + "K"

    @pytest.mark.parametrize(
        ("start", "end", "reference", "alternate"),
        [
            pytest.param(52, 74, "N" * 23, "-", id="stop_codon_and_utr_deletion"),
            pytest.param(52, 55, "TAAN", "A" * 10, id="stop_codon_delins"),
        ],
    )
    def test_alternate_cds_truncates_at_stop_codon_edit(
        self, transcript_with_utr: Transcript, start: int, end: int, reference: str, alternate: str
    ) -> None:
        """Test that an edit across the stop codon excludes sequence beyond the CDS end."""
        assert transcript_with_utr.spliced_seq is not None
        expected_alternate_cds = transcript_with_utr.spliced_seq[9:42]  # ATG and the poly-G run
        variation_allele = _variation_allele(transcript_with_utr, start, end, reference, alternate)
        transcript_variation_allele = _transcript_variation_allele(transcript_with_utr, variation_allele)

        assert transcript_variation_allele.alternate_cds == expected_alternate_cds
        assert transcript_variation_allele.alternate_cds_peptide == "M" + "G" * 10

    def test_alternate_cds_is_empty_after_whole_cds_deletion(self, transcript_with_utr: Transcript) -> None:
        """Test that deleting the complete CDS produces an empty alternate CDS."""
        variation_allele = _variation_allele(transcript_with_utr, 10, 54, CODING_SEQUENCE, "-")
        transcript_variation_allele = _transcript_variation_allele(transcript_with_utr, variation_allele)

        assert transcript_variation_allele.alternate_cds == ""
        assert transcript_variation_allele.alternate_cds_peptide == ""

    def test_alternate_cds_peptide_gains_residue_after_inframe_insertion(
        self, transcript_with_utr: Transcript
    ) -> None:
        """Test that an inframe CDS insertion adds a residue to the alternate peptide."""
        variation_allele = _variation_allele(transcript_with_utr, 15, 14, "-", "GGG")
        transcript_variation_allele = _transcript_variation_allele(transcript_with_utr, variation_allele)

        assert transcript_variation_allele.alternate_cds_peptide == "M" + "G" * 11 + "*"

    def test_alternate_cds_peptide_is_unchanged_by_reverse_strand_utr_snv(
        self, reverse_transcript: Transcript
    ) -> None:
        """Test that a reverse-strand 5' UTR SNV leaves the peptide unchanged."""
        variation_allele = _variation_allele(reverse_transcript, 72, 72, "N", "A")
        transcript_variation_allele = _transcript_variation_allele(reverse_transcript, variation_allele)

        assert transcript_variation_allele.alternate_cds_peptide == REFERENCE_PEPTIDE

    def test_alternate_cds_peptide_loses_start_methionine_on_reverse_strand(
        self, reverse_transcript: Transcript
    ) -> None:
        """Test that a reverse-strand deletion across the start codon removes the start methionine."""
        variation_allele = _variation_allele(reverse_transcript, 64, 67, "NNNN", "-")
        transcript_variation_allele = _transcript_variation_allele(reverse_transcript, variation_allele)

        assert transcript_variation_allele.alternate_cds_peptide is not None
        assert not transcript_variation_allele.alternate_cds_peptide.startswith("M")

    def test_cds_alleles_are_none_without_spliced_exons(self, transcript_without_exons: Transcript) -> None:
        """Test that a transcript without exons has no reference or alternate CDS."""
        variation_allele = _variation_allele(transcript_without_exons, 5, 5, "N", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript_without_exons, variation_allele)

        assert transcript_variation_allele.reference_cds is None
        assert transcript_variation_allele.alternate_cds is None
        assert transcript_variation_allele.reference_cds_peptide is None
        assert transcript_variation_allele.alternate_cds_peptide is None


class TestFullCdsAllelesReadthrough:
    """Test CDS and peptide properties across readthrough levels."""

    @pytest.mark.parametrize(
        ("start", "end", "reference", "alternate"),
        [
            pytest.param(15, 15, "G", "A", id="cds_snv"),
            pytest.param(15, 17, "GGG", "-", id="inframe_cds_deletion"),
            pytest.param(5, 5, "N", "A", id="five_prime_utr_snv"),
        ],
    )
    def test_alternate_cds_is_readthrough_invariant_when_stop_codon_survives(
        self, readthrough_transcript: Transcript, start: int, end: int, reference: str, alternate: str
    ) -> None:
        """Test that readthrough levels have no effect when the stop codon survives."""
        variation_allele = _variation_allele(readthrough_transcript, start, end, reference, alternate)
        alternate_cds_by_level = {
            _transcript_variation_allele(
                readthrough_transcript, variation_allele, readthrough=level
            ).alternate_cds
            for level in Readthrough
        }

        assert len(alternate_cds_by_level) == 1

    def test_utr_readthrough_finds_stop_codon_in_three_prime_utr(
        self, readthrough_transcript: Transcript
    ) -> None:
        """Test that UTR readthrough extends the alternate CDS to the next in-frame stop codon."""
        variation_allele = _variation_allele(readthrough_transcript, 52, 52, "T", "A")
        default = _transcript_variation_allele(readthrough_transcript, variation_allele)
        readthrough = _transcript_variation_allele(
            readthrough_transcript, variation_allele, readthrough=Readthrough.UTR
        )

        assert default.alternate_cds_peptide == "M" + "G" * 10 + "K"
        assert readthrough.alternate_cds_peptide == "M" + "G" * 10 + "K*"
        assert default.alternate_cds is not None
        assert readthrough.alternate_cds == default.alternate_cds + "TAA"

    def test_utr_readthrough_stops_at_transcript_end(self) -> None:
        """Test that UTR readthrough does not read past the spliced transcript."""
        transcript = _forward_transcript(three_prime_utr="N" * 20, downstream_extra="N" + "TAA" + "AAAA")
        variation_allele = _variation_allele(transcript, 52, 52, "T", "A")
        readthrough = _transcript_variation_allele(transcript, variation_allele, readthrough=Readthrough.UTR)

        assert readthrough.alternate_cds is not None
        assert len(readthrough.alternate_cds) == 56  # anchor 9 to the end of the 65-base cDNA
        assert readthrough.alternate_cds_peptide is not None
        assert "*" not in readthrough.alternate_cds_peptide

    def test_beyond_transcript_readthrough_reads_downstream_flanking_sequence(self) -> None:
        """Test that beyond-transcript readthrough finds a stop codon in downstream sequence."""
        transcript = _forward_transcript(three_prime_utr="N" * 20, downstream_extra="N" + "TAA" + "AAAA")
        variation_allele = _variation_allele(transcript, 52, 52, "T", "A")
        beyond = _transcript_variation_allele(
            transcript, variation_allele, readthrough=Readthrough.BEYOND_TRANSCRIPT
        )

        assert beyond.alternate_cds is not None
        assert len(beyond.alternate_cds) == 60  # 56 in-transcript + "N" completing the codon + TAA
        assert beyond.alternate_cds_peptide is not None
        assert beyond.alternate_cds_peptide.endswith("*")

    def test_beyond_transcript_readthrough_when_downstream_flanking_sequence_not_available(self) -> None:
        """Test that beyond-transcript readthrough finds a stop codon in downstream sequence."""
        transcript = _forward_transcript()
        variation_allele = _variation_allele(transcript, 52, 52, "T", "A")
        beyond = _transcript_variation_allele(
            transcript, variation_allele, readthrough=Readthrough.BEYOND_TRANSCRIPT
        )

        assert beyond.alternate_cds is not None
        assert len(beyond.alternate_cds) == 56  # 56 length of the whole transcript
        assert beyond.alternate_cds_peptide is not None
        assert beyond.alternate_cds_peptide.endswith("X")

    def test_beyond_transcript_readthrough_is_bounded(self, long_readthrough_transcript: Transcript) -> None:
        """Test that beyond-transcript readthrough is limited to the maximum distance."""
        variation_allele = _variation_allele(long_readthrough_transcript, 52, 52, "T", "A")
        beyond = _transcript_variation_allele(
            long_readthrough_transcript, variation_allele, readthrough=Readthrough.BEYOND_TRANSCRIPT
        )

        assert beyond.alternate_cds is not None
        assert len(beyond.alternate_cds) <= 56 + MAX_READTHROUGH_BASES
        assert beyond.alternate_cds_peptide is not None
        assert "*" not in beyond.alternate_cds_peptide

    def test_reference_cds_is_readthrough_invariant(self, readthrough_transcript: Transcript) -> None:
        """Test that the reference CDS is unchanged by the readthrough level."""
        variation_allele = _variation_allele(readthrough_transcript, 52, 52, "T", "A")
        reference_cds_by_level = {
            _transcript_variation_allele(
                readthrough_transcript, variation_allele, readthrough=level
            ).reference_cds
            for level in Readthrough
        }

        assert reference_cds_by_level == {CODING_SEQUENCE}

    def test_cds_peptides_translate_cds_at_non_default_readthrough(
        self, readthrough_transcript: Transcript
    ) -> None:
        """Test that nucleotide and peptide accessors use the same readthrough level."""
        variation_allele = _variation_allele(readthrough_transcript, 52, 52, "T", "A")
        transcript_variation_allele = _transcript_variation_allele(
            readthrough_transcript, variation_allele, readthrough=Readthrough.UTR
        )

        assert transcript_variation_allele.alternate_cds is not None
        assert transcript_variation_allele.alternate_cds_peptide == translate_dna(
            transcript_variation_allele.alternate_cds
        )
        assert transcript_variation_allele.reference_cds is not None
        assert transcript_variation_allele.reference_cds_peptide == translate_dna(
            transcript_variation_allele.reference_cds
        )

    def test_readthrough_reassignment_invalidates_cached_alleles(
        self, readthrough_transcript: Transcript
    ) -> None:
        """Test that reassigning readthrough invalidates cached allele sequences."""
        variation_allele = _variation_allele(readthrough_transcript, 52, 52, "T", "A")
        transcript_variation_allele = _transcript_variation_allele(readthrough_transcript, variation_allele)

        assert transcript_variation_allele.alternate_cds_peptide == "M" + "G" * 10 + "K"
        transcript_variation_allele.readthrough = Readthrough.UTR
        assert transcript_variation_allele.alternate_cds_peptide == "M" + "G" * 10 + "K*"

    def test_downstream_flanking_sequence_is_reverse_complemented_on_reverse_strand(self) -> None:
        """Test that reverse-strand downstream sequence is reverse-complemented."""
        genomic_slice = DummySlice("ACGTACGTACGTACGTACG" + "N" * 30)
        transcript = Transcript(
            start=20, end=30, strand=Strand.REVERSE, stable_id="ENST_REV", slice=genomic_slice
        )
        transcript.spliced_exons = [
            Exon(start=20, end=30, strand=Strand.REVERSE, stable_id="ENST_REV_e0", slice=genomic_slice)
        ]
        variation_allele = _variation_allele(transcript, 20, 20, "T", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)
        # pylint: disable-next=protected-access
        flank = transcript_variation_allele._downstream_flanking_sequence(10)

        assert genomic_slice.seq is not None
        assert flank == reverse_complement(genomic_slice.seq[9:19])

    def test_downstream_flanking_sequence_is_empty_when_length_is_zero(
        self, transcript_with_utr: Transcript
    ) -> None:
        """Test that missing downstream slice sequence returns an empty string."""
        variation_allele = _variation_allele(transcript_with_utr, 52, 52, "T", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript_with_utr, variation_allele)
        # pylint: disable-next=protected-access
        flank = transcript_variation_allele._downstream_flanking_sequence(0)
        assert flank == ""

    def test_downstream_flanking_sequence_is_empty_when_slice_is_none(self) -> None:
        """Test that a transcript with no slice returns empty string"""
        # create a transcript with no slice
        transcript = Transcript(start=1, end=74, strand=Strand.FORWARD, stable_id="ENST_FWD")
        transcript.spliced_exons = _forward_exons(None)
        _seed_coding_boundaries(transcript)
        variation_allele = _variation_allele(transcript, 52, 52, "T", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)
        # pylint: disable-next=protected-access
        flank = transcript_variation_allele._downstream_flanking_sequence(50)
        assert flank == ""

    def test_downstream_flanking_sequence_is_empty_when_slice_seq_is_none(self) -> None:
        """Test that a transcript with no slice sequence returns empty string"""
        # create a slice with no sequence source
        slice_no_seq = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="chrPin",
            start=1,
            end=74,
            strand=Strand.FORWARD,
            sequence_source=None,
        )
        transcript = Transcript(
            start=1, end=74, strand=Strand.FORWARD, stable_id="ENST_FWD", slice=slice_no_seq
        )
        transcript.spliced_exons = _forward_exons(slice_no_seq)
        _seed_coding_boundaries(transcript)
        variation_allele = _variation_allele(transcript, 52, 52, "T", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)
        # pylint: disable-next=protected-access
        flank = transcript_variation_allele._downstream_flanking_sequence(50)
        assert flank == ""

    def test_downstream_flanking_sequence_is_empty_past_slice(self, transcript_with_utr: Transcript) -> None:
        """Test that missing downstream slice sequence returns an empty string."""
        variation_allele = _variation_allele(transcript_with_utr, 52, 52, "T", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript_with_utr, variation_allele)
        # pylint: disable-next=protected-access
        flank = transcript_variation_allele._downstream_flanking_sequence(50)

        assert flank == ""

    def test_downstream_flanking_sequence_reverse_transcript_is_empty_past_slice(
        self, reverse_transcript: Transcript
    ) -> None:
        """Test that missing downstream slice sequence returns an empty string with
        a transcript on the reverse strand."""
        variation_allele = _variation_allele(reverse_transcript, 52, 52, "T", "A")
        transcript_variation_allele = _transcript_variation_allele(reverse_transcript, variation_allele)
        # pylint: disable-next=protected-access
        flank = transcript_variation_allele._downstream_flanking_sequence(50)

        assert flank == ""


class TestReferenceCdsPeptideReuse:
    """Test reuse of the transcript's own sequences for the reference CDS and peptide."""

    def test_reference_cds_peptide_reuses_translated_seq_before_translation_resolves(
        self, transcript: Transcript
    ) -> None:
        """Test that the reference peptide is the transcript's translated sequence."""
        assert not is_resolved(transcript, "_translation")
        variation_allele = _variation_allele(transcript, 5, 5, "N", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert (
            transcript_variation_allele.reference_cds_peptide
            == transcript.translated_seq
            == REFERENCE_PEPTIDE
        )

    def test_reference_cds_peptide_reuses_translated_seq_for_seeded_translation_without_seq(
        self, transcript: Transcript
    ) -> None:
        """Test that a seeded translation without its own sequence, as GFF/GTF builds, is reused."""
        derived_translation = transcript.translation
        assert derived_translation is not None
        seeded_translation = Translation(
            start=derived_translation.start,
            end=derived_translation.end,
            stable_id="ENSP_SEEDED",
            exon_coding_start=derived_translation.exon_coding_start,
            exon_coding_end=derived_translation.exon_coding_end,
        )
        transcript.seed_attrs(translation=seeded_translation)
        assert transcript.translation is seeded_translation
        assert transcript.translation.seq is None

        variation_allele = _variation_allele(transcript, 5, 5, "N", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert (
            transcript_variation_allele.reference_cds_peptide
            == transcript.translated_seq
            == REFERENCE_PEPTIDE
        )

    def test_reference_cds_peptide_retranslates_seeded_curated_peptide(self, transcript: Transcript) -> None:
        """Test that a seeded curated peptide is not trusted as the reference peptide."""
        curated_translation = Translation(start=1, end=36, stable_id="ENSP_CURATED")
        curated_translation.seq = "CURATED_PEPTIDE_STRING"
        transcript.seed_attrs(translation=curated_translation)
        assert transcript.translated_seq == "CURATED_PEPTIDE_STRING"

        variation_allele = _variation_allele(transcript, 5, 5, "N", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert transcript_variation_allele.reference_cds_peptide == REFERENCE_PEPTIDE

    def test_reference_cds_alleles_ignore_cds_seq_of_wrong_length(self, transcript: Transcript) -> None:
        """Test that a cached CDS sequence disagreeing with ``cds_length`` is not reused."""
        cds = transcript.cds
        assert cds is not None and cds.seq is not None
        cds.seq = cds.seq[:-3]

        variation_allele = _variation_allele(transcript, 5, 5, "N", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert transcript_variation_allele.reference_cds == CODING_SEQUENCE
        assert transcript_variation_allele.reference_cds_peptide == REFERENCE_PEPTIDE

    def test_reference_cds_reuses_cds_seq(self, transcript: Transcript) -> None:
        """Test that the reference CDS is the transcript's own ``cds.seq`` object."""
        cds = transcript.cds
        assert cds is not None and cds.seq is not None
        variation_allele = _variation_allele(transcript, 5, 5, "N", "A")
        transcript_variation_allele = _transcript_variation_allele(transcript, variation_allele)

        assert transcript_variation_allele.reference_cds is cds.seq

    def test_reference_cds_peptide_reuse_is_independent_of_resolution_order(
        self, transcript: Transcript
    ) -> None:
        """Test that a later overlap reuses a translation resolved by an earlier one."""
        first_variation_allele = _variation_allele(transcript, 5, 5, "N", "A")
        first_overlap = _transcript_variation_allele(transcript, first_variation_allele)
        assert first_overlap.reference_cds_peptide == REFERENCE_PEPTIDE
        assert transcript.translation is not None
        assert transcript.translation.derived_from_coding_sequence is True

        second_variation_allele = _variation_allele(transcript, 15, 15, "G", "A")
        second_overlap = _transcript_variation_allele(transcript, second_variation_allele)

        assert second_overlap.reference_cds_peptide is transcript.translated_seq


class TestFullCdsAllelesInvalidation:
    """Test full CDS alleles after the transcript's exons are reassigned."""

    def test_reference_cds_alleles_follow_reassigned_exons(self) -> None:
        """Test that a fresh overlap sees the new CDS when the translation was never resolved."""
        transcript = _forward_transcript()
        variation_allele = _variation_allele(transcript, 5, 5, "N", "A")
        before = _transcript_variation_allele(transcript, variation_allele)
        assert before.reference_cds == CODING_SEQUENCE
        assert not is_resolved(transcript, "_translation")

        new_coding_sequence = "ATG" + "GGG" * 4 + "GCG" + "GGG" * 5 + "TAA"
        _reassign_forward_coding_exons(transcript, new_coding_sequence)
        after = _transcript_variation_allele(transcript, variation_allele)

        assert after.reference_cds == new_coding_sequence
        assert after.reference_cds_peptide == transcript.translated_seq == translate_dna(new_coding_sequence)

    def test_reference_cds_peptide_is_stale_once_translation_resolved(self) -> None:
        """Test that a resolved translation survives an exon reassignment, while the CDS does not."""
        # Known trade-off: _translation is not in Transcript._CODING_DERIVATIONS, so translated_seq
        # itself is already stale. Production never reassigns spliced_exons after build.
        transcript = _forward_transcript()
        variation_allele = _variation_allele(transcript, 5, 5, "N", "A")
        before = _transcript_variation_allele(transcript, variation_allele)
        assert before.reference_cds_peptide == REFERENCE_PEPTIDE
        stale_translated_seq = transcript.translated_seq

        new_coding_sequence = "ATG" + "GGG" * 4 + "GCG" + "GGG" * 5 + "TAA"
        _reassign_forward_coding_exons(transcript, new_coding_sequence)
        after = _transcript_variation_allele(transcript, variation_allele)

        assert transcript.translated_seq == stale_translated_seq
        assert after.reference_cds == new_coding_sequence
        assert after.reference_cds_peptide == stale_translated_seq
