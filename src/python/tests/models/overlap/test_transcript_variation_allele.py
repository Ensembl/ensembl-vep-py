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
from ensembl.vep.models.overlap.transcript_variation_allele import MAX_READTHROUGH_BASES, _edited_boundary
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import VariationAllele
from ensembl.vep.models.genic import Exon, Transcript, Translation
from ensembl.vep.models.sequence import CoordSystem, Slice

CODING_SEQUENCE = "ATG" + "GGG" * 10 + "TAA"  # (M + 10*G + stop), transcript-sense


class DummySlice(Slice):
    """Slice with a sequence supplied through the unsafe test hook."""

    def __init__(self, sequence: str) -> None:
        """Initialize the slice with a cached sequence."""
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


def _build_forward_transcript(
    *,
    five_prime_utr: str = "N" * 9,
    three_prime_utr: str = "",
    downstream_extra: str = "",
    genome_length: int = 74,
) -> Transcript:
    """Build a forward-strand transcript with two exons and configurable flanking sequence.

    Exon one contains the 5' UTR and first 11 CDS bases. Exon two contains the remaining CDS and
    3' UTR. Downstream sequence is placed beyond the transcript's genomic end.
    """
    exon_one_cds = CODING_SEQUENCE[:11]
    exon_two_rest = CODING_SEQUENCE[11:]
    chars = list("N" * max(genome_length, 74 + len(downstream_extra)))
    _write_sequence_at_position(chars, 1, five_prime_utr)
    _write_sequence_at_position(chars, 10, exon_one_cds)
    _write_sequence_at_position(chars, 30, exon_two_rest)
    _write_sequence_at_position(chars, 55, three_prime_utr)
    _write_sequence_at_position(chars, 75, downstream_extra)
    slice = DummySlice("".join(chars))

    transcript = Transcript(start=1, end=74, strand=Strand.FORWARD, stable_id="ENST_FWD", slice=slice)
    transcript.spliced_exons = [
        Exon(start=1, end=20, strand=Strand.FORWARD, stable_id="ENST_FWD_e0", slice=slice),
        Exon(start=30, end=74, strand=Strand.FORWARD, stable_id="ENST_FWD_e1", slice=slice),
    ]
    transcript.seed_attrs(
        first_coding_exon_index=0, last_coding_exon_index=1, exon_coding_start=10, exon_coding_end=25
    )
    return transcript


def _build_reverse_transcript(
    *,
    five_prime_utr: str = "N" * 9,
    three_prime_utr: str = "",
    downstream_extra: str = "",
    genome_length: int = 74,
) -> Transcript:
    """Build a reverse-strand transcript with two exons and configurable flanking sequence.

    Exon one carries the 5' UTR and first 11 CDS bases; exon two carries the remaining CDS and
    3' UTR. Downstream sequence is placed below the transcript's genomic start.
    """
    exon_one_cds = CODING_SEQUENCE[:11]
    exon_two_rest = CODING_SEQUENCE[11:]
    transcript_offset = len(downstream_extra)
    transcript_start = transcript_offset + 1
    transcript_end = transcript_offset + 74
    exon_one_start, exon_one_end = transcript_offset + 55, transcript_end
    exon_two_start, exon_two_end = transcript_start, transcript_offset + 45

    chars = list("N" * max(genome_length, transcript_end))
    _write_sequence_at_position(chars, 1, reverse_complement(downstream_extra))
    _write_sequence_at_position(
        chars,
        exon_one_start,
        reverse_complement(five_prime_utr + exon_one_cds),
    )
    exon_two_sequence = exon_two_rest + three_prime_utr
    _write_sequence_at_position(
        chars,
        exon_two_end - len(exon_two_sequence) + 1,
        reverse_complement(exon_two_sequence),
    )
    slice = DummySlice("".join(chars))

    transcript = Transcript(
        start=transcript_start,
        end=transcript_end,
        strand=Strand.REVERSE,
        stable_id="ENST_REV",
        slice=slice,
    )
    transcript.spliced_exons = [
        Exon(
            start=exon_one_start,
            end=exon_one_end,
            strand=Strand.REVERSE,
            stable_id="ENST_REV_e0",
            slice=slice,
        ),
        Exon(
            start=exon_two_start,
            end=exon_two_end,
            strand=Strand.REVERSE,
            stable_id="ENST_REV_e1",
            slice=slice,
        ),
    ]
    transcript.seed_attrs(
        first_coding_exon_index=0, last_coding_exon_index=1, exon_coding_start=10, exon_coding_end=25
    )
    return transcript


def _create_variation_allele(
    transcript: Transcript,
    start: int,
    end: int,
    reference: str,
    alternate: str,
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


def _create_transcript_variation_allele(
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
    """Create a forward-strand coding transcript fixture."""
    return _build_forward_transcript()


@pytest.fixture
def reverse_coding_transcript() -> Transcript:
    """Create a reverse-strand coding transcript fixture."""
    return _build_reverse_transcript()


@pytest.fixture
def forward_transcript_with_utr() -> Transcript:
    """Create a forward-strand transcript fixture with a 20-base 3' UTR."""
    return _build_forward_transcript(three_prime_utr="N" * 20)


@pytest.fixture
def forward_readthrough_transcript() -> Transcript:
    """Create a transcript fixture with an in-frame stop in its 3' UTR."""
    return _build_forward_transcript(three_prime_utr="TAA" + "N" * 17)


@pytest.fixture
def forward_long_readthrough_transcript() -> Transcript:
    """Create a transcript fixture on a long slice without a downstream stop codon."""
    return _build_forward_transcript(three_prime_utr="N" * 20, genome_length=4000)


@pytest.fixture
def variation_allele(transcript: Transcript) -> VariationAllele:
    """Create a variation allele fixture on the coding transcript."""
    return _create_variation_allele(transcript, 50, 150, "A", "T")


@pytest.fixture
def transcript_variation_allele(
    transcript: Transcript, variation_allele: VariationAllele
) -> TranscriptVariationAllele:
    """Create a transcript variation allele fixture."""
    return _create_transcript_variation_allele(transcript, variation_allele)


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
        """Test that feature and overlapping_feature properties are set correctly."""
        assert transcript_variation_allele.feature == transcript
        assert transcript_variation_allele.overlapping_feature == variation_allele

    def test_default_properties(self, transcript_variation_allele: TranscriptVariationAllele) -> None:
        """Test that default properties are set correctly."""
        assert transcript_variation_allele.is_reference is False

    def test_cdna_locus(self, transcript_variation_allele: TranscriptVariationAllele) -> None:
        """Test cdna_locus property."""


class TestFullCdsAlleles:
    """Test reference/alternate CDS and peptide properties."""

    def test_forward_reference_cds_and_peptide(self, transcript: Transcript) -> None:
        """Test that a forward transcript exposes its unedited CDS and peptide."""
        variation_allele = _create_variation_allele(transcript, 5, 5, "N", "A")
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert tva.reference_cds == CODING_SEQUENCE
        assert tva.reference_cds_peptide == "M" + "G" * 10 + "*"

    @pytest.mark.parametrize(
        ("start", "end", "reference", "alternate", "label"),
        [
            (5, 5, "N", "A", "5' UTR SNV"),
            (3, 6, "NNNN", "-", "5' UTR deletion"),
            (10, 9, "-", "GGG", "insertion immediately before the start codon"),
        ],
    )
    def test_forward_no_op_upstream_edits(
        self,
        transcript: Transcript,
        start: int,
        end: int,
        reference: str,
        alternate: str,
        label: str,
    ) -> None:
        """Test that an edit upstream of the CDS leaves the peptide unchanged."""
        variation_allele = _create_variation_allele(transcript, start, end, reference, alternate)
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert tva.alternate_cds_peptide == tva.reference_cds_peptide, label

    def test_forward_start_codon_destroyed(self, transcript: Transcript) -> None:
        """Test that deleting UTR and start-codon bases removes the start methionine."""
        variation_allele = _create_variation_allele(transcript, 8, 11, "NNAT", "-")
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert not tva.alternate_cds_peptide.startswith("M")

    def test_forward_whole_start_codon_deletion(self, transcript: Transcript) -> None:
        """Test that deleting the complete start codon removes the start methionine."""
        variation_allele = _create_variation_allele(transcript, 10, 12, "ATG", "-")
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert not tva.alternate_cds_peptide.startswith("M")

    def test_forward_start_codon_snv_retained(self, transcript: Transcript) -> None:
        """Test that an insertion preserving the start codon retains the start methionine."""
        variation_allele = _create_variation_allele(transcript, 11, 10, "A", "AT")
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert tva.alternate_cds_peptide.startswith("M")

    def test_forward_start_codon_delins_anchor_reads_through_the_replacement(
        self, transcript: Transcript
    ) -> None:
        """Test that a start-spanning delins anchors the alternate CDS at the edit start."""
        # Replace 2bp of 5' UTR + the ATG's own first base with a 10bp run of G's.
        expected_tail = "G" * 10 + transcript.spliced_seq[10:45]  # replacement, then untouched CDS remainder
        variation_allele = _create_variation_allele(transcript, 8, 10, "NNA", "G" * 10)
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert not tva.alternate_cds_peptide.startswith("M")
        assert tva.alternate_cds == expected_tail

    def test_forward_intron_only_deletion_is_a_no_op(self, transcript: Transcript) -> None:
        """Test that an intron-only deletion leaves the CDS and peptide unchanged."""
        variation_allele = _create_variation_allele(transcript, 23, 26, "NNNN", "-")
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert tva.reference_cds == CODING_SEQUENCE
        assert tva.alternate_cds == CODING_SEQUENCE
        assert tva.alternate_cds_peptide == tva.reference_cds_peptide

    def test_forward_stop_codon_snv_lost(self, forward_transcript_with_utr: Transcript) -> None:
        """Test that a stop-codon SNV remains bounded by the annotated CDS end."""
        transcript = forward_transcript_with_utr
        variation_allele = _create_variation_allele(transcript, 52, 52, "T", "A")
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert tva.reference_cds_peptide == "M" + "G" * 10 + "*"
        assert tva.alternate_cds_peptide == "M" + "G" * 10 + "K"

    def test_forward_stop_plus_3utr_deletion_truncates_at_the_deletion(
        self, forward_transcript_with_utr: Transcript
    ) -> None:
        """Test that deleting the stop codon and 3' UTR truncates the alternate CDS."""
        transcript = forward_transcript_with_utr
        expected_head = transcript.spliced_seq[9:42]  # ATG + the poly-G run, unedited
        variation_allele = _create_variation_allele(transcript, 52, 74, "N" * 23, "-")
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert tva.reference_cds_peptide == "M" + "G" * 10 + "*"
        assert tva.alternate_cds == expected_head
        assert tva.alternate_cds_peptide == "M" + "G" * 10

    def test_forward_stop_codon_delins_end_bound_excludes_the_replacement(
        self, forward_transcript_with_utr: Transcript
    ) -> None:
        """Test that a stop-spanning delins excludes sequence beyond the disrupted CDS end."""
        transcript = forward_transcript_with_utr
        expected_head = transcript.spliced_seq[9:42]  # ATG + the poly-G run, unedited
        variation_allele = _create_variation_allele(transcript, 52, 55, "TAAN", "A" * 10)
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert tva.alternate_cds == expected_head

    def test_forward_whole_cds_deletion(self, forward_transcript_with_utr: Transcript) -> None:
        """Test that deleting the complete CDS produces an empty alternate CDS."""
        transcript = forward_transcript_with_utr
        variation_allele = _create_variation_allele(transcript, 10, 54, CODING_SEQUENCE, "-")
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert tva.alternate_cds == ""
        assert not tva.alternate_cds_peptide.startswith("M")

    def test_insertion_inside_cds(self, forward_transcript_with_utr: Transcript) -> None:
        """Test that an in-frame CDS insertion adds a residue to the alternate peptide."""
        transcript = forward_transcript_with_utr
        variation_allele = _create_variation_allele(transcript, 15, 14, "-", "GGG")
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert tva.alternate_cds_peptide == "M" + "G" * 11 + "*"
        assert len(tva.alternate_cds_peptide) > len(tva.reference_cds_peptide)

    def test_purely_intronic_snv_is_a_no_op(self, transcript: Transcript) -> None:
        """Test that an intronic SNV leaves the CDS unchanged."""
        variation_allele = _create_variation_allele(transcript, 24, 24, "N", "A")
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert tva.alternate_cds == tva.reference_cds == CODING_SEQUENCE

    def test_variant_entirely_outside_the_transcript_is_a_no_op(self) -> None:
        """Test that a variant outside the transcript leaves the CDS unchanged."""
        transcript = _build_forward_transcript(genome_length=120)
        variation_allele = _create_variation_allele(transcript, 100, 100, "N", "A")
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert tva.alternate_cds == tva.reference_cds == CODING_SEQUENCE

    def test_reverse_no_op_5utr_snv(self, reverse_coding_transcript: Transcript) -> None:
        """Test that a reverse-strand 5' UTR SNV leaves the peptide unchanged."""
        transcript = reverse_coding_transcript
        utr_position = transcript.spliced_exons[0].end - 2
        variation_allele = _create_variation_allele(transcript, utr_position, utr_position, "N", "A")
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert tva.alternate_cds_peptide == tva.reference_cds_peptide

    def test_reverse_start_codon_destroyed(self, reverse_coding_transcript: Transcript) -> None:
        """Test that a reverse-strand deletion spanning the start codon removes methionine."""
        transcript = reverse_coding_transcript
        coding_start = transcript.coding_start
        assert coding_start is not None
        variation_allele = _create_variation_allele(
            transcript, coding_start - 1, coding_start + 2, "NNNN", "-"
        )
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert not tva.alternate_cds_peptide.startswith("M")

    def test_no_spliced_exons(self) -> None:
        """Test that a transcript without exons has no reference or alternate CDS."""
        transcript = Transcript(start=1, end=10, strand=Strand.FORWARD, stable_id="EMPTY")
        variation_allele = _create_variation_allele(transcript, 5, 5, "N", "A")
        tva = _create_transcript_variation_allele(transcript, variation_allele)
        assert tva.reference_cds is None
        assert tva.alternate_cds is None


class TestFullCdsAllelesReadthrough:
    """Test CDS and peptide properties with readthrough options."""

    @pytest.mark.parametrize(
        ("start", "end", "reference", "alternate", "label"),
        [
            (15, 15, "G", "A", "mid-CDS SNV"),
            (15, 17, "GGG", "-", "in-frame mid-CDS deletion"),
            (5, 5, "N", "A", "5' UTR SNV"),
        ],
    )
    def test_readthrough_is_a_no_op_when_the_stop_survives(
        self,
        forward_readthrough_transcript: Transcript,
        start: int,
        end: int,
        reference: str,
        alternate: str,
        label: str,
    ) -> None:
        """Test that readthrough levels have no effect when the stop codon survives."""
        transcript = forward_readthrough_transcript
        variation_allele = _create_variation_allele(transcript, start, end, reference, alternate)
        results = [
            _create_transcript_variation_allele(transcript, variation_allele, readthrough=level).alternate_cds
            for level in Readthrough
        ]
        assert results[0] == results[1] == results[2], label

    def test_utr_readthrough_finds_a_stop_in_the_three_prime_utr(
        self, forward_readthrough_transcript: Transcript
    ) -> None:
        """Test that UTR readthrough finds the next in-frame stop in the 3' UTR."""
        transcript = forward_readthrough_transcript
        variation_allele = _create_variation_allele(transcript, 52, 52, "T", "A")

        default = _create_transcript_variation_allele(transcript, variation_allele)
        assert default.alternate_cds_peptide == "M" + "G" * 10 + "K"
        assert "*" not in default.alternate_cds_peptide

        readthrough = _create_transcript_variation_allele(
            transcript, variation_allele, readthrough=Readthrough.UTR
        )
        assert readthrough.alternate_cds_peptide == "M" + "G" * 10 + "K*"
        assert readthrough.alternate_cds == default.alternate_cds + "TAA"

    def test_utr_readthrough_does_not_read_past_the_transcript(self) -> None:
        """Test that UTR readthrough stops at the end of the spliced transcript."""
        transcript = _build_forward_transcript(
            three_prime_utr="N" * 20, downstream_extra="N" + "TAA" + "AAAA"
        )
        variation_allele = _create_variation_allele(transcript, 52, 52, "T", "A")
        readthrough = _create_transcript_variation_allele(
            transcript, variation_allele, readthrough=Readthrough.UTR
        )

        assert len(readthrough.alternate_cds) == 56  # anchor 9 to the end of the 65-base cDNA
        assert "*" not in readthrough.alternate_cds_peptide

    def test_beyond_transcript_reads_flanking_sequence(self) -> None:
        """Test that beyond-transcript readthrough finds a stop in downstream sequence."""
        transcript = _build_forward_transcript(
            three_prime_utr="N" * 20, downstream_extra="N" + "TAA" + "AAAA"
        )
        variation_allele = _create_variation_allele(transcript, 52, 52, "T", "A")
        beyond = _create_transcript_variation_allele(
            transcript, variation_allele, readthrough=Readthrough.BEYOND_TRANSCRIPT
        )

        assert beyond.alternate_cds_peptide.endswith("*")
        assert len(beyond.alternate_cds) == 60  # 56 in-transcript + "N" completing the codon + TAA

    def test_beyond_transcript_is_bounded(self, forward_long_readthrough_transcript: Transcript) -> None:
        """Test that beyond-transcript readthrough is limited to the maximum distance."""
        transcript = forward_long_readthrough_transcript
        variation_allele = _create_variation_allele(transcript, 52, 52, "T", "A")
        beyond = _create_transcript_variation_allele(
            transcript, variation_allele, readthrough=Readthrough.BEYOND_TRANSCRIPT
        )

        assert "*" not in beyond.alternate_cds_peptide
        assert len(beyond.alternate_cds) <= 56 + MAX_READTHROUGH_BASES

    def test_reference_cds_is_level_invariant(self, forward_readthrough_transcript: Transcript) -> None:
        """Test that the reference CDS is unchanged by the readthrough level."""
        transcript = forward_readthrough_transcript
        variation_allele = _create_variation_allele(transcript, 52, 52, "T", "A")
        references = {
            _create_transcript_variation_allele(transcript, variation_allele, readthrough=level).reference_cds
            for level in Readthrough
        }
        assert references == {CODING_SEQUENCE}

    def test_every_accessor_agrees_at_a_non_default_level(
        self, forward_readthrough_transcript: Transcript
    ) -> None:
        """Test that nucleotide and peptide accessors use the same readthrough level."""
        transcript = forward_readthrough_transcript
        variation_allele = _create_variation_allele(transcript, 52, 52, "T", "A")
        tva = _create_transcript_variation_allele(transcript, variation_allele, readthrough=Readthrough.UTR)

        assert tva.alternate_cds_peptide == translate_dna(tva.alternate_cds)
        assert tva.reference_cds_peptide == translate_dna(tva.reference_cds)

    def test_reassigning_readthrough_invalidates_cached_alleles(
        self, forward_readthrough_transcript: Transcript
    ) -> None:
        """Test that reassigning readthrough invalidates cached allele sequences."""
        transcript = forward_readthrough_transcript
        variation_allele = _create_variation_allele(transcript, 52, 52, "T", "A")
        tva = _create_transcript_variation_allele(transcript, variation_allele)

        assert tva.alternate_cds_peptide == "M" + "G" * 10 + "K"  # cache the NONE answer first
        tva.readthrough = Readthrough.UTR
        assert tva.alternate_cds_peptide == "M" + "G" * 10 + "K*"

    def test_downstream_flanking_sequence_reverse_strand(self) -> None:
        """Test that reverse-strand downstream sequence is reverse-complemented."""
        slice = DummySlice("ACGTACGTACGTACGTACG" + "N" * 30)
        transcript = Transcript(start=20, end=30, strand=Strand.REVERSE, stable_id="RT", slice=slice)
        transcript.spliced_exons = [
            Exon(start=20, end=30, strand=Strand.REVERSE, stable_id="RT_e0", slice=slice)
        ]
        variation_allele = _create_variation_allele(transcript, 20, 20, "T", "A")
        tva = _create_transcript_variation_allele(transcript, variation_allele)

        flank = tva._downstream_flanking_sequence(10)  # pylint: disable=protected-access
        assert flank == reverse_complement(slice.seq[9:19])

    def test_downstream_flanking_sequence_is_empty_past_the_slice(
        self, forward_transcript_with_utr: Transcript
    ) -> None:
        """Test that missing downstream slice sequence returns an empty string."""
        transcript = forward_transcript_with_utr
        variation_allele = _create_variation_allele(transcript, 52, 52, "T", "A")
        tva = _create_transcript_variation_allele(transcript, variation_allele)

        assert tva._downstream_flanking_sequence(50) == ""  # pylint: disable=protected-access


class TestSequenceEditGeometry:
    """Test CDS and cDNA edit geometry."""

    def test_length_change_signs(self, forward_transcript_with_utr: Transcript) -> None:
        """Test length-change signs for substitutions, deletions, and insertions."""
        transcript = forward_transcript_with_utr
        substitution = _create_variation_allele(transcript, 15, 15, "G", "A")
        deletion = _create_variation_allele(transcript, 15, 17, "GGG", "-")
        insertion = _create_variation_allele(transcript, 15, 14, "-", "GGG")

        assert _create_transcript_variation_allele(transcript, substitution).cds_length_change == 0
        assert _create_transcript_variation_allele(transcript, deletion).cds_length_change == -3
        assert _create_transcript_variation_allele(transcript, insertion).cds_length_change == 3

    def test_cds_and_cdna_edits_are_not_interchangeable(
        self, forward_transcript_with_utr: Transcript
    ) -> None:
        """Test that a CDS-start variant has distinct CDS and cDNA edit footprints."""
        transcript = forward_transcript_with_utr
        # 5' UTR is genomic 1-9 and the CDS starts at 10, so this removes 2 UTR + 3 CDS bases.
        variation_allele = _create_variation_allele(transcript, 8, 12, "NNATG", "-")
        tva = _create_transcript_variation_allele(transcript, variation_allele)

        assert tva._cds_edit.footprint == (0, 2)  # pylint: disable=protected-access
        assert tva._cdna_edit.footprint == (7, 11)  # pylint: disable=protected-access
        assert tva._cds_edit.footprint != tva._cdna_edit.footprint  # pylint: disable=protected-access

    def test_cds_length_change_reports_coding_bases_not_cdna_bases(
        self, forward_transcript_with_utr: Transcript
    ) -> None:
        """Test that CDS length change counts only coding bases in a boundary deletion."""
        transcript = forward_transcript_with_utr
        variation_allele = _create_variation_allele(transcript, 8, 12, "NNATG", "-")
        tva = _create_transcript_variation_allele(transcript, variation_allele)

        assert tva.cds_length_change == -3
        assert tva._cdna_edit.length_change == -5  # pylint: disable=protected-access
        assert tva.cds_length_change % 3 == 0  # in frame, no abs() required

    def test_no_cds_overlap_gives_none(self, forward_transcript_with_utr: Transcript) -> None:
        """Test that variants outside the CDS have no coding length change."""
        transcript = forward_transcript_with_utr
        utr_variation_allele = _create_variation_allele(transcript, 3, 3, "N", "A")
        intron_variation_allele = _create_variation_allele(transcript, 24, 24, "N", "A")

        assert _create_transcript_variation_allele(transcript, utr_variation_allele).cds_length_change is None
        assert (
            _create_transcript_variation_allele(transcript, intron_variation_allele).cds_length_change is None
        )

    def test_negative_length_change_needs_no_abs_for_the_inframe_test(self) -> None:
        """Test that negative in-frame length changes are divisible by three without normalization."""
        for length_change in (-6, -3, 0, 3, 6):
            assert length_change % 3 == 0
        for length_change in (-5, -4, -2, -1, 1, 2):
            assert length_change % 3 != 0


# --------------------------------------------------------------------------
# Stage 2 parity: the reuse-based _full_cds_alleles/_full_peptide_alleles must match the
# pre-Stage-2 always-slice, always-retranslate implementation kept here as an oracle.
# --------------------------------------------------------------------------


def _oracle_nucleotide_sequence_up_to_stop(nucleotide_sequence: str) -> str:
    """Pre-Stage-2 truncate-only helper: always retranslates rather than returning the peptide
    the scan already computed."""
    peptide = translate_dna(nucleotide_sequence)
    stop_index = peptide.find("*")
    if stop_index == -1:
        return nucleotide_sequence
    return nucleotide_sequence[: (stop_index + 1) * 3]


def _oracle_full_cds_alleles(tva: TranscriptVariationAllele) -> tuple[str, str] | None:
    """Reimplementation of the pre-Stage-2 ``_full_cds_alleles``: always slices ``spliced_seq``
    for the reference half (never ``transcript.cds.seq``) and always retranslates the alternate
    tail to check the BEYOND_TRANSCRIPT stop condition. The oracle the parity tests below check
    the real implementation against.
    """
    transcript = tva.feature
    cdna_mapper = transcript.cdna_mapper
    spliced_seq = transcript.spliced_seq
    cds_length = transcript.cds_length
    cdna_coding_start = transcript.cdna_coding_start
    if cdna_mapper is None or spliced_seq is None or not cds_length or cdna_coding_start is None:
        return None

    cdna_coding_start0 = cdna_coding_start - 1
    cdna_coding_end0 = cdna_coding_start0 + cds_length
    reference_cds = spliced_seq[cdna_coding_start0:cdna_coding_end0]

    edit = tva._cdna_edit  # pylint: disable=protected-access
    if edit is None:
        return reference_cds, reference_cds
    edit_start, edit_end, transcript_sense_alternate = edit.start, edit.end, edit.replacement

    alternate_spliced_seq = spliced_seq[:edit_start] + transcript_sense_alternate + spliced_seq[edit_end:]
    length_change = edit.length_change

    anchor = _edited_boundary(cdna_coding_start0, edit_start, edit_end, length_change)
    if tva.readthrough >= Readthrough.UTR:
        end_bound = len(alternate_spliced_seq)
    else:
        end_bound = _edited_boundary(cdna_coding_end0, edit_start, edit_end, length_change)

    alternate_tail = alternate_spliced_seq[anchor:end_bound]
    alternate_cds = _oracle_nucleotide_sequence_up_to_stop(alternate_tail)

    if tva.readthrough >= Readthrough.BEYOND_TRANSCRIPT and "*" not in translate_dna(alternate_cds):
        flank = tva._downstream_flanking_sequence(MAX_READTHROUGH_BASES)  # pylint: disable=protected-access
        if flank:
            alternate_cds = _oracle_nucleotide_sequence_up_to_stop(alternate_tail + flank)

    return reference_cds, alternate_cds


def _oracle_full_peptide_alleles(tva: TranscriptVariationAllele) -> tuple[str, str] | None:
    """Pre-Stage-2 ``_full_peptide_alleles``: always translates both halves from scratch."""
    alleles = _oracle_full_cds_alleles(tva)
    if alleles is None:
        return None
    reference_cds, alternate_cds = alleles
    return translate_dna(reference_cds), translate_dna(alternate_cds)


def _assert_matches_oracle(tva: TranscriptVariationAllele) -> None:
    """Assert every full-CDS/peptide accessor agrees with the pre-Stage-2 oracle above."""
    assert (tva.reference_cds, tva.alternate_cds) == _oracle_full_cds_alleles(tva)
    assert (tva.reference_cds_peptide, tva.alternate_cds_peptide) == _oracle_full_peptide_alleles(tva)


class TestFullCdsAllelesParityWithOracle:
    """The reuse-based implementation must be indistinguishable from the pre-Stage-2 oracle."""

    @pytest.mark.parametrize(
        ("start", "end", "reference", "alternate"),
        [
            pytest.param(15, 15, "G", "A", id="mid_cds_snv"),
            pytest.param(10, 12, "ATG", "-", id="start_codon_deleted"),
            pytest.param(52, 52, "T", "A", id="stop_codon_destroyed_snv"),
            pytest.param(15, 17, "GGG", "-", id="inframe_mid_cds_deletion"),
            pytest.param(15, 15, "G", "-", id="frameshift_deletion"),
        ],
    )
    def test_parity_across_readthrough_levels(
        self,
        forward_readthrough_transcript: Transcript,
        start: int,
        end: int,
        reference: str,
        alternate: str,
    ) -> None:
        """Every readthrough level agrees with the oracle, for a range of edit shapes."""
        transcript = forward_readthrough_transcript
        variation_allele = _create_variation_allele(transcript, start, end, reference, alternate)
        for level in Readthrough:
            tva = _create_transcript_variation_allele(transcript, variation_allele, readthrough=level)
            _assert_matches_oracle(tva)

    def test_parity_with_no_downstream_stop(self, forward_long_readthrough_transcript: Transcript) -> None:
        """A frameshift with no stop anywhere downstream (not even past the transcript) is where
        the old and new BEYOND_TRANSCRIPT "*" check would most plausibly diverge."""
        transcript = forward_long_readthrough_transcript
        variation_allele = _create_variation_allele(transcript, 15, 15, "G", "-")
        for level in Readthrough:
            tva = _create_transcript_variation_allele(transcript, variation_allele, readthrough=level)
            _assert_matches_oracle(tva)

    def test_parity_reverse_strand(self, reverse_coding_transcript: Transcript) -> None:
        """The reverse-strand fixture agrees too, where cDNA offsets run against genomic order."""
        transcript = reverse_coding_transcript
        coding_start = transcript.coding_start
        assert coding_start is not None
        variation_allele = _create_variation_allele(
            transcript, coding_start - 1, coding_start + 2, "NNNN", "-"
        )
        for level in Readthrough:
            tva = _create_transcript_variation_allele(transcript, variation_allele, readthrough=level)
            _assert_matches_oracle(tva)

    def test_parity_n_filled_exon(self) -> None:
        """An N-filled exon (no attached sequence) still agrees with the oracle."""
        exon_one_seq, exon_three_seq = "ATGGGA", "CCCTAA"
        slice_ = DummySlice(exon_one_seq + "N" * 9 + exon_three_seq)
        transcript = Transcript(start=1, end=21, strand=Strand.FORWARD, stable_id="ENST_NFILL", slice=slice_)
        transcript.spliced_exons = [
            Exon(start=1, end=6, strand=Strand.FORWARD, stable_id="ENST_NFILL_e0", slice=slice_),
            Exon(start=10, end=12, strand=Strand.FORWARD, stable_id="ENST_NFILL_e1"),  # no slice -> N-filled
            Exon(start=16, end=21, strand=Strand.FORWARD, stable_id="ENST_NFILL_e2", slice=slice_),
        ]
        transcript.seed_attrs(
            first_coding_exon_index=0, last_coding_exon_index=2, exon_coding_start=1, exon_coding_end=6
        )
        variation_allele = _create_variation_allele(transcript, 2, 2, "T", "C")
        for level in Readthrough:
            tva = _create_transcript_variation_allele(transcript, variation_allele, readthrough=level)
            _assert_matches_oracle(tva)

    def test_parity_partial_last_codon(self) -> None:
        """A CDS whose length is not a multiple of three still agrees with the oracle."""
        genomic_sequence = "ATGGGANNNCCCT"  # 13 bases: 4 whole codons plus a dangling base
        slice_ = DummySlice(genomic_sequence)
        transcript = Transcript(
            start=1, end=len(genomic_sequence), strand=Strand.FORWARD, stable_id="ENST_PARTIAL", slice=slice_
        )
        transcript.spliced_exons = [
            Exon(start=1, end=13, strand=Strand.FORWARD, stable_id="ENST_PARTIAL_e0", slice=slice_),
        ]
        transcript.seed_attrs(
            first_coding_exon_index=0, last_coding_exon_index=0, exon_coding_start=1, exon_coding_end=13
        )
        variation_allele = _create_variation_allele(transcript, 5, 5, "G", "A")
        for level in Readthrough:
            tva = _create_transcript_variation_allele(transcript, variation_allele, readthrough=level)
            _assert_matches_oracle(tva)


class TestReferenceCdsPeptideReuse:
    """Guard behaviour of ``TranscriptVariationAllele._reference_cds_peptide``."""

    def test_reuses_translated_seq_when_translation_is_not_yet_resolved(self, transcript: Transcript) -> None:
        """The common case: nothing has touched ``transcript.translation`` yet."""
        assert not is_resolved(transcript, "_translation")
        variation_allele = _create_variation_allele(transcript, 5, 5, "N", "A")
        tva = _create_transcript_variation_allele(transcript, variation_allele)

        assert tva.reference_cds_peptide == transcript.translated_seq == "M" + "G" * 10 + "*"

    def test_reuses_translated_seq_when_the_seeded_translation_has_no_seq(
        self, transcript: Transcript
    ) -> None:
        """Mirrors the GFF/GTF builder: a seeded translation whose own ``seq`` is None."""
        derived = transcript.translation  # resolve _translation the ordinary way first
        assert derived is not None
        seeded = Translation(
            start=derived.start,
            end=derived.end,
            stable_id="seeded",
            exon_coding_start=derived.exon_coding_start,
            exon_coding_end=derived.exon_coding_end,
        )
        transcript.seed_attrs(translation=seeded)
        assert transcript.translation is seeded
        assert transcript.translation.seq is None

        variation_allele = _create_variation_allele(transcript, 5, 5, "N", "A")
        tva = _create_transcript_variation_allele(transcript, variation_allele)

        assert tva.reference_cds_peptide == transcript.translated_seq == "M" + "G" * 10 + "*"

    def test_falls_back_to_retranslating_for_a_seeded_curated_peptide(self, transcript: Transcript) -> None:
        """A future source could seed a curated peptide of its own; reuse must not trust it."""
        curated = Translation(start=1, end=36, stable_id="curated")
        curated.seq = "CURATED_PEPTIDE_STRING"
        transcript.seed_attrs(translation=curated)
        assert transcript.translated_seq == "CURATED_PEPTIDE_STRING"

        variation_allele = _create_variation_allele(transcript, 5, 5, "N", "A")
        tva = _create_transcript_variation_allele(transcript, variation_allele)

        assert tva.reference_cds_peptide == "M" + "G" * 10 + "*"
        assert tva.reference_cds_peptide != transcript.translated_seq

    def test_falls_back_when_cached_cds_length_disagrees(self, transcript: Transcript) -> None:
        """A CDS whose cached ``seq`` disagrees in length with ``cds_length`` cannot be trusted,
        for either the reference nucleotides or the reference peptide."""
        cds = transcript.cds
        assert cds is not None and cds.seq is not None
        cds.seq = cds.seq[:-3]  # now len(cds.seq) != transcript.cds_length

        variation_allele = _create_variation_allele(transcript, 5, 5, "N", "A")
        tva = _create_transcript_variation_allele(transcript, variation_allele)

        assert tva.reference_cds == CODING_SEQUENCE
        assert tva.reference_cds_peptide == "M" + "G" * 10 + "*"

    def test_reuses_transcript_cds_seq_for_reference_cds(self, transcript: Transcript) -> None:
        """When the guard holds, ``reference_cds`` is the very string ``transcript.cds.seq`` is,
        not a freshly sliced copy."""
        cds = transcript.cds
        assert cds is not None and cds.seq is not None
        variation_allele = _create_variation_allele(transcript, 5, 5, "N", "A")
        tva = _create_transcript_variation_allele(transcript, variation_allele)

        assert tva.reference_cds is cds.seq

    def test_reuse_does_not_depend_on_translation_resolution_order(self, transcript: Transcript) -> None:
        """A transcript whose ``_translation`` was already resolved by an earlier overlap still
        gets the reuse benefit for a later one.

        The guard checks ``Translation.derived_from_coding_sequence`` rather than *when*
        ``_translation`` resolved, so reuse no longer depends on which overlap happens to touch it
        first -- unlike a resolution-order-based guard, which would only ever benefit the first
        overlap on a transcript and force every later one to retranslate.
        """
        first_variant = _create_variation_allele(transcript, 5, 5, "N", "A")
        first_tva = _create_transcript_variation_allele(transcript, first_variant)
        assert first_tva.reference_cds_peptide == "M" + "G" * 10 + "*"  # resolves _translation
        assert is_resolved(transcript, "_translation")
        assert transcript.translation is not None
        assert transcript.translation.derived_from_coding_sequence is True

        second_variant = _create_variation_allele(transcript, 15, 15, "G", "A")
        second_tva = _create_transcript_variation_allele(transcript, second_variant)

        assert second_tva.reference_cds_peptide is transcript.translated_seq


class TestFullCdsAllelesInvalidation:
    """After the transcript's coding derivations change, a fresh TVA reflects the new state."""

    def test_picks_up_recomputed_cds_when_translation_was_never_resolved(self) -> None:
        """``_cds``/``_translated_seq`` are in ``Transcript._CODING_DERIVATIONS`` and invalidated
        by reassigning ``spliced_exons``; when nothing had resolved ``_translation`` beforehand,
        the fresh TVA's peptide is the freshly (re)derived ``transcript.translated_seq``."""
        transcript = _build_forward_transcript()
        variation_allele = _create_variation_allele(transcript, 5, 5, "N", "A")

        before = _create_transcript_variation_allele(transcript, variation_allele)
        assert before.reference_cds == CODING_SEQUENCE
        assert not is_resolved(transcript, "_translation")

        new_coding_sequence = "ATG" + "GGG" * 4 + "GCG" + "GGG" * 5 + "TAA"  # one residue swapped
        _reassign_forward_coding_exons(transcript, new_coding_sequence)

        after = _create_transcript_variation_allele(transcript, variation_allele)
        assert after.reference_cds == new_coding_sequence
        assert after.reference_cds != before.reference_cds
        assert after.reference_cds_peptide == transcript.translated_seq
        assert after.reference_cds_peptide == translate_dna(new_coding_sequence)

    def test_reuses_stale_translated_seq_after_spliced_exons_change_once_translation_resolved(
        self,
    ) -> None:
        """A known trade-off of order-independent reuse (``Translation.derived_from_coding_sequence``,
        not *when* ``_translation`` resolved): once resolved the plain way, it stays trusted even
        after a later ``spliced_exons`` reassignment, because ``_translation`` is deliberately absent
        from ``Transcript._CODING_DERIVATIONS`` (its own docstring: "an explicitly seeded translation
        must survive a later cds seed and vice versa").

        This is not something the reuse fix introduces: ``transcript.translated_seq`` itself already
        reads this same stale ``translation.seq`` regardless of any TVA-level reuse (asserted below).
        It also never arises in production, where a transcript's ``spliced_exons`` is assigned once,
        at build time, and never reassigned afterwards. ``reference_cds`` itself (``cds.seq``, which
        *is* invalidated by the reassignment) is never stale.
        """
        transcript = _build_forward_transcript()
        variation_allele = _create_variation_allele(transcript, 5, 5, "N", "A")

        before = _create_transcript_variation_allele(transcript, variation_allele)
        assert before.reference_cds_peptide == "M" + "G" * 10 + "*"  # resolves _translation
        stale_translated_seq = transcript.translated_seq

        new_coding_sequence = "ATG" + "GGG" * 4 + "GCG" + "GGG" * 5 + "TAA"  # one residue swapped
        _reassign_forward_coding_exons(transcript, new_coding_sequence)

        assert transcript.translated_seq == stale_translated_seq  # translated_seq is stale already

        after = _create_transcript_variation_allele(transcript, variation_allele)
        assert after.reference_cds == new_coding_sequence  # cds.seq itself is never stale
        assert after.reference_cds_peptide == stale_translated_seq  # reused despite being stale


def _reassign_forward_coding_exons(transcript: Transcript, coding_sequence: str) -> None:
    """Rebuild ``transcript``'s two coding exons over fresh, never-before-accessed Exon/Slice
    objects carrying ``coding_sequence``, matching ``_build_forward_transcript``'s own layout.

    A fresh Exon/Slice pair (rather than mutating the existing slice's characters) sidesteps each
    Exon's own ``_seq`` cache, which a slice-content mutation alone would not invalidate.
    """
    genomic_chars = list("N" * 74)
    _write_sequence_at_position(genomic_chars, 10, coding_sequence[:11])
    _write_sequence_at_position(genomic_chars, 30, coding_sequence[11:])
    new_slice = DummySlice("".join(genomic_chars))
    transcript.spliced_exons = [
        Exon(start=1, end=20, strand=Strand.FORWARD, stable_id="ENST_FWD_e0", slice=new_slice),
        Exon(start=30, end=74, strand=Strand.FORWARD, stable_id="ENST_FWD_e1", slice=new_slice),
    ]
    # Reassigning spliced_exons invalidates every seeded coding-boundary slot too (they are all in
    # _CODING_DERIVATIONS), so they must be re-seeded exactly as _build_forward_transcript does.
    transcript.seed_attrs(
        first_coding_exon_index=0, last_coding_exon_index=1, exon_coding_start=10, exon_coding_end=25
    )


class TestShiftedProteinLocus:
    """``_shifted_protein_locus`` must not recurse into itself."""

    def test_returns_none_without_recursing(
        self, transcript_variation_allele: TranscriptVariationAllele
    ) -> None:
        """The stub returns None, matching every other not-yet-implemented shifted locus."""
        assert transcript_variation_allele.shifted_protein_locus is None


class TestNotYetImplementedShiftedAndOriginalLoci:
    """The sibling not-yet-implemented locus stubs all follow the same pattern as
    ``shifted_protein_locus``: a public accessor delegating to a cached property that returns
    None. Exercised directly since nothing else in the suite reaches them."""

    def test_shifted(self, transcript_variation_allele: TranscriptVariationAllele) -> None:
        """Test that shifted is False by default."""
        assert transcript_variation_allele.shifted is False

    def test_shifted_cdna_locus(self, transcript_variation_allele: TranscriptVariationAllele) -> None:
        """Test that shifted_cdna_locus is None by default."""
        assert transcript_variation_allele.shifted_cdna_locus is None

    def test_shifted_cds_locus(self, transcript_variation_allele: TranscriptVariationAllele) -> None:
        """Test that shifted_cds_locus is None by default."""
        assert transcript_variation_allele.shifted_cds_locus is None

    def test_original_cdna_locus(self, transcript_variation_allele: TranscriptVariationAllele) -> None:
        """Test that original_cdna_locus is None by default."""
        assert transcript_variation_allele.original_cdna_locus is None

    def test_original_cds_locus(self, transcript_variation_allele: TranscriptVariationAllele) -> None:
        """Test that original_cds_locus is None by default."""
        assert transcript_variation_allele.original_cds_locus is None
