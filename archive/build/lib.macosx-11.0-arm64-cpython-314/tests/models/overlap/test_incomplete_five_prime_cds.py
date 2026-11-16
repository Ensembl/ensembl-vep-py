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
"""Tests for a CDS whose 5' end starts part-way through a codon."""

from __future__ import annotations

import pytest

from ensembl.vep.common import ExonPhase, Strand
from ensembl.vep.common.codons import translate_dna
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import VariationAllele
from tests.models.overlap.conftest import (
    SINGLE_EXON_CDS_START,
    SINGLE_EXON_CODING_SEQUENCE,
    SingleExonTranscript,
)


def _overlap(
    transcript: Transcript, position: int, reference: str, alternate: str
) -> TranscriptVariationAllele:
    """Return an SNV at ``position`` over ``transcript``."""
    variation_allele = VariationAllele(
        start=position,
        end=position,
        strand=Strand.FORWARD,
        slice=transcript.slice,
        reference_allele=reference,
        alternate_allele=alternate,
    )
    return TranscriptVariationAllele(feature=transcript, overlapping_feature=variation_allele)


class TestPaddedCodingSequence:
    """The CDS read in VEP's frame: padded with ``N`` to a whole first codon."""

    def test_one_base_of_padding_for_phase_one(self, single_exon_transcript: SingleExonTranscript) -> None:
        """Ensembl phase 1: one base of the first codon lies before the CDS."""
        transcript = single_exon_transcript(ExonPhase.ONE)

        assert transcript.cds_start_phase == 1
        assert transcript.coding_sequence == "N" + SINGLE_EXON_CODING_SEQUENCE
        assert transcript.genomic_to_cds_offset(SINGLE_EXON_CDS_START) == 1

    def test_only_the_padded_lengths_count_the_padding(
        self, single_exon_transcript: SingleExonTranscript
    ) -> None:
        """Which lengths include VEP's ``N`` padding and which count real bases only.

        ``cds.seq`` is genomic sequence, so it has no padding, while ``cds.nucleotide_length``
        comes from the padded ``cds_mapper``: the two differ by the start phase for a 5' incomplete
        CDS, and anything comparing them has to say which it means.
        """
        transcript = single_exon_transcript(ExonPhase.ONE)
        cds = transcript.cds

        assert cds is not None and cds.seq == SINGLE_EXON_CODING_SEQUENCE
        assert transcript.cds_length == len(SINGLE_EXON_CODING_SEQUENCE)
        assert transcript.cds_mapper is not None
        assert transcript.cds_mapper.length == len(SINGLE_EXON_CODING_SEQUENCE) + 1
        assert cds.nucleotide_length == len(cds.seq) + transcript.cds_start_phase
        assert len(transcript.coding_sequence or "") == cds.nucleotide_length

    @pytest.mark.parametrize("phase", [None, ExonPhase.ZERO, ExonPhase.UNKNOWN])
    def test_no_phase_is_a_complete_five_prime_end(
        self, single_exon_transcript: SingleExonTranscript, phase: ExonPhase | None
    ) -> None:
        """Only phase 1 or 2 on the first coding exon pads the CDS."""
        transcript = single_exon_transcript(phase)

        assert transcript.cds_start_phase == 0
        assert transcript.coding_sequence == SINGLE_EXON_CODING_SEQUENCE

    def test_the_translation_reads_the_padded_frame(
        self, single_exon_transcript: SingleExonTranscript
    ) -> None:
        """The partial first codon translates to ``X``, and the rest in the shifted frame."""
        transcript = single_exon_transcript(ExonPhase.TWO)

        assert transcript.translated_seq == translate_dna("NN" + SINGLE_EXON_CODING_SEQUENCE)
        assert transcript.translated_seq is not None and transcript.translated_seq.startswith("X")


class TestPaddedOverlap:
    """Codons and the full CDS of an overlap, in the padded frame."""

    def test_codons_are_read_in_the_padded_frame(self, single_exon_transcript: SingleExonTranscript) -> None:
        """The start codon's third base opens codon 1 once the frame shifts by one."""
        overlap = _overlap(single_exon_transcript(ExonPhase.ONE), SINGLE_EXON_CDS_START + 2, "G", "A")

        assert (overlap.reference_codon, overlap.alternate_codon) == ("Ggg", "Agg")
        assert (overlap.reference_peptide, overlap.alternate_peptide) == ("G", "R")

    def test_the_full_cds_carries_the_same_padding(
        self, single_exon_transcript: SingleExonTranscript
    ) -> None:
        """Reference and alternate CDS stay aligned with the codon window."""
        overlap = _overlap(single_exon_transcript(ExonPhase.TWO), SINGLE_EXON_CDS_START + 5, "G", "A")

        assert overlap.reference_cds == "NN" + SINGLE_EXON_CODING_SEQUENCE
        assert overlap.alternate_cds is not None
        assert overlap.alternate_cds.startswith("NN")
        assert len(overlap.alternate_cds) == len(overlap.reference_cds)
        assert overlap.reference_cds_peptide == translate_dna("NN" + SINGLE_EXON_CODING_SEQUENCE)
