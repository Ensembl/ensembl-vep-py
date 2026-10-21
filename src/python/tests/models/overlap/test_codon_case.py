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
"""Tests for how an overlap cases its codons for display."""

from __future__ import annotations

import pytest

from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.facts.base import get_fact
from ensembl.vep.common import Strand
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import VariationAllele
from tests.models.overlap.conftest import SINGLE_EXON_CDS_START, SingleExonTranscript


def _overlap(
    single_exon_transcript: SingleExonTranscript, offset: int, end_offset: int, reference: str, alternate: str
) -> TranscriptVariationAllele:
    """Return a variant spanning CDS offsets ``offset`` to ``end_offset`` of the single-exon transcript."""
    transcript = single_exon_transcript(None)
    variation_allele = VariationAllele(
        start=SINGLE_EXON_CDS_START + offset,
        end=SINGLE_EXON_CDS_START + end_offset,
        strand=Strand.FORWARD,
        slice=transcript.slice,
        reference_allele=reference,
        alternate_allele=alternate,
    )
    return TranscriptVariationAllele(feature=transcript, overlapping_feature=variation_allele)


class TestDisplayCase:
    """VEP's ``display_codon``: lower case but for each allele's own bases."""

    @pytest.mark.parametrize(
        ("offset", "end_offset", "reference", "alternate", "expected"),
        [
            pytest.param(4, 4, "G", "A", ("gGg", "gAg"), id="snv"),
            pytest.param(3, 4, "GG", "CA", ("GGg", "CAg"), id="two_base_substitution"),
            pytest.param(3, 5, "GGG", "-", ("GGG", ""), id="whole_codon_deletion"),
            pytest.param(4, 4, "G", "-", ("gGg", "gg"), id="frameshift_deletion"),
            pytest.param(4, 3, "-", "TTT", ("ggg", "gTTTgg"), id="insertion_inside_a_codon"),
            pytest.param(3, 2, "-", "TTT", ("", "TTT"), id="insertion_between_codons"),
        ],
    )
    def test_each_allele_is_upper_case_in_its_own_codon(
        self,
        single_exon_transcript: SingleExonTranscript,
        offset: int,
        end_offset: int,
        reference: str,
        alternate: str,
        expected: tuple[str, str],
    ) -> None:
        """A deleted base is upper case on the reference side and simply absent on the other."""
        overlap = _overlap(single_exon_transcript, offset, end_offset, reference, alternate)

        assert (overlap.reference_codon, overlap.alternate_codon) == expected

    def test_an_insertion_between_codons_has_no_reference_window(
        self, single_exon_transcript: SingleExonTranscript
    ) -> None:
        """VEP's window closes before it opens there, leaving the inserted bases on their own."""
        overlap = _overlap(single_exon_transcript, 3, 2, "-", "TTT")

        assert (overlap.reference_codon, overlap.alternate_codon) == ("", "TTT")
        assert (overlap.reference_peptide, overlap.alternate_peptide) == ("", "F")

    def test_peptides_are_translated_from_the_bases_alone(
        self, single_exon_transcript: SingleExonTranscript
    ) -> None:
        """Case is presentation only: ``gGg`` is still glycine."""
        overlap = _overlap(single_exon_transcript, 4, 4, "G", "A")

        assert (overlap.reference_peptide, overlap.alternate_peptide) == ("G", "E")

    def test_the_codon_fact_is_upper_case(self, single_exon_transcript: SingleExonTranscript) -> None:
        """Facts compare bases, so the fact layer upper-cases the displayed codons first."""
        overlap = _overlap(single_exon_transcript, 4, 4, "G", "A")

        assert get_fact(FactContext(overlap=overlap), "_codon_alleles") == ("GGG", "GAG")
