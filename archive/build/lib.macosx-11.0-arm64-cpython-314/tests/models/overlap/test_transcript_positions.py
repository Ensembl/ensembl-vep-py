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
"""Tests for an overlap's cDNA, CDS and protein positions."""

from __future__ import annotations

import pytest

from ensembl.vep.common import CoordSystemLevel, ExonPhase, Strand
from ensembl.vep.models.genic import Exon, Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.overlap.transcript_variation_allele import TranscriptPosition
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.models.variation import VariationAllele
from ensembl.vep.output.factory.write.base import OutputFactory


def _transcript(strand: Strand, phase: ExonPhase | None = None) -> Transcript:
    """Return a two-exon transcript, exons at 1-20 and 31-60, 3 bases of UTR at each end.

    Forward it codes 4-20 and 31-57; reverse, 57-31 and 20-4. ``phase`` goes on the first coding
    exon.
    """
    genomic_slice = Slice(
        coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
        seq_region_name="chrPositions",
        start=1,
        end=100,
        strand=Strand.FORWARD,
    )
    low_exon = Exon(start=1, end=20, strand=strand, stable_id="ENST_POS_low", slice=genomic_slice)
    high_exon = Exon(start=31, end=60, strand=strand, stable_id="ENST_POS_high", slice=genomic_slice)
    exons = [low_exon, high_exon] if strand == Strand.FORWARD else [high_exon, low_exon]
    exons[0].phase = phase
    transcript = Transcript(start=1, end=60, strand=strand, stable_id="ENST_POS", slice=genomic_slice)
    transcript.spliced_exons = exons
    transcript.seed_attrs(
        first_coding_exon_index=0, last_coding_exon_index=1, exon_coding_start=4, exon_coding_end=17
    )
    return transcript


def _overlap(transcript: Transcript, start: int, end: int) -> TranscriptVariationAllele:
    """Return a variant over ``start``..``end``; ``end == start - 1`` is an insertion."""
    length = end - start + 1
    variation_allele = VariationAllele(
        start=start,
        end=end,
        strand=Strand.FORWARD,
        slice=transcript.slice,
        reference_allele="A" * length if length else "-",
        alternate_allele="T" * length if length else "TT",
    )
    return TranscriptVariationAllele(feature=transcript, overlapping_feature=variation_allele)


class TestForwardPositions:
    """Positions on a forward-strand transcript, 5' end first."""

    @pytest.mark.parametrize(
        ("start", "end", "cdna", "cds", "protein"),
        [
            pytest.param(10, 10, (10, 10), (7, 7), (3, 3), id="first_exon_snv"),
            pytest.param(35, 35, (25, 25), (22, 22), (8, 8), id="second_exon_snv"),
            pytest.param(18, 33, (18, 23), (15, 20), (5, 7), id="deletion_over_the_intron"),
            pytest.param(25, 33, (None, 23), (None, 20), (None, 7), id="five_prime_end_in_the_intron"),
            pytest.param(11, 10, (10, 11), (7, 8), (3, 3), id="insertion"),
            pytest.param(2, 2, (2, 2), None, None, id="five_prime_utr"),
            pytest.param(25, 25, None, None, None, id="intronic"),
        ],
    )
    def test_each_end_is_placed_or_marked_missing(
        self,
        start: int,
        end: int,
        cdna: TranscriptPosition | None,
        cds: TranscriptPosition | None,
        protein: TranscriptPosition | None,
    ) -> None:
        """A sequence the variant misses gives no position; one end off it gives None there."""
        overlap = _overlap(_transcript(Strand.FORWARD), start, end)

        assert overlap.cdna_position == cdna
        assert overlap.cds_position == cds
        assert overlap.protein_position == protein


class TestReversePositions:
    """On the reverse strand the 5' end is the high genomic coordinate."""

    @pytest.mark.parametrize(
        ("start", "end", "cdna", "cds"),
        [
            pytest.param(57, 57, (4, 4), (1, 1), id="first_coding_base"),
            pytest.param(55, 57, (4, 6), (1, 3), id="start_codon_deletion"),
            pytest.param(18, 33, (28, 33), (25, 30), id="deletion_over_the_intron"),
        ],
    )
    def test_positions_count_from_the_high_end(
        self, start: int, end: int, cdna: TranscriptPosition, cds: TranscriptPosition
    ) -> None:
        """The 5' end comes first whichever strand the transcript is on."""
        overlap = _overlap(_transcript(Strand.REVERSE), start, end)

        assert overlap.cdna_position == cdna
        assert overlap.cds_position == cds


class TestIncompleteFivePrime:
    """VEP counts CDS and protein positions from a whole first codon, the cDNA as it is."""

    def test_the_start_phase_shifts_cds_and_protein_positions(self) -> None:
        """Phase 1: one padding base before the CDS, so its first base is CDS position 2."""
        overlap = _overlap(_transcript(Strand.FORWARD, ExonPhase.ONE), 6, 6)

        assert overlap.cdna_position == (6, 6)
        assert overlap.cds_position == (4, 4)
        assert overlap.protein_position == (2, 2)


class TestOutput:
    """How the positions are written."""

    @pytest.mark.parametrize(
        ("start", "end", "expected"),
        [
            pytest.param(
                18,
                33,
                {"cDNA_position": "18-23", "CDS_position": "15-20", "Protein_position": "5-7"},
                id="range",
            ),
            pytest.param(
                25,
                33,
                {"cDNA_position": "?-23", "CDS_position": "?-20", "Protein_position": "?-7"},
                id="open_ended",
            ),
            pytest.param(2, 2, {"cDNA_position": "2"}, id="utr_only"),
            pytest.param(25, 25, {}, id="intronic"),
        ],
    )
    def test_positions_are_written_as_vep_formats_them(
        self, start: int, end: int, expected: dict[str, str]
    ) -> None:
        """Only the sequences the variant touches get a column; a missing end prints as ``?``."""
        transcript = _transcript(Strand.FORWARD)
        transcript.seed_attrs(translation=None)
        output_hash: dict[str, object] = {}

        # pylint: disable-next=protected-access
        OutputFactory._position_fields(output_hash, _overlap(transcript, start, end))

        assert output_hash == expected
