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
"""Tests for synthetic transcript construction."""

from __future__ import annotations

import pytest

from ensembl.vep.common import CoordSystemLevel, Strand
from ensembl.vep.corpus.synthetic.components.transcript import (
    build_cds_genomic_bounds,
    build_transcript,
    build_translation,
    calculate_cds_length,
    patch_terminal_codons,
)
from ensembl.vep.corpus.synthetic.default import build_genome
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.sequence import CoordSystem, Slice


class TestTranscriptComponent:
    """Transcript relationships created for the default genome."""

    def test_reverse_transcript_exons_follow_transcript_order(self) -> None:
        """Reverse-strand construction reverses genomic exons for transcript traversal."""
        _sequence_region, transcripts = build_genome()
        reverse = next(transcript for transcript in transcripts if transcript.strand == Strand.REVERSE)
        starts = [exon.start for exon in reverse.spliced_exons]
        assert starts == sorted(starts, reverse=True)

    def test_build_translation_uses_the_exons_containing_its_boundaries(self) -> None:
        """Build exon-relative offsets when the CDS starts after the first exon."""
        synthetic_slice = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="test",
            start=1,
            end=400,
            strand=Strand.FORWARD,
        )
        transcript = build_transcript(
            stable_id="ENST_TEST",
            synthetic_slice=synthetic_slice,
            strand=Strand.FORWARD,
            exons_genomic_ascending=((100, 119), (200, 299)),
            cds_genomic_bounds=None,
            gene=None,
        )

        translation = build_translation(transcript, (230, 289))

        assert translation.transcript is transcript
        assert translation.stable_id == "ENST_TEST_translation"
        assert translation.start == 131
        assert translation.end == 190
        assert translation.exon_coding_start == 31
        assert translation.exon_coding_end == 90

    def test_terminal_codon_outside_the_cds_is_rejected(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Reject a terminal codon when a transcript cannot map its CDS offsets."""
        _sequence_region, transcripts = build_genome()

        def outside_cds(
            _transcript: Transcript,
            _cds_offset: int,
            *,
            accelerator: object | None = None,
        ) -> None:
            """Return no genomic position for any CDS offset."""
            del accelerator

        monkeypatch.setattr(Transcript, "cds_offset_to_genomic", outside_cds)

        with pytest.raises(ValueError, match="terminal codon falls outside"):
            patch_terminal_codons(["A"] * 60000, transcripts)


class TestCdsBounds:
    """Reusable CDS coordinate placement."""

    @pytest.mark.parametrize(
        ("exons", "strand", "expected_start"),
        [
            (((100, 119), (200, 299)), Strand.FORWARD, 230),
            (((100, 199), (300, 319)), Strand.REVERSE, 169),
        ],
    )
    def test_start_offset_carries_into_the_next_exon(
        self,
        exons: tuple[tuple[int, int], ...],
        strand: Strand,
        expected_start: int,
    ) -> None:
        """A short 5-prime exon never places the start codon inside its intron."""
        low, high = build_cds_genomic_bounds(
            exons,
            strand=strand,
            frame_index=0,
            start_offset=50,
            stop_offset=10,
        )
        transcript_start = low if strand == Strand.FORWARD else high
        assert transcript_start == expected_start
        assert calculate_cds_length(exons, low, high) % 3 == 0

    @pytest.mark.parametrize(
        ("exons", "strand", "expected_start"),
        [
            (((100, 119), (200, 299)), Strand.FORWARD, 230),
            (((100, 199), (300, 319)), Strand.REVERSE, 169),
        ],
    )
    def test_short_first_exon_builds_a_valid_coding_transcript(
        self,
        exons: tuple[tuple[int, int], ...],
        strand: Strand,
        expected_start: int,
    ) -> None:
        """A CDS carried into the next exon retains valid model offsets and terminal codons."""
        chromosome_chars = ["C"] * 400
        synthetic_slice = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="test",
            start=1,
            end=len(chromosome_chars),
            strand=Strand.FORWARD,
        )
        with synthetic_slice.unsafe() as unsafe_slice:
            setattr(unsafe_slice, "_seq", "".join(chromosome_chars))
        genomic_bounds = build_cds_genomic_bounds(
            exons,
            strand=strand,
            frame_index=0,
            start_offset=50,
            stop_offset=10,
        )

        transcript = build_transcript(
            stable_id="ENST_TEST",
            synthetic_slice=synthetic_slice,
            strand=strand,
            exons_genomic_ascending=exons,
            cds_genomic_bounds=genomic_bounds,
            gene=None,
        )
        patch_terminal_codons(chromosome_chars, [transcript])
        with synthetic_slice.unsafe() as unsafe_slice:
            setattr(unsafe_slice, "_seq", "".join(chromosome_chars))

        assert transcript.first_coding_exon_index == 1
        assert transcript.coding_start == expected_start
        assert transcript.cds_length is not None
        assert transcript.cds_length % 3 == 0
        assert transcript.coding_sequence is not None
        assert transcript.coding_sequence.startswith("ATG")
        assert transcript.coding_sequence.endswith("TAA")

    def test_build_transcript_rejects_a_cds_boundary_in_an_intron(self) -> None:
        """Translation construction rejects CDS terminal bounds outside exons."""
        synthetic_slice = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="test",
            start=1,
            end=300,
            strand=Strand.FORWARD,
        )
        with pytest.raises(ValueError, match="must fall within transcript exons"):
            build_transcript(
                stable_id="ENST_TEST",
                synthetic_slice=synthetic_slice,
                strand=Strand.FORWARD,
                exons_genomic_ascending=((100, 119), (200, 219)),
                cds_genomic_bounds=(150, 219),
                gene=None,
            )

    def test_stop_offset_is_clamped_after_the_start(self) -> None:
        """An excessive 3-prime offset clamps after a start carried into the next exon."""
        exons = ((100, 109), (200, 209))
        low, high = build_cds_genomic_bounds(
            exons,
            strand=Strand.FORWARD,
            frame_index=1,
            start_offset=12,
            stop_offset=100,
        )
        assert (low, high) == (202, 204)
        assert calculate_cds_length(exons, low, high) == 3

    def test_too_short_transcript_is_rejected(self) -> None:
        """A layout that cannot hold a codon fails with a domain error."""
        with pytest.raises(ValueError, match="too short"):
            build_cds_genomic_bounds(
                ((100, 101),),
                strand=Strand.FORWARD,
                frame_index=0,
                start_offset=0,
                stop_offset=0,
            )

    def test_empty_exons_are_rejected(self) -> None:
        """A coding layout needs an exon before terminal offsets can be resolved."""
        with pytest.raises(ValueError, match="at least one exon"):
            build_cds_genomic_bounds((), strand=Strand.FORWARD, frame_index=0)

    @pytest.mark.parametrize(
        ("frame_index", "start_offset", "stop_offset"),
        [(-1, 0, 0), (0, -1, 0), (0, 0, -1)],
    )
    def test_negative_terminal_parameters_are_rejected(
        self,
        frame_index: int,
        start_offset: int,
        stop_offset: int,
    ) -> None:
        """Frame and terminal offsets are non-negative transcript distances."""
        with pytest.raises(ValueError, match="must not be negative"):
            build_cds_genomic_bounds(
                ((100, 109),),
                strand=Strand.FORWARD,
                frame_index=frame_index,
                start_offset=start_offset,
                stop_offset=stop_offset,
            )
