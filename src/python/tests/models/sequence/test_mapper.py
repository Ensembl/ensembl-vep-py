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
"""Unit tests for the Mapper class."""

from __future__ import annotations
import logging

import pytest

from ensembl.vep.common import Strand
from ensembl.vep.models.sequence import Mapper

# Two 10-base exon blocks with a 90-base gap between them: cds_length 20, offsets 0-9 in the
# first block, 10-19 in the second. genomic_spans is transcript order (Mapper's actual contract),
# which coincides with genomic-ascending order on the forward strand only.
FORWARD_GENOMIC_SPANS = [(100, 109), (200, 209)]
# Reverse strand: transcript order puts the genomically-*higher* block first (it is the 5'-most
# one in transcript sense) -- the mirror image of FORWARD_SPANS, not the same list.
REVERSE_GENOMIC_SPANS = [(200, 209), (100, 109)]


@pytest.fixture
def forward_mapper() -> Mapper:
    """Create a forward-strand mapper fixture."""
    return Mapper(
        genomic_spans=FORWARD_GENOMIC_SPANS,
        strand=Strand.FORWARD,
        coding_start_offset=0,
        coding_end_offset=19,
    )


@pytest.fixture
def reverse_mapper() -> Mapper:
    """Create a reverse-strand mapper fixture."""
    return Mapper(
        genomic_spans=REVERSE_GENOMIC_SPANS,
        strand=Strand.REVERSE,
        coding_start_offset=0,
        coding_end_offset=19,
    )


@pytest.fixture
def unknown_mapper() -> Mapper:
    """Create an unknow-strand mapper fixture."""
    return Mapper(
        genomic_spans=FORWARD_GENOMIC_SPANS,
        strand=Strand.UNKNOWN,
        coding_start_offset=0,
        coding_end_offset=19,
    )


@pytest.fixture
def junction_spanning_mapper() -> Mapper:
    """Create a mapper whose first codon spans two genomic regions."""
    return Mapper(
        genomic_spans=[(98, 99), (200, 217)],
        strand=Strand.FORWARD,
        coding_start_offset=0,
        coding_end_offset=19,
    )


@pytest.fixture
def multi_exon_mapper() -> Mapper:
    """Create a mapper over three forward-strand genomic regions."""
    return Mapper(
        genomic_spans=[(1, 10), (20, 30), (40, 50)],
        strand=Strand.FORWARD,
        coding_start_offset=0,
        coding_end_offset=30,
    )


@pytest.fixture
def partial_codon_mapper() -> Mapper:
    """Create a mapper whose sequence contains only a partial codon."""
    return Mapper(
        genomic_spans=[(100, 101)], strand=Strand.FORWARD, coding_start_offset=0, coding_end_offset=2
    )


class TestCoordinateMapping:
    """Test genomic and sequence coordinate mapping."""

    @pytest.mark.parametrize("cds_offset", range(20))
    def test_forward_mapper_round_trips_each_offset(self, forward_mapper: Mapper, cds_offset: int) -> None:
        """Test that forward-strand offsets round trip through genomic coordinates."""
        genomic_position = forward_mapper.offset_to_genomic(cds_offset)
        assert genomic_position is not None
        assert forward_mapper.genomic_to_offset(genomic_position) == cds_offset

    @pytest.mark.parametrize("cds_offset", range(20))
    def test_reverse_mapper_round_trips_each_offset(self, reverse_mapper: Mapper, cds_offset: int) -> None:
        """Test that reverse-strand offsets round trip through genomic coordinates."""
        genomic_position = reverse_mapper.offset_to_genomic(cds_offset)
        assert genomic_position is not None
        assert reverse_mapper.genomic_to_offset(genomic_position) == cds_offset

    @pytest.mark.parametrize("cds_offset", range(20))
    def test_unknown_mapper_round_trips_each_offset(
        self,
        unknown_mapper: Mapper,
        forward_mapper: Mapper,
        cds_offset: int,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Test that unknown-strand is treated as forward strand with warning."""
        with caplog.at_level(logging.WARNING):
            genomic_position = unknown_mapper.offset_to_genomic(cds_offset)
            assert "Strand is unknown" in caplog.text, "Incorrect warning in offset_to_genomic"
            assert genomic_position == forward_mapper.offset_to_genomic(cds_offset)

        genomic_position = forward_mapper.offset_to_genomic(cds_offset)
        assert genomic_position is not None
        with caplog.at_level(logging.WARNING):
            offset = unknown_mapper.genomic_to_offset(genomic_position)
            assert "Strand is unknown" in caplog.text, "Incorrect warning in genomic_to_offset"
            assert offset == cds_offset

    def test_reverse_offset_0_is_the_high_genomic_block(self, reverse_mapper: Mapper) -> None:
        """Test that reverse-strand offset zero maps to the higher genomic region."""
        assert reverse_mapper.offset_to_genomic(0) == 209

    def test_out_of_range_positions_map_to_none(self, forward_mapper: Mapper, reverse_mapper: Mapper) -> None:
        """Test that positions and offsets outside the mapped regions return None."""
        for position in (50, 110, 150, 199, 210, 300):
            assert forward_mapper.genomic_to_offset(position) is None
            assert reverse_mapper.genomic_to_offset(position) is None
        for offset in (-1, 20, 1000):
            assert forward_mapper.offset_to_genomic(offset) is None


class TestGenomicSpanToOffsets:
    """Test mapping codon offsets to genomic spans."""

    def test_forward_mapper_genomic_span_to_offsets(self, forward_mapper: Mapper) -> None:
        """Test forward_mapper.genomic_span_to_offsets gives the correct offsets"""
        assert forward_mapper.genomic_span_to_offsets(105, 109) == (5, 9)
        assert forward_mapper.genomic_span_to_offsets(105, 205) == (5, 15)
        assert forward_mapper.genomic_span_to_offsets(205, 105) is None
        assert forward_mapper.genomic_span_to_offsets(1, 10) is None
        # test span that starts before first mapped region is clamped
        assert forward_mapper.genomic_span_to_offsets(50, 100) == (0, 0)
        # test span start in intron advances to next mapped region
        assert forward_mapper.genomic_span_to_offsets(150, 200) == (10, 10)
        # test span entirely between mapped regions returns none
        assert forward_mapper.genomic_span_to_offsets(150, 199) is None

    def test_reverse_mapper_genomic_span_to_offsets(self, reverse_mapper: Mapper) -> None:
        """Test reverse_mapper.genomic_span_to_offsets gives the correct offsets"""
        assert reverse_mapper.genomic_span_to_offsets(105, 109) == (10, 14)
        assert reverse_mapper.genomic_span_to_offsets(105, 205) == (4, 14)

    def test_unknown_mapper_gives_warning(
        self, unknown_mapper: Mapper, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Test unkown_mapper gives a warning"""
        with caplog.at_level(logging.WARNING):
            offsets = unknown_mapper.genomic_span_to_offsets(105, 109)
            assert "Strand is unknown" in caplog.text, "Incorrect warning in genomic_to_offset"
            assert offsets == (5, 9)


class TestCodonGenomicSpans:
    """Test mapping codon offsets to genomic spans."""

    def test_non_straddling_codon_is_one_span(self, forward_mapper: Mapper) -> None:
        """Test that a codon within one genomic region maps to one span."""
        assert forward_mapper.codon_genomic_spans(0) == [(100, 102)]
        # cds_length is 20, not a multiple of 3 -- the last (partial) codon is offsets 18/19 only.
        assert forward_mapper.codon_genomic_spans(forward_mapper.length - 1) == [(208, 209)]

    def test_reverse_start_codon_is_the_high_genomic_block(self, reverse_mapper: Mapper) -> None:
        """Test that a reverse-strand start codon maps to the higher genomic region."""
        assert reverse_mapper.codon_genomic_spans(0) == [(207, 209)]
        assert reverse_mapper.codon_genomic_spans(reverse_mapper.length - 1) == [(100, 101)]

    def test_straddling_start_codon_is_two_spans(self, junction_spanning_mapper: Mapper) -> None:
        """Test that a start codon spanning a junction maps to two genomic spans."""
        assert junction_spanning_mapper.codon_genomic_spans(0) == [(98, 99), (200, 200)]

    def test_normalizes_to_the_codon_aligned_start(self, multi_exon_mapper: Mapper) -> None:
        """Test that each base in a codon returns the same genomic spans."""
        expected_spans = multi_exon_mapper.codon_genomic_spans(3)
        assert expected_spans == [(4, 6)]
        assert multi_exon_mapper.codon_genomic_spans(4) == expected_spans
        assert multi_exon_mapper.codon_genomic_spans(5) == expected_spans

    def test_partial_and_out_of_range_codon_offsets(self, partial_codon_mapper: Mapper) -> None:
        """Test mapping a partial codon and offsets outside the mapped sequence."""
        assert partial_codon_mapper.length == 2
        assert partial_codon_mapper.codon_genomic_spans(0) == [(100, 101)]
        # -1's codon-aligned start is -3 (floor division): offsets -3/-2/-1 are all out of range.
        assert partial_codon_mapper.codon_genomic_spans(-1) is None
        assert partial_codon_mapper.codon_genomic_spans(5) is None


class TestCodonWindow:
    """Test mapping codon offsets to codon windows."""

    def test_codon_window_one_span(self, forward_mapper: Mapper) -> None:
        """Test offsets within a codon maps to one codon"""
        assert forward_mapper.codon_window((4, 5)) == (3, 6)
        assert forward_mapper.codon_window((18, 20)) == (18, 20)
        assert forward_mapper.codon_window((18, 21)) == (18, 20)
        assert forward_mapper.codon_window((21, 22)) is None
        assert forward_mapper.codon_window((-1, 2)) == (0, 3)
        assert forward_mapper.codon_window((-5, -4)) is None

    def test_codon_window_two_spans(self, forward_mapper: Mapper) -> None:
        """Test offsets over two codons maps to two codons"""
        assert forward_mapper.codon_window((4, 7)) == (3, 9)


class TestProteinMapping:
    """Test mapping between sequence and protein offsets."""

    def test_cds_offset_to_protein_offset(self, forward_mapper: Mapper) -> None:
        """Test that sequence offsets map to the expected protein offsets."""
        assert forward_mapper.offset_to_protein_offset(0) == 0
        assert forward_mapper.offset_to_protein_offset(2) == 0
        assert forward_mapper.offset_to_protein_offset(3) == 1
        assert forward_mapper.offset_to_protein_offset(19) == 6
        assert forward_mapper.offset_to_protein_offset(20) is None
        assert forward_mapper.offset_to_protein_offset(-1) is None

    def test_protein_offset_to_cds_span(self, forward_mapper: Mapper) -> None:
        """Test that protein offsets map to the expected sequence spans."""
        assert forward_mapper.protein_offset_to_span(0) == (0, 2)
        assert forward_mapper.protein_offset_to_span(1) == (3, 5)
        # cds_length 20 is not a multiple of 3 -- the final codon (offset 18-20) is clipped to 19.
        assert forward_mapper.protein_offset_to_span(6) == (18, 19)
        assert forward_mapper.protein_offset_to_span(7) is None
        assert forward_mapper.protein_offset_to_span(-1) is None


class TestAccelerator:
    """Test indexed coordinate mapping."""

    def test_accelerator_matches_bisect_at_every_position(
        self, forward_mapper: Mapper, reverse_mapper: Mapper
    ) -> None:
        """Test that indexed and non-indexed coordinate mappings agree."""
        for mapper in (forward_mapper, reverse_mapper):
            accelerator = mapper.build_index()
            for offset in range(20):
                genomic_position = mapper.offset_to_genomic(offset)
                indexed_genomic_position = mapper.offset_to_genomic(offset, accelerator=accelerator)
                assert genomic_position == indexed_genomic_position
            for position in range(95, 215):
                sequence_offset = mapper.genomic_to_offset(position)
                indexed_sequence_offset = mapper.genomic_to_offset(position, accelerator=accelerator)
                assert sequence_offset == indexed_sequence_offset
