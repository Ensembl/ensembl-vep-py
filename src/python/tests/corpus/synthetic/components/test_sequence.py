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
"""Tests for synthetic sequence primitives."""

from __future__ import annotations

import random

import pytest

from ensembl.vep.corpus.synthetic.components.sequence import fill_locus_sequence, random_sequence


class TestRandomSequence:
    """Seeded synthetic nucleotide generation."""

    def test_is_seeded_and_contains_only_dna_bases(self) -> None:
        """A seed makes generated references repeatable across test runs."""
        assert random_sequence(20, random.Random(7)) == random_sequence(20, random.Random(7))
        assert set(random_sequence(100, random.Random(11))) <= {"A", "C", "G", "T"}


class TestLocusSequence:
    """Inclusive locus filling."""

    def test_fills_only_the_requested_locus(self) -> None:
        """Coordinates are one-based and the caller controls the nucleotide alphabet."""
        chromosome = list("NNNNNN")
        fill_locus_sequence(chromosome, 2, 5, random.Random(7), bases="AC")
        assert chromosome[0] == chromosome[-1] == "N"
        assert set(chromosome[1:5]) <= {"A", "C"}

    @pytest.mark.parametrize(("start", "end"), [(0, 2), (3, 2), (1, 7)])
    def test_rejects_coordinates_outside_the_chromosome(self, start: int, end: int) -> None:
        """Invalid locus bounds fail before mutating sequence."""
        with pytest.raises(ValueError, match="outside"):
            fill_locus_sequence(list("NNNNNN"), start, end, random.Random(7))

    def test_rejects_an_empty_alphabet(self) -> None:
        """Every filled position needs at least one possible base."""
        with pytest.raises(ValueError, match="alphabet"):
            fill_locus_sequence(list("NNNNNN"), 1, 2, random.Random(7), bases="")
