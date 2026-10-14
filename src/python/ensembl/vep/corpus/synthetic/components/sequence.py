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
"""Build synthetic sequences."""

from __future__ import annotations

import random

from ensembl.vep.common.constants import CANONICAL_DNA_BASES

__all__ = ["fill_locus_sequence", "random_sequence"]


def random_sequence(length: int, random_generator: random.Random) -> str:
    """Return a random nucleotide sequence.

    Args:
        length (int): Number of bases to generate.
        random_generator (random.Random): Seeded random generator, for reproducibility.

    Returns:
        str: A random sequence of ``ACGT`` bases.
    """
    return "".join(random_generator.choice(CANONICAL_DNA_BASES) for _ in range(length))


def fill_locus_sequence(
    chromosome_chars: list[str],
    start: int,
    end: int,
    random_generator: random.Random,
    *,
    bases: str = "ACG",
) -> None:
    """Fill an inclusive locus with bases drawn from ``bases``.

    The default alphabet cannot produce a stop codon in either strand or reading frame. Passing
    ``ACGT`` produces an unconstrained reference sequence while retaining the same reusable
    coordinate handling.

    Args:
        chromosome_chars (list[str]): Mutable whole-chromosome sequence.
        start (int): One-based inclusive locus start.
        end (int): One-based inclusive locus end.
        random_generator (random.Random): Seeded random generator.
        bases (str): Non-empty alphabet to draw from.

    Raises:
        ValueError: If the locus is outside the chromosome or ``bases`` is empty.
    """
    if start < 1 or end < start or end > len(chromosome_chars):
        raise ValueError(f"locus {start}-{end} is outside a chromosome of {len(chromosome_chars)} bases")
    if not bases:
        raise ValueError("locus sequence alphabet must not be empty")
    for position in range(start, end + 1):
        chromosome_chars[position - 1] = random_generator.choice(bases)
