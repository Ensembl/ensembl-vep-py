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
"""Codon table and DNA translation utilities."""

from __future__ import annotations

from itertools import product

__all__ = ["CODON_TABLE", "STOP_CODONS", "translate_dna"]

_CODON_BASES = "TCAG"
_AMINO_ACIDS_BY_CODON_ORDER = (
    "FFLLSSSSYY**CC*W"
    + "LLLLPPPPHHQQRRRR"
    + "IIIMTTTTNNKKSSRR"
    + "VVVVAAAADDEEGGGG"
)


def _build_codon_table() -> dict[str, str]:
    """Build the standard nuclear codon table.

    Returns:
        dict[str, str]: Codons mapped to one-letter amino acids or ``*`` for stop codons.
    """
    codons = ("".join(bases) for bases in product(_CODON_BASES, repeat=3))
    return dict(zip(codons, _AMINO_ACIDS_BY_CODON_ORDER, strict=True))


CODON_TABLE: dict[str, str] = _build_codon_table()
"""Standard nuclear codons mapped to one-letter amino acids or ``*`` for stop codons."""

STOP_CODONS: frozenset[str] = frozenset(
    codon for codon, amino_acid in CODON_TABLE.items() if amino_acid == "*"
)
"""Standard nuclear stop codons."""


def translate_dna(coding_sequence: str) -> str:
    """Translate a coding-orientation DNA sequence.

    Args:
        coding_sequence (str): Uppercase DNA sequence in transcript-sense orientation.

    Returns:
        str: One-letter amino acid sequence, using ``*`` for stop codons and ``X`` for unknown
            or incomplete codons.
    """
    return "".join(
        CODON_TABLE.get(coding_sequence[position : position + 3], "X")
        for position in range(0, len(coding_sequence), 3)
    )
