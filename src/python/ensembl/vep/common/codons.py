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

try:
    from ensembl.vep.ensembl_vep_core import translate_dna as _translate_dna_native
except ImportError:
    _translate_dna_native = None

__all__ = [
    "CODON_TABLE",
    "DEFAULT_GENETIC_CODE_TABLE_ID",
    "GENETIC_CODE_TABLE_IDS",
    "STOP_CODONS",
    "get_stop_codons",
    "translate_dna",
]

_CODON_BASES = "TCAG"
DEFAULT_GENETIC_CODE_TABLE_ID = 1
"""Default NCBI genetic-code table identifier."""

_AMINO_ACIDS_BY_TABLE_ID = {
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    2: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    3: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    5: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    6: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    9: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    10: "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    12: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    13: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    14: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    15: "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    16: "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    21: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    22: "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    23: "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    24: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
    25: "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    26: "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    27: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    28: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    29: "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    30: "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    31: "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    32: "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    33: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
}


def _build_codon_table(table_id: int = DEFAULT_GENETIC_CODE_TABLE_ID) -> dict[str, str]:
    """Build an NCBI genetic-code table.

    This uses the NCBI ``ncbieaa`` row. Alternative initiation codons belong to the separate
    ``sncbieaa`` row and do not change translation within a coding sequence.

    Args:
        table_id (int): NCBI genetic-code table identifier.

    Returns:
        dict[str, str]: Codons mapped to one-letter amino acids or ``*`` for stop codons.

    Raises:
        ValueError: If ``table_id`` is not an NCBI genetic-code table identifier.
    """
    try:
        amino_acids = _AMINO_ACIDS_BY_TABLE_ID[table_id]
    except KeyError as error:
        raise ValueError(f"Unsupported NCBI genetic-code table ID: {table_id}") from error

    codons = ("".join(bases) for bases in product(_CODON_BASES, repeat=3))
    return dict(zip(codons, amino_acids, strict=True))


GENETIC_CODE_TABLE_IDS: frozenset[int] = frozenset(_AMINO_ACIDS_BY_TABLE_ID)
"""Supported NCBI genetic-code table identifiers."""

_CODON_TABLES = {table_id: _build_codon_table(table_id) for table_id in GENETIC_CODE_TABLE_IDS}
_STOP_CODONS_BY_TABLE_ID = {
    table_id: frozenset(codon for codon, amino_acid in codon_table.items() if amino_acid == "*")
    for table_id, codon_table in _CODON_TABLES.items()
}


CODON_TABLE: dict[str, str] = _CODON_TABLES[DEFAULT_GENETIC_CODE_TABLE_ID]
"""Standard nuclear codons mapped to one-letter amino acids or ``*`` for stop codons."""

STOP_CODONS: frozenset[str] = _STOP_CODONS_BY_TABLE_ID[DEFAULT_GENETIC_CODE_TABLE_ID]
"""Standard nuclear stop codons."""


_NATIVE_SIGNATURE = getattr(_translate_dna_native, "__text_signature__", "") or ""
_NATIVE_SUPPORTS_TABLE_ID = "table_id" in _NATIVE_SIGNATURE


def _codon_table(table_id: int) -> dict[str, str]:
    """Return a supported codon table."""
    try:
        return _CODON_TABLES[table_id]
    except KeyError as error:
        raise ValueError(f"Unsupported NCBI genetic-code table ID: {table_id}") from error


def get_stop_codons(
    table_id: int = DEFAULT_GENETIC_CODE_TABLE_ID,
) -> frozenset[str]:
    """Return stop codons for an NCBI genetic-code table."""
    _codon_table(table_id)
    return _STOP_CODONS_BY_TABLE_ID[table_id]


def _translate_dna_python(coding_sequence: str, codon_table: dict[str, str]) -> str:
    """Translate DNA using a resolved Python codon table."""
    return "".join(
        codon_table.get(coding_sequence[position : position + 3], "X")
        for position in range(0, len(coding_sequence), 3)
    )


def translate_dna(
    coding_sequence: str,
    table_id: int = DEFAULT_GENETIC_CODE_TABLE_ID,
) -> str:
    """Translate a coding-orientation DNA sequence.

    Args:
        coding_sequence (str): Uppercase DNA sequence in transcript-sense orientation.
        table_id (int): NCBI genetic-code table identifier.

    Returns:
        str: One-letter amino acid sequence, using ``*`` for stop codons and ``X`` for unknown
            or incomplete codons.

    Raises:
        ValueError: If ``table_id`` is not an NCBI genetic-code table identifier.
    """
    codon_table = _codon_table(table_id)

    if _translate_dna_native is not None:
        if table_id == DEFAULT_GENETIC_CODE_TABLE_ID:
            return _translate_dna_native(coding_sequence)
        if _NATIVE_SUPPORTS_TABLE_ID:
            return _translate_dna_native(coding_sequence, table_id)
    return _translate_dna_python(coding_sequence, codon_table)
