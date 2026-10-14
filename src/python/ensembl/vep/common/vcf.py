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
"""How VCF spells a set of alleles: the trimming VEP applies, and the padding base an indel needs."""

from __future__ import annotations

from collections.abc import Callable, Sequence

from ensembl.vep.common.constants import EMPTY_ALLELE

__all__ = ["trim_right", "vcf_record_alleles"]

_BASES = frozenset("ACGTN")
_PLAIN = frozenset("ACGTN.")
"""What a plain, non-symbolic allele may contain. ``.`` counts, as it does in ensembl-variation."""


def trim_right(alleles: Sequence[str]) -> list[str]:
    """Remove the bases every allele ends with, never trimming an allele away entirely.

    A port of ensembl-variation's ``Utils::Sequence::trim_right``: trimming stops as soon as any
    allele is down to one base, or the alleles no longer share their last base.
    """
    trimmed = list(alleles)
    if len(trimmed) < 2:
        return trimmed
    while all(len(allele) > 1 for allele in trimmed) and len({allele[-1] for allele in trimmed}) == 1:
        trimmed = [allele[:-1] for allele in trimmed]
    return trimmed


def _is_bases(allele: str) -> bool:
    """Return whether an allele is empty or plain bases, as opposed to symbolic."""
    return allele == "" or set(allele) <= _BASES


def vcf_record_alleles(
    start: int, alleles: Sequence[str], base_at: Callable[[int], str]
) -> tuple[int, list[str]]:
    """Return the ``POS`` and the alleles a VCF record spells a short variant with.

    A port of ensembl-variation's ``VariationFeature::to_VCF_record``. Alleles of one length --
    SNVs and MNVs -- are written as they are. Otherwise the common trailing bases are trimmed, and a
    padding base is added unless every allele already starts with the same base. A symbolic allele
    is wrapped in angle brackets.

    One deliberate difference: at position 1 there is no preceding base, and Perl takes whatever the
    slice returns there. The VCF specification pads with the *following* base instead, and so does
    this.

    Args:
        start (int): 1-based start, in VEP's trimmed convention -- an insertion sits before
            ``start``.
        alleles (Sequence[str]): The reference first, then the alternates; ``-`` or an empty string
            for an empty allele.
        base_at (Callable[[int], str]): Returns the reference base at a 1-based position, and is
            only called when a padding base is needed.

    Returns:
        tuple[int, list[str]]: The record's position, and its alleles with the reference first.
    """
    cleaned = [allele.replace(EMPTY_ALLELE, "") for allele in alleles]
    symbolic = any(allele and not set(allele) <= _PLAIN for allele in cleaned)
    if not symbolic and len({len(allele) for allele in cleaned}) == 1:
        return start, cleaned

    if not symbolic:
        cleaned = trim_right(cleaned)
    if len(cleaned) == 1:
        cleaned.insert(0, "")

    first_bases = {allele[:1] for allele in cleaned if "*" not in allele}
    position = start
    before = ""
    after = ""
    if len(first_bases) > 1:
        if start > 1:
            position = start - 1
            before = base_at(position)
        else:
            after = base_at(start + len(cleaned[0]))

    spelled = []
    for allele in cleaned:
        spelled.append(before + allele + after if _is_bases(allele) else f"<{allele}>")
    return position, spelled
