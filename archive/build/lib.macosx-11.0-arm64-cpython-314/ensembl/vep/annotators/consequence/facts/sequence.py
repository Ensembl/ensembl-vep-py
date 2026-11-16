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
"""Sequence facts: arithmetic or string operations over the codon/peptide window, or allele
string lengths. The expensive tier -- ``_codon_alleles``/``_peptide_alleles`` splice and
translate real sequence, and everything else here reads their output or does length arithmetic.
"""

from __future__ import annotations


from ensembl.vep.annotators.consequence.facts.base import fact, get_fact
from ensembl.vep.common import CostTier
from ensembl.vep.annotators.consequence.facts.helper import (
    get_codon_alleles,
    get_peptide_alleles,
    get_full_cds_peptide_alleles,
)
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.models.overlap import TranscriptVariationAllele

__all__ = [
    "alt_matches_ref_codon_edge",
    "alt_peptide_has_ter",
    "alt_peptide_matches_ref_edge",
    "alt_peptide_is_ter",
    "cds_length_change_inframe",
    "cds_start_lost",
    "cds_stop_gain",
    "cds_start_retained",
    "cds_stop_lost",
    "cds_stop_retained",
    "insertion_after_start",
    "insertion_before_start",
    "peptide_changed",
    "peptide_has_X",
    "peptide_same_length",
    "ref_peptide_has_ter",
    "ref_peptide_is_ter",
    "ref_peptide_startswith_ter",
]

EMPTY_ALLELE = "-"
"""How an empty codon window can be spelled, besides the empty string."""


def _is_empty_codon(codon: str) -> bool:
    """Whether a codon window holds no bases at all."""
    return codon in ("", EMPTY_ALLELE)


def _matches_codon_edge(container: str, contained: str) -> bool:
    r"""Whether ``contained`` is the start of ``container``, its end, or whole codons inside it.

    VEP's two-branch test (``VariationEffect::inframe_deletion``): the string match
    ``$container =~ /^\Q$contained\E/ || /\Q$contained\E$/``, then, failing that,
    ``trim_sequences`` on the pair and nothing left of ``contained`` with what remains of
    ``container`` a whole number of codons. An empty ``contained`` matches either edge already,
    as the Perl regex does.
    """
    if _is_empty_codon(contained):
        return True
    if container.startswith(contained) or container.endswith(contained):
        return True

    leading = 0
    while leading < len(container) and leading < len(contained) and container[leading] == contained[leading]:
        leading += 1
    remaining_container, remaining_contained = container[leading:], contained[leading:]

    trailing = 0
    while (
        trailing < len(remaining_container)
        and trailing < len(remaining_contained)
        and remaining_container[-1 - trailing] == remaining_contained[-1 - trailing]
    ):
        trailing += 1
    if trailing:
        remaining_container = remaining_container[:-trailing]
        remaining_contained = remaining_contained[:-trailing]
    return not remaining_contained and len(remaining_container) % 3 == 0


# --------------------------------------------------------------------------
# Comparisons between only codon/peptide alleles variant overlaps.
# --------------------------------------------------------------------------


@fact(CostTier.SEQUENCE)
def _codon_alleles(fact_context: FactContext[TranscriptVariationAllele]) -> tuple[str, str] | None:
    """(reference_codon, alternate_codon) spanning whole codons around the variant, upper case.

    The overlap cases its codons for display; every comparison here is on the bases alone.
    """
    reference_codon = fact_context.overlap.reference_codon
    alternate_codon = fact_context.overlap.alternate_codon
    if reference_codon is None or alternate_codon is None:
        return None
    return reference_codon.upper(), alternate_codon.upper()


@fact(CostTier.SEQUENCE)
def _peptide_alleles(fact_context: FactContext[TranscriptVariationAllele]) -> tuple[str, str] | None:
    """(reference_peptide, alternate_peptide) translated from the codon window."""
    reference_peptide = fact_context.overlap.reference_peptide
    alternate_peptide = fact_context.overlap.alternate_peptide
    if reference_peptide is None or alternate_peptide is None:
        return None
    return reference_peptide, alternate_peptide


@fact(
    CostTier.SEQUENCE,
    deps=("_codon_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
    ),
)
def alt_matches_ref_codon_edge(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """What the deletion removed is whole codons from one edge of the window, or from inside it.

    VEP's own test for an inframe deletion (``VariationEffect::inframe_deletion``): see
    ``_matches_codon_edge`` for its two branches.
    """
    codons = get_codon_alleles(fact_context)
    if codons is None:
        return False
    reference_codon, alternate_codon = codons
    return _matches_codon_edge(reference_codon, alternate_codon)


@fact(
    CostTier.SEQUENCE,
    deps=("_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
    ),
)
def alt_peptide_matches_ref_edge(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The alternate peptide window still opens or closes with the whole reference one.

    VEP's test for an inframe insertion (``VariationEffect::inframe_insertion``), which it makes on
    the peptides rather than the codons: residues added at one end of the window leave the
    reference residues intact at the other. The alternate window is cut at its first stop first, so
    that an insertion which also introduces a stop is judged on what precedes it.
    """
    peptides = get_peptide_alleles(fact_context)
    if peptides is None:
        return False
    reference_peptide, alternate_peptide = peptides
    if not reference_peptide:
        return True
    stop_index = alternate_peptide.find("*")
    if stop_index != -1:
        alternate_peptide = alternate_peptide[: stop_index + 1]
    return alternate_peptide.startswith(reference_peptide) or alternate_peptide.endswith(reference_peptide)


@fact(
    CostTier.SEQUENCE,
    deps=("_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
    ),
)
def peptide_changed(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The reference and alternate peptide windows differ."""
    peptides = get_peptide_alleles(fact_context)
    return peptides is not None and peptides[0] != peptides[1]


@fact(
    CostTier.SEQUENCE,
    deps=("_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
    ),
)
def peptide_same_length(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The reference and alternate peptide windows have the same length."""
    peptides = get_peptide_alleles(fact_context)
    return peptides is not None and len(peptides[0]) == len(peptides[1])


@fact(
    CostTier.SEQUENCE,
    deps=("_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
    ),
)
def alt_peptide_has_ter(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The alternate peptide window contains a stop codon."""
    peptides = get_peptide_alleles(fact_context)
    return peptides is not None and "*" in peptides[1]


@fact(
    CostTier.SEQUENCE,
    deps=("_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
    ),
)
def ref_peptide_has_ter(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The reference peptide window contains a stop codon."""
    peptides = get_peptide_alleles(fact_context)
    return peptides is not None and "*" in peptides[0]


@fact(
    CostTier.SEQUENCE,
    deps=("_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
    ),
)
def peptide_has_X(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """Either peptide window contains an ambiguous residue ("X").

    "X" comes from an unresolvable codon (incomplete, or an ambiguity code).
    """
    peptides = get_peptide_alleles(fact_context)
    return peptides is not None and ("X" in peptides[0] or "X" in peptides[1])


@fact(
    CostTier.SEQUENCE,
    deps=("_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
    ),
)
def ref_peptide_is_ter(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The reference peptide window is exactly a single stop codon, no more, no less."""
    peptides = get_peptide_alleles(fact_context)
    return peptides is not None and peptides[0] == "*"


@fact(
    CostTier.SEQUENCE,
    deps=("_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
    ),
)
def alt_peptide_is_ter(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The alternate peptide window is exactly a single stop codon, no more, no less."""
    peptides = get_peptide_alleles(fact_context)
    return peptides is not None and peptides[1] == "*"


@fact(
    CostTier.SEQUENCE,
    deps=("_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
    ),
)
def ref_peptide_startswith_ter(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The reference peptide window begins with a stop codon."""
    peptides = get_peptide_alleles(fact_context)
    return peptides is not None and peptides[0].startswith("*")


@fact(
    CostTier.SEQUENCE,
    deps=("_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
    ),
)
def insertion_before_start(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The insertion happens before the start codon keeping it intact.

    Because variant overlaps start codon (verified by overlap_start_codon fact requirement) the reference
    peptide should include the start peptide. If the alternate peptide window ends with the reference peptide
    then it also include the start codon with inserted amino acids before it.

    An empty reference peptide is excluded because ``str.endswith("")`` is trivially True for any string,
    which would make this fact vacuously True for every overlap.
    """
    peptides = get_peptide_alleles(fact_context)
    if peptides is None:
        return False
    reference_peptide, alternate_peptide = peptides
    return bool(reference_peptide) and alternate_peptide.endswith(reference_peptide)


@fact(
    CostTier.SEQUENCE,
    deps=("_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
    ),
)
def insertion_after_start(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The insertion happens after the start codon keeping it intact."""
    peptides = get_peptide_alleles(fact_context)
    if peptides is None:
        return False
    reference_peptide, alternate_peptide = peptides
    return bool(reference_peptide) and alternate_peptide.startswith(reference_peptide)


# --------------------------------------------------------------------------
# Comparisons between full CDS/protein sequences.
#
# Unlike _codon_alleles/_peptide_alleles above, a few bases around the variant, the facts below need
# the real downstream sequence. For example, an indel elsewhere in the CDS that shifts the
# reading frame isn't visible in a narrow window.
#
# Deliberately not routed through codon_cache as that cache's key is the narrow codon window's genomic
# positions specifically (see its own module docstring on the "frame-vs-footprint failure family")
# --------------------------------------------------------------------------


@fact(CostTier.TRANSCRIPT_SEQUENCE)
def _full_cds_peptide_alleles(fact_context: FactContext[TranscriptVariationAllele]) -> tuple[str, str] | None:
    """(reference_cds_peptide, alternate_cds_peptide): the whole CDS translated, reference vs.
    the variant applied.

    ``data/fact_costs.json`` (the opt-in measured-cost override, see facts/base.py) has no entry
    for this fact, so its declared ``CostTier.TRANSCRIPT_SEQUENCE`` applies whenever that override
    is enabled too -- not regenerated here since ``vep catalogue costs`` profiles against a live
    corpus and picking one is a judgement call outside this change's scope.
    """
    reference = fact_context.overlap.reference_cds_peptide
    alternate = fact_context.overlap.alternate_cds_peptide
    if reference is None or alternate is None:
        return None
    return reference, alternate


@fact(
    CostTier.FREE,
    deps=("_full_cds_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
        "overlap_start_codon",
    ),
)
def cds_start_lost(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The initiator start codon is no longer encoded."""
    peptides = get_full_cds_peptide_alleles(fact_context)
    if peptides is None:
        return False
    reference_peptide, alternate_peptide = peptides
    # We need slice here to avoid IndexError on empty string
    return reference_peptide[0:1] != alternate_peptide[0:1]


@fact(
    CostTier.FREE,
    deps=("_full_cds_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
    ),
)
def cds_stop_lost(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The stop codon is no longer encoded.

    Requires the alternate peptide to match the reference up to (but not including) the reference's own stop
    position, so a frameshift several codons upstream that happens to remove the stop incidentally doesn't
    get conflated with an edit that actually changes the stop codon itself.
    """
    peptides = get_full_cds_peptide_alleles(fact_context)
    if peptides is None:
        return False
    reference_peptide, alternate_peptide = peptides
    reference_stop = reference_peptide.find("*")
    if reference_stop == -1:
        return False
    if alternate_peptide[:reference_stop] != reference_peptide[:reference_stop]:
        return False
    return len(alternate_peptide) <= reference_stop or alternate_peptide[reference_stop] != "*"


@fact(
    CostTier.FREE,
    deps=("_full_cds_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
        "overlap_start_codon",
    ),
)
def cds_start_retained(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The initiator methionine is still encoded.

    Computed from the peptides directly rather than as ``not cds_start_lost``. Negating a *guarded*
    fact is unsound: when a required atom is false, ``cds_start_lost`` short-circuits to False, and
    ``not False`` claims the start codon is retained for a variant that never touches it -- which
    is exactly what ``verify_requires`` catches by running this body unguarded. Deriving it here
    also returns False (not True) when the peptides are unavailable, matching ``cds_stop_retained``.
    Caching is unaffected: ``_full_cds_peptide_alleles`` is the memoised dep both facts read.
    """

    # Explicit guard here to enforce the design a required False mean the fact short-circuits to False;
    # otherwise, the fact would be True for a variant that never touches the start codon.
    if get_fact(fact_context, "overlap_start_codon") is False:
        return False

    peptides = get_full_cds_peptide_alleles(fact_context)
    if peptides is None:
        return False
    reference_peptide, alternate_peptide = peptides
    # We need slice here to avoid IndexError on empty string
    return reference_peptide[0:1] == alternate_peptide[0:1]


@fact(
    CostTier.FREE,
    deps=("_full_cds_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
        "overlap_stop_codon",
    ),
)
def cds_stop_retained(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The stop codon is unchanged.

    Requires stop at the exact same reading-frame position, with the same prefix leading up to it."""

    # Explicit guard here to enforce the design a required False mean the fact short-circuits to False;
    # otherwise, the fact would be True for a variant that never touches the stop codon. Except for a
    # unknown peptide (generally caused by frameshift) that may cause the stop codon to appear retained.
    if (
        get_fact(fact_context, "overlap_stop_codon") is False
        and get_fact(fact_context, "peptide_has_X") is False
    ):
        return False

    peptides = get_full_cds_peptide_alleles(fact_context)
    if peptides is None:
        return False
    reference_peptide, alternate_peptide = peptides
    reference_stop = reference_peptide.find("*")
    if reference_stop == -1:
        return False
    if alternate_peptide[:reference_stop] != reference_peptide[:reference_stop]:
        return False
    return len(alternate_peptide) > reference_stop and alternate_peptide[reference_stop] == "*"


@fact(
    CostTier.FREE,
    deps=("_full_cds_peptide_alleles",),
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
    ),
)
def cds_stop_gain(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """A new premature stop codon is introduced."""
    peptides = get_full_cds_peptide_alleles(fact_context)
    if peptides is None:
        return False
    reference_peptide, alternate_peptide = peptides
    reference_stop = reference_peptide.find("*")
    if reference_stop == -1:
        return False
    alternate_stop = alternate_peptide.find("*")
    return 0 <= alternate_stop < reference_stop


# --------------------------------------------------------------------------
# Sequence length arithmetic.
# --------------------------------------------------------------------------


@fact(
    CostTier.MAPPING,
    requires=(
        "is_transcript",
        "overlap_feature",
        "overlap_cds",
        "has_translation",
    ),
)
def cds_length_change_inframe(fact_context: FactContext[TranscriptVariationAllele]) -> bool:
    """The variant's net change in coding length is a multiple of three.

    False when the variant does not overlap the CDS. There is then no coding length change to be in
    frame, and reporting True would contradict this fact's own ``requires``.
    """
    length_change = fact_context.overlap.cds_length_change
    if length_change is None:
        return False
    return length_change % 3 == 0
