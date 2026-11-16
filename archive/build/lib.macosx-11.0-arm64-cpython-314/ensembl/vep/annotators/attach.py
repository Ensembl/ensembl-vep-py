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
"""An annotator's answers onto the model: consequences onto overlaps, then overlaps onto their allele.

The order is load-bearing. ``VariationAllele.transcript_variation_alleles`` carries
``on_setattr=invalidates("_overlap_consequences")``, and that aggregate is built from the overlaps'
own consequences, lazily. Assign every overlap its consequences first and attach the overlaps
second, and the aggregate is always computed from complete overlaps. The other way round, anything
that read the aggregate in between would keep a stale, empty answer cached.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

from ensembl.vep.annotators.base import Annotator
from ensembl.vep.models.overlap import OverlapConsequence, TranscriptVariationAllele
from ensembl.vep.models.variation import Variation, VariationAllele

__all__ = ["annotate_and_attach", "annotate_variations"]


def annotate_and_attach(
    annotator: Annotator[OverlapConsequence],
    allele: VariationAllele,
    overlaps: Sequence[TranscriptVariationAllele],
) -> None:
    """Annotate an allele's overlaps in one batch and attach them, consequences first.

    Args:
        annotator (Annotator[OverlapConsequence]): The annotator to run.
        allele (VariationAllele): The allele the overlaps belong to. An allele with no overlaps is
            attached an empty list, which is what makes it intergenic downstream.
        overlaps (Sequence[TranscriptVariationAllele]): The allele's overlaps, each already built
            against ``allele``.
    """
    batch = list(overlaps)
    results = annotator.annotate_feature_overlaps(batch) if batch else []
    for overlap, consequences in zip(batch, results, strict=True):
        overlap.overlap_consequences = list(consequences)
    allele.transcript_variation_alleles = batch


def annotate_variations(
    annotator: Annotator[OverlapConsequence],
    named_overlaps: Iterable[tuple[str, TranscriptVariationAllele]],
    variations: Iterable[Variation] | None = None,
) -> list[Variation]:
    """Annotate a source's overlaps and return the variations they belong to.

    The overlaps of a file-backed source are built against alleles that already belong to real
    ``Variation`` objects, so the variations are recovered rather than rebuilt: a multi-allelic
    record stays one variation with two alleles, which is what an output line has to be.

    Every allele of every variation is attached a list, including the alleles that overlap nothing.
    That emptiness is the signal for intergenic downstream, so an allele left unattached would be
    silently dropped instead of reported.

    Pass ``variations`` whenever the full set is known. A variant that overlaps no feature at all
    produces no overlaps, so it cannot be recovered from them -- and dropping it would quietly
    shorten the output, losing exactly the variants that should be reported as intergenic.

    Args:
        annotator (Annotator[OverlapConsequence]): The annotator to run.
        named_overlaps (Iterable[tuple[str, TranscriptVariationAllele]]): The source's overlaps,
            each paired with its variant's name, which this ignores.
        variations (Iterable[Variation] | None): Every variation the input holds, in input order.
            Defaults to only those the overlaps reach.

    Returns:
        list[Variation]: Each variation once, annotated. In input order when ``variations`` was
            given, otherwise in the order each variation's first overlap appeared.
    """
    by_allele: dict[int, list[TranscriptVariationAllele]] = {}
    alleles: dict[int, VariationAllele] = {}
    found: dict[int, Variation] = {}
    for _name, overlap in named_overlaps:
        allele = overlap.overlapping_feature
        assert isinstance(allele, VariationAllele)
        by_allele.setdefault(id(allele), []).append(overlap)
        alleles.setdefault(id(allele), allele)
        variation = allele.variation
        if variation is not None:
            found.setdefault(id(variation), variation)

    universe = list(found.values()) if variations is None else list(variations)
    ordered: dict[int, Variation] = {}
    for variation in universe:
        ordered.setdefault(id(variation), variation)
    for key, variation in found.items():
        ordered.setdefault(key, variation)

    attached: set[int] = set()
    for variation in ordered.values():
        for allele in variation.variation_alleles or ():
            annotate_and_attach(annotator, allele, by_allele.get(id(allele), []))
            attached.add(id(allele))

    # An allele whose variation never resolved is still worth annotating; it just cannot be
    # reported as part of one.
    for key, allele in alleles.items():
        if key not in attached:
            annotate_and_attach(annotator, allele, by_allele[key])

    return list(ordered.values())
