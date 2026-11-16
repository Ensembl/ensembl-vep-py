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
"""A whole corpus in the shape an annotated file has: one variation per variant, its overlaps attached.

A corpus holds variants plus overlaps against their standalone alleles. A run -- and anything
reading a run's output, like the summary -- works on ``Variation`` objects whose alleles carry
their overlaps. This rebuilds that shape, including variants with no overlaps, so a corpus can be
summarised or written exactly as an annotated file would be.
"""

from __future__ import annotations

from collections.abc import Iterator

from ensembl.vep.annotators.attach import annotate_and_attach
from ensembl.vep.annotators.base import Annotator
from ensembl.vep.corpus.base import Corpus, as_corpus
from ensembl.vep.models.overlap import OverlapConsequence, TranscriptVariationAllele
from ensembl.vep.models.variation import Variation, VariationAllele

__all__ = ["annotated_variations", "variation_for"]


def variation_for(name: str, allele: VariationAllele) -> Variation:
    """Return a single-allele variation with a corpus allele's coordinates and alleles.

    An empty allele stays an empty string, the corpus's own spelling, so the rebuilt allele is the
    same allele the corpus overlap was annotated against.
    """
    return Variation(
        start=allele.start,
        end=allele.end,
        strand=allele.strand,
        slice=allele.slice,
        stable_id=name,
        allele_string=f"{allele.reference_allele or ''}/{allele.alternate_allele or ''}",
    )


def annotated_variations(corpus: Corpus, annotator: Annotator[OverlapConsequence]) -> Iterator[Variation]:
    """Yield one annotated variation per generated variant, in corpus order.

    A variant's overlaps are contiguous in a corpus -- one per feature in its flank window -- so
    grouping by name preserves their order. Variants with no overlaps are still yielded with an
    empty attached-overlap collection.

    Args:
        corpus (Corpus): The corpus to annotate.
        annotator (Annotator[OverlapConsequence]): The annotator to run.

    Yields:
        Variation: Each variant, its allele's overlaps rebuilt against it and annotated.
    """
    corpus = as_corpus(corpus)
    grouped: dict[str, list[TranscriptVariationAllele]] = {}
    for name, overlap in corpus.named_overlaps:
        grouped.setdefault(name, []).append(overlap)
    for name, source_allele in corpus.named_variants:
        variation = variation_for(name, source_allele)
        allele = variation.variation_alleles[0]
        rebuilt = [
            TranscriptVariationAllele(feature=overlap.feature, overlapping_feature=allele)
            for overlap in grouped.get(name, [])
        ]
        annotate_and_attach(annotator, allele, rebuilt)
        yield variation
