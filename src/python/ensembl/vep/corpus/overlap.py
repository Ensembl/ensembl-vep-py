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
"""Construct corpus transcript overlaps."""

from __future__ import annotations

from collections.abc import Sequence

from ensembl.vep.corpus.base import NamedOverlap, NamedVariant
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele

__all__ = ["build_transcript_overlaps"]


def build_transcript_overlaps(
    named_variants: Sequence[NamedVariant],
    transcripts: Sequence[Transcript],
    *,
    window: int = 5000,
) -> tuple[NamedOverlap, ...]:
    """Return same-region variant-transcript pairs within ``window`` bases.

    Args:
        named_variants (Sequence[NamedVariant]): Named alleles in stable corpus order.
        transcripts (Sequence[Transcript]): Candidate transcripts in stable corpus order.
        window (int): Flanking bases included on either side of each transcript.

    Returns:
        tuple[NamedOverlap, ...]: Matching pairs in variant then transcript order.

    Raises:
        ValueError: If ``window`` is negative.
    """
    if window < 0:
        raise ValueError("window must be zero or greater")

    named_overlaps = []
    for name, variant_allele in named_variants:
        if variant_allele.slice is None:
            continue
        variant_min = min(variant_allele.start, variant_allele.end)
        variant_max = max(variant_allele.start, variant_allele.end)
        for transcript in transcripts:
            if transcript.slice is not variant_allele.slice:
                continue
            if variant_max < transcript.start - window or variant_min > transcript.end + window:
                continue
            named_overlaps.append(
                (
                    name,
                    TranscriptVariationAllele(
                        feature=transcript,
                        overlapping_feature=variant_allele,
                    ),
                )
            )
    return tuple(named_overlaps)
