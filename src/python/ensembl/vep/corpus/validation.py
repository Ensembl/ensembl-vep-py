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
"""Corpus consistency validation."""

from __future__ import annotations

from collections.abc import Sequence

from ensembl.vep.corpus.base import NamedOverlap, NamedVariant
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.sequence import Slice

__all__ = ["validate_corpus"]


def validate_corpus(
    sequence_regions: Sequence[Slice],
    named_variants: Sequence[NamedVariant],
    transcripts: Sequence[Transcript],
    named_overlaps: Sequence[NamedOverlap],
) -> None:
    """Reject inconsistent corpus membership and relationships."""
    if not sequence_regions:
        raise ValueError("a corpus must contain at least one sequence region")
    sequence_region_identities = {id(sequence_region) for sequence_region in sequence_regions}
    if len(sequence_region_identities) != len(sequence_regions):
        raise ValueError("duplicate corpus sequence region object")
    sequence_region_names = [sequence_region.seq_region_name for sequence_region in sequence_regions]
    if len(set(sequence_region_names)) != len(sequence_region_names):
        raise ValueError("duplicate corpus sequence region name")

    variants_by_name = {}
    for name, allele in named_variants:
        if name in variants_by_name:
            raise ValueError(f"duplicate corpus variant name {name!r}")
        if allele.slice is None or id(allele.slice) not in sequence_region_identities:
            raise ValueError(f"corpus variant {name!r} belongs to another sequence region")
        variants_by_name[name] = allele

    transcript_identities = _validate_transcripts(transcripts, sequence_region_identities)
    _validate_overlaps(named_overlaps, variants_by_name, transcript_identities)


def _validate_transcripts(
    transcripts: Sequence[Transcript], sequence_region_identities: set[int]
) -> set[int]:
    transcript_identities: set[int] = set()
    transcript_stable_ids: set[str] = set()
    for transcript in transcripts:
        if id(transcript) in transcript_identities:
            raise ValueError(f"duplicate corpus transcript object {transcript.stable_id!r}")
        if transcript.stable_id in transcript_stable_ids:
            raise ValueError(f"duplicate corpus transcript stable id {transcript.stable_id!r}")
        if transcript.slice is None or id(transcript.slice) not in sequence_region_identities:
            raise ValueError(f"corpus transcript {transcript.stable_id!r} belongs to another sequence region")
        transcript_identities.add(id(transcript))
        transcript_stable_ids.add(transcript.stable_id)

    return transcript_identities


def _validate_overlaps(
    named_overlaps: Sequence[NamedOverlap], variants_by_name: dict, transcript_identities: set[int]
) -> None:
    overlap_identities: set[tuple[str, int]] = set()
    for name, overlap in named_overlaps:
        defined_allele = variants_by_name.get(name)
        if defined_allele is None:
            raise ValueError(f"overlap {name!r} has no corpus variant")
        if overlap.overlapping_feature is not defined_allele:
            raise ValueError(f"overlap {name!r} does not refer to its corpus variant allele")
        if id(overlap.feature) not in transcript_identities:
            raise ValueError(f"overlap {name!r} refers to a transcript outside the corpus")
        if overlap.feature.slice is not overlap.overlapping_feature.slice:
            raise ValueError(
                f"overlap {name!r} refers to a transcript and variant on different sequence regions"
            )
        overlap_identity = (name, id(overlap.feature))
        if overlap_identity in overlap_identities:
            raise ValueError(
                f"duplicate overlap for variant {name!r} and transcript {overlap.feature.stable_id!r}"
            )
        overlap_identities.add(overlap_identity)
