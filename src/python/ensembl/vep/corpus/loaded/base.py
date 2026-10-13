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
"""A corpus materialised from external sources."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from contextlib import ExitStack
from typing import Protocol, runtime_checkable

from attrs import define, field

from ensembl.vep.corpus.base import Corpus, NamedOverlap, NamedVariant
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.corpus.validation import validate_corpus
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.sequence import Slice
from ensembl.vep.sources.sequence_source import SequenceSource

__all__ = ["CorpusResource", "LoadedCorpus"]


@runtime_checkable
class CorpusResource(Protocol):
    """A resource retained by a materialised corpus and released with it."""

    def close(self) -> None:
        """Release the resource."""


def _slice_tuple(values: Iterable[Slice]) -> tuple[Slice, ...]:
    """Freeze sequence regions in source order."""
    return tuple(values)


def _variant_tuple(values: Iterable[NamedVariant]) -> tuple[NamedVariant, ...]:
    """Freeze named variants in source order."""
    return tuple(values)


def _transcript_tuple(values: Iterable[Transcript]) -> tuple[Transcript, ...]:
    """Freeze transcripts in source order."""
    return tuple(values)


def _overlap_tuple(values: Iterable[NamedOverlap]) -> tuple[NamedOverlap, ...]:
    """Freeze named overlaps in source order."""
    return tuple(values)


def _resource_tuple(values: Iterable[CorpusResource]) -> tuple[CorpusResource, ...]:
    """Freeze owned resources in acquisition order."""
    return tuple(values)


@define(kw_only=True, auto_attribs=True)
class LoadedCorpus(Corpus):
    """Hold a validated corpus and any external resources that back it.

    Attributes:
        spec (CorpusSpec): Construction specification.
        sequence_regions (Sequence[Slice]): Referenced sequence regions.
        named_variants (Sequence[NamedVariant]): Stable variant names and alleles.
        transcripts (Sequence[Transcript]): Loaded transcripts.
        named_overlaps (Sequence[NamedOverlap]): Variant-transcript overlaps.
        resources (Sequence[CorpusResource]): Lazy sources or connections owned by the corpus.
    """

    _spec: CorpusSpec = field(alias="spec")
    _sequence_regions: tuple[Slice, ...] = field(converter=_slice_tuple, alias="sequence_regions")
    _named_variants: tuple[NamedVariant, ...] = field(converter=_variant_tuple, alias="named_variants")
    _transcripts: tuple[Transcript, ...] = field(converter=_transcript_tuple, alias="transcripts")
    _named_overlaps: tuple[NamedOverlap, ...] = field(converter=_overlap_tuple, alias="named_overlaps")
    resources: tuple[CorpusResource, ...] = field(factory=tuple, converter=_resource_tuple, repr=False)
    _closed: bool = field(default=False, init=False, repr=False)

    def __attrs_post_init__(self) -> None:
        """Validate relationships after freezing input collections."""
        validate_corpus(
            self.sequence_regions,
            self.named_variants,
            self.transcripts,
            self.named_overlaps,
        )

    @property
    def spec(self) -> CorpusSpec:
        """Return the specification that loaded the corpus."""
        return self._spec

    @property
    def sequence_regions(self) -> Sequence[Slice]:
        """Return every referenced sequence region."""
        return self._sequence_regions

    @property
    def named_variants(self) -> Sequence[NamedVariant]:
        """Return every named variation allele."""
        return self._named_variants

    @property
    def transcripts(self) -> Sequence[Transcript]:
        """Return every loaded transcript."""
        return self._transcripts

    @property
    def named_overlaps(self) -> Sequence[NamedOverlap]:
        """Return every named transcript overlap."""
        return self._named_overlaps

    @property
    def description(self) -> str:
        """Return a short loaded-corpus description."""
        return f"loaded corpus ({self.spec.kind})"

    @property
    def sequence_source(self) -> SequenceSource | None:
        """Return the retained sequence source, when this loader uses one."""
        return next(
            (resource for resource in self.resources if isinstance(resource, SequenceSource)),
            None,
        )

    def close(self) -> None:
        """Release every owned source or connection once."""
        if self._closed:
            return
        self._closed = True
        with ExitStack() as resources:
            for resource in self.resources:
                resources.callback(resource.close)


def collect_sequence_regions(features: Iterable[object]) -> tuple[Slice, ...]:
    """Return referenced slices once each, preserving encounter order."""
    sequence_regions = []
    identities: set[int] = set()
    for feature in features:
        sequence_region = getattr(feature, "slice", None)
        if sequence_region is None or id(sequence_region) in identities:
            continue
        identities.add(id(sequence_region))
        sequence_regions.append(sequence_region)
    return tuple(sequence_regions)
