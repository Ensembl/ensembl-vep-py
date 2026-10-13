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
"""Combine corpora into one consumer-facing corpus."""

from __future__ import annotations

from collections.abc import Sequence

from attrs import define, field

from ensembl.vep.corpus.base import Corpus, NamedOverlap, NamedVariant
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.corpus.validation import validate_corpus
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.sequence import Slice

__all__ = ["CombinedCorpus", "combine_corpora"]


def _corpus_tuple(corpora: Sequence[Corpus]) -> tuple[Corpus, ...]:
    """Freeze component corpora for stable traversal."""
    return tuple(corpora)


@define(kw_only=True, auto_attribs=True)
class CombinedCorpus(Corpus):
    """Present several corpora through one flattened corpus view.

    Attributes:
        spec (CorpusSpec): Specification describing the combined construction.
        corpora (tuple[Corpus, ...]): Component corpora in traversal and cleanup order.

    Sequence regions, variants, transcripts, and overlaps are flattened and validated when the
    combined corpus is constructed.
    """

    _spec: CorpusSpec = field(alias="spec")
    corpora: tuple[Corpus, ...] = field(converter=_corpus_tuple)
    _sequence_regions: tuple[Slice, ...] = field(init=False)
    _named_variants: tuple[NamedVariant, ...] = field(init=False)
    _transcripts: tuple[Transcript, ...] = field(init=False)
    _named_overlaps: tuple[NamedOverlap, ...] = field(init=False)

    def __attrs_post_init__(self) -> None:
        """Flatten and validate component membership."""
        if not self.corpora:
            raise ValueError("a combined corpus needs at least one component corpus")
        self._sequence_regions = tuple(
            sequence_region for corpus in self.corpora for sequence_region in corpus.sequence_regions
        )
        self._named_variants = tuple(
            named_variant for corpus in self.corpora for named_variant in corpus.named_variants
        )
        self._transcripts = tuple(transcript for corpus in self.corpora for transcript in corpus.transcripts)
        self._named_overlaps = tuple(
            named_overlap for corpus in self.corpora for named_overlap in corpus.named_overlaps
        )
        validate_corpus(
            self.sequence_regions,
            self.named_variants,
            self.transcripts,
            self.named_overlaps,
        )

    @property
    def spec(self) -> CorpusSpec:
        """Return the combined construction specification."""
        return self._spec

    @property
    def sequence_regions(self) -> Sequence[Slice]:
        """Return every component sequence region."""
        return self._sequence_regions

    @property
    def named_variants(self) -> Sequence[NamedVariant]:
        """Return every component variant."""
        return self._named_variants

    @property
    def transcripts(self) -> Sequence[Transcript]:
        """Return every component transcript."""
        return self._transcripts

    @property
    def named_overlaps(self) -> Sequence[NamedOverlap]:
        """Return every component transcript overlap."""
        return self._named_overlaps

    def close(self) -> None:
        """Release every component corpus."""
        for corpus in self.corpora:
            corpus.close()


def combine_corpora(specification: CorpusSpec, corpora: Sequence[Corpus]) -> CombinedCorpus:
    """Return one validated corpus over the supplied component corpora.

    Args:
        specification (CorpusSpec): Specification identifying the combined corpus.
        corpora (Sequence[Corpus]): Components in traversal and cleanup order.

    Returns:
        CombinedCorpus: Flattened view of all component corpora.

    Raises:
        ValueError: If no components were supplied or their memberships conflict.
    """
    return CombinedCorpus(spec=specification, corpora=corpora)
