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
"""The corpus contract."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence
from typing import Protocol, TypeAlias

from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.sequence import Slice
from ensembl.vep.models.variation import VariationAllele

__all__ = ["Corpus", "CorpusLike", "CorpusProvider", "NamedOverlap", "NamedVariant", "as_corpus"]

NamedVariant = tuple[str, VariationAllele]
"""A stable corpus name and its variation allele."""

NamedOverlap = tuple[str, TranscriptVariationAllele]
"""A corpus variant name and one transcript overlap for its allele."""


class Corpus(ABC):
    """Provide inputs for transcript consequence annotation across sequence regions.

    Implementations may generate or load their data from any backend. The biological contract is
    intentionally transcript-oriented for the annotators available today; adding regulatory or
    motif annotation will require a corresponding feature-level extension to this interface.
    """

    @property
    @abstractmethod
    def spec(self) -> CorpusSpec:
        """Return the specification that constructed the corpus."""

    @property
    @abstractmethod
    def sequence_regions(self) -> Sequence[Slice]:
        """Return every sequence region referenced by the corpus."""

    @property
    @abstractmethod
    def named_variants(self) -> Sequence[NamedVariant]:
        """Return every variant in stable corpus order."""

    @property
    @abstractmethod
    def transcripts(self) -> Sequence[Transcript]:
        """Return every transcript in stable corpus order."""

    @property
    @abstractmethod
    def named_overlaps(self) -> Sequence[NamedOverlap]:
        """Return every named transcript overlap in stable corpus order."""

    @property
    def description(self) -> str:
        """Return a short description of the corpus."""
        return f"corpus ({self.spec.kind})"

    def close(self) -> None:
        """Release resources owned by the corpus."""

    def __enter__(self) -> Corpus:
        """Return the corpus for use as a context manager."""
        return self

    def __exit__(self, *_exc_info: object) -> None:
        """Release resources when leaving a context manager."""
        self.close()


class CorpusProvider(Protocol):
    """Compatibility object that exposes a consumer-facing corpus."""

    @property
    def corpus(self) -> Corpus:
        """Return the consumer-facing corpus."""


CorpusLike: TypeAlias = Corpus | CorpusProvider
"""A corpus or a compatibility object that provides one."""


def as_corpus(value: object) -> Corpus:
    """Return ``value`` as a corpus, accepting the retired synthetic facade."""
    if isinstance(value, Corpus):
        return value
    adapted = getattr(value, "corpus", None)
    if isinstance(adapted, Corpus):
        return adapted
    raise TypeError(f"expected a Corpus, got {type(value).__name__}")
