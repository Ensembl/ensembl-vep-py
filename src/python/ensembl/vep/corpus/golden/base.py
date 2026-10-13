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
"""The golden corpus model."""

from __future__ import annotations

from collections.abc import Sequence
from contextlib import AbstractContextManager

from attrs import define, field

from ensembl.vep.corpus.base import Corpus, NamedOverlap, NamedVariant
from ensembl.vep.corpus.golden.terms import GoldenTerms
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.sequence import Slice

__all__ = ["GoldenCorpus"]


@define(kw_only=True, auto_attribs=True)
class GoldenCorpus(Corpus):
    """Pair a corpus with pinned expected consequence terms.

    Attributes:
        spec (CorpusSpec): Specification identifying the golden dataset.
        corpus (Corpus): Loaded or generated corpus whose overlaps the terms describe.
        expected_terms (GoldenTerms): Pinned consequences for the wrapped corpus.
        resource_context (AbstractContextManager[object] | None): Context retaining packaged
            files for the lifetime of the corpus.
    """

    _spec: CorpusSpec = field(alias="spec")
    corpus: Corpus
    expected_terms: GoldenTerms
    resource_context: AbstractContextManager[object] | None = field(default=None, repr=False)
    _closed: bool = field(default=False, init=False, repr=False)

    @property
    def spec(self) -> CorpusSpec:
        """Return the golden corpus specification."""
        return self._spec

    @property
    def sequence_regions(self) -> Sequence[Slice]:
        """Return the wrapped sequence regions."""
        return self.corpus.sequence_regions

    @property
    def named_variants(self) -> Sequence[NamedVariant]:
        """Return the wrapped variants."""
        return self.corpus.named_variants

    @property
    def transcripts(self) -> Sequence[Transcript]:
        """Return the wrapped transcripts."""
        return self.corpus.transcripts

    @property
    def named_overlaps(self) -> Sequence[NamedOverlap]:
        """Return the wrapped overlaps."""
        return self.corpus.named_overlaps

    @property
    def description(self) -> str:
        """Return a short golden-corpus description."""
        return f"golden corpus ({self.spec.parameter('dataset', 'default')})"

    def close(self) -> None:
        """Release the wrapped corpus and packaged-resource context once."""
        if self._closed:
            return
        self._closed = True
        try:
            self.corpus.close()
        finally:
            if self.resource_context is not None:
                self.resource_context.__exit__(None, None, None)
