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
"""The synthetic corpus model."""

from __future__ import annotations

from attrs import field, frozen
from attrs.validators import instance_of

from ensembl.vep.corpus.base import Corpus, NamedOverlap, NamedVariant
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.corpus.validation import validate_corpus
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.sequence import Slice
from ensembl.vep.models.variation import VariationAllele

__all__ = [
    "SyntheticCorpus",
]


@frozen(kw_only=True, auto_attribs=True)
class SyntheticCorpus(Corpus):
    """A synthetic chromosome slice, its transcripts, variants, and named variant overlaps.

    Attributes:
        spec (CorpusSpec): Specification used to generate the corpus.
        synthetic_slice (Slice): The whole-chromosome slice every feature and variant is
            expressed against.
        named_variants (tuple[tuple[str, VariationAllele], ...]): Every generated variant, including
            variants that do not overlap any transcript.
        transcripts (tuple[Transcript, ...]): Every transcript in the corpus.
        named_overlaps (tuple[tuple[str, TranscriptVariationAllele], ...]): (name, overlap) pairs,
            suitable for ``pytest.mark.parametrize`` IDs.
    """

    _synthetic_slice: Slice = field(alias="synthetic_slice")
    _named_variants: tuple[tuple[str, VariationAllele], ...] = field(
        validator=instance_of(tuple), alias="named_variants"
    )
    _transcripts: tuple[Transcript, ...] = field(
        factory=tuple, validator=instance_of(tuple), alias="transcripts"
    )
    _named_overlaps: tuple[tuple[str, TranscriptVariationAllele], ...] = field(
        factory=tuple, validator=instance_of(tuple), alias="named_overlaps"
    )
    _spec: CorpusSpec = field(factory=lambda: CorpusSpec(kind="synthetic.default"), alias="spec")

    def __attrs_post_init__(self) -> None:
        """Validate identities shared by the variant and overlap views of the corpus."""
        validate_corpus(
            self.sequence_regions,
            self.named_variants,
            self.transcripts,
            self.named_overlaps,
        )

    @property
    def spec(self) -> CorpusSpec:
        """Return the specification that generated the corpus."""
        return self._spec

    @property
    def synthetic_slice(self) -> Slice:
        """Return the synthetic sequence region."""
        return self._synthetic_slice

    @property
    def sequence_regions(self) -> tuple[Slice, ...]:
        """Return the synthetic sequence region as a one-item sequence."""
        return (self._synthetic_slice,)

    @property
    def named_variants(self) -> tuple[NamedVariant, ...]:
        """Return every generated variant."""
        return self._named_variants

    @property
    def transcripts(self) -> tuple[Transcript, ...]:
        """Return every generated transcript."""
        return self._transcripts

    @property
    def named_overlaps(self) -> tuple[NamedOverlap, ...]:
        """Return every generated transcript overlap."""
        return self._named_overlaps

    @property
    def description(self) -> str:
        """Return a short synthetic-corpus description."""
        return f"synthetic corpus ({self.spec.kind})"
