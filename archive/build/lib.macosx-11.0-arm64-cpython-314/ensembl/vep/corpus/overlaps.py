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
"""Overlaps from a generated corpus."""

from __future__ import annotations

from collections.abc import Iterator
from functools import cached_property

from attrs import define, field

from ensembl.vep.corpus.base import Corpus
from ensembl.vep.corpus.combine import combine_corpora
from ensembl.vep.corpus.spec import KINDS, CorpusSpec
from ensembl.vep.corpus.synthetic import SyntheticCorpus, default_corpus, heavy_locus_corpus
from ensembl.vep.models.overlap import TranscriptVariationAllele

__all__ = ["KINDS", "SyntheticOverlapSource"]


@define(init=False, frozen=True)
class SyntheticOverlapSource:
    """Overlaps generated from one or both synthetic corpora.

    The corpora are built on first use and kept, since a caller that walks the overlaps twice --
    a conformance gate then a timing pass, say -- must see the same ones both times. ``CorpusSpec``
    is the single source of configuration. The individual keyword arguments remain as a
    compatibility convenience and override the corresponding value in ``spec`` when supplied.

    Attributes:
        kind (str): Which corpus to draw from, one of ``KINDS``.
        fuzz (int): Random variant count for the boundary corpus.
        step (int): Position increment for the heavy-locus sweep.
        genome_seed (int): Seed for the default corpus's chromosome.
        boundary_seed (int): Seed for the boundary variants.
        fuzz_seed (int): Seed for the fuzz variants.
        heavy_genome_seed (int): Seed for the heavy locus's chromosome.
        sweep_seed (int): Seed for the heavy-locus sweep.
    """

    _spec: CorpusSpec = field(alias="spec")

    def __init__(
        self,
        *,
        spec: CorpusSpec | None = None,
        kind: str | None = None,
        fuzz: int | None = None,
        step: int | None = None,
        genome_seed: int | None = None,
        boundary_seed: int | None = None,
        fuzz_seed: int | None = None,
        heavy_genome_seed: int | None = None,
        sweep_seed: int | None = None,
    ) -> None:
        """Create a source from a spec, optionally overriding individual parameters.

        The flat keywords preserve the original public constructor while keeping defaults and
        validation in ``CorpusSpec`` rather than duplicating them here.
        """
        base = spec or CorpusSpec()
        resolved_kind = kind if kind is not None else base.kind
        if resolved_kind not in (*KINDS, "synthetic.default", "synthetic.heavy", "synthetic.all"):
            raise ValueError(f"unknown synthetic corpus kind {resolved_kind!r}")
        resolved_spec = CorpusSpec(
            kind=resolved_kind,
            fuzz=fuzz if fuzz is not None else base.fuzz,
            step=step if step is not None else base.step,
            genome_seed=genome_seed if genome_seed is not None else base.genome_seed,
            boundary_seed=boundary_seed if boundary_seed is not None else base.boundary_seed,
            fuzz_seed=fuzz_seed if fuzz_seed is not None else base.fuzz_seed,
            heavy_genome_seed=(
                heavy_genome_seed if heavy_genome_seed is not None else base.heavy_genome_seed
            ),
            sweep_seed=sweep_seed if sweep_seed is not None else base.sweep_seed,
        )
        object.__setattr__(self, "_spec", resolved_spec)

    @classmethod
    def from_spec(cls, spec: CorpusSpec) -> SyntheticOverlapSource:
        """Return the source a spec describes."""
        return cls(spec=spec)

    @property
    def spec(self) -> CorpusSpec:
        """The spec this source builds from."""
        return self._spec

    @property
    def kind(self) -> str:
        """Which corpus kinds to build."""
        return self.spec.kind

    @property
    def fuzz(self) -> int:
        """Random variant count for the default corpus."""
        return self.spec.fuzz

    @property
    def step(self) -> int:
        """Position increment for the heavy-locus sweep."""
        return self.spec.step

    @property
    def genome_seed(self) -> int:
        """Seed for the default corpus chromosome."""
        return self.spec.genome_seed

    @property
    def boundary_seed(self) -> int:
        """Seed for the boundary variants."""
        return self.spec.boundary_seed

    @property
    def fuzz_seed(self) -> int:
        """Seed for the fuzz variants."""
        return self.spec.fuzz_seed

    @property
    def heavy_genome_seed(self) -> int:
        """Seed for the heavy-locus chromosome."""
        return self.spec.heavy_genome_seed

    @property
    def sweep_seed(self) -> int:
        """Seed for the heavy-locus variant sweep."""
        return self.spec.sweep_seed

    @property
    def description(self) -> str:
        """Describe the corpora this draws from."""
        return f"synthetic corpus ({self.kind})"

    @property
    def corpora(self) -> tuple[SyntheticCorpus, ...]:
        """Public accessor for the built corpora."""
        return self._corpora

    @cached_property
    def corpus(self) -> Corpus:
        """Return the compatibility source as one consumer-facing corpus."""
        if len(self.corpora) == 1:
            return self.corpora[0]
        return combine_corpora(self.spec, self.corpora)

    @cached_property
    def _corpora(self) -> tuple[SyntheticCorpus, ...]:
        """Build the corpora this source draws from.

        Raises:
            ValueError: If ``kind`` names no corpus.
        """
        kinds = self.spec.kinds
        built = []
        if "default" in kinds:
            built.append(
                default_corpus(
                    fuzz_count=self.fuzz,
                    seed=self.genome_seed,
                    boundary_seed=self.boundary_seed,
                    fuzz_seed=self.fuzz_seed,
                    spec=self.spec,
                )
            )
        if "heavy" in kinds:
            built.append(
                heavy_locus_corpus(
                    step=self.step,
                    seed=self.heavy_genome_seed,
                    sweep_seed=self.sweep_seed,
                    spec=self.spec,
                )
            )
        return tuple(built)

    def named_overlaps(self) -> Iterator[tuple[str, TranscriptVariationAllele]]:
        """Yield every overlap the corpora hold, named as the generator named it."""
        for corpus in self.corpora:
            yield from corpus.named_overlaps
