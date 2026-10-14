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
"""Tests for corpus factory registration."""

from __future__ import annotations

import pytest

from ensembl.vep.corpus import Corpus, available_corpora, build_corpus, register_corpus
from ensembl.vep.corpus.base import NamedOverlap, NamedVariant
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.sequence import Slice


class EmptyCorpus(Corpus):
    """Small concrete corpus used to exercise the registry contract."""

    def __init__(self, specification: CorpusSpec) -> None:
        self._specification = specification

    @property
    def spec(self) -> CorpusSpec:
        """Return the construction specification."""
        return self._specification

    @property
    def sequence_regions(self) -> tuple[Slice, ...]:
        """Return no sequence regions."""
        return ()

    @property
    def named_variants(self) -> tuple[NamedVariant, ...]:
        """Return no variants."""
        return ()

    @property
    def transcripts(self) -> tuple[Transcript, ...]:
        """Return no transcripts."""
        return ()

    @property
    def named_overlaps(self) -> tuple[NamedOverlap, ...]:
        """Return no overlaps."""
        return ()


class TestCorpusRegistry:
    """Canonical registered names and aliases."""

    def test_registers_and_builds_a_corpus(self) -> None:
        """The registry passes the complete specification to its factory."""
        name = "test.registry.build"

        @register_corpus(name)
        def factory(specification: CorpusSpec) -> Corpus:
            return EmptyCorpus(specification)

        specification = CorpusSpec(kind=name, parameters={"seed": 7})
        corpus = build_corpus(specification)
        assert corpus.spec == specification
        assert name in available_corpora()

    def test_unknown_name_lists_available_corpora(self) -> None:
        """Lookup failures tell callers how to correct the requested name."""
        name = "test.registry.available"
        register_corpus(name)(EmptyCorpus)
        with pytest.raises(KeyError, match=f"Unknown corpus 'missing.corpus'.*{name}"):
            build_corpus(CorpusSpec(kind="missing.corpus"))

    def test_registration_is_atomic_when_an_alias_conflicts(self) -> None:
        """A failed registration does not leave a partially registered factory."""
        name = "test.registry.canonical"
        register_corpus(name, aliases=("test.registry.alias",))(EmptyCorpus)
        with pytest.raises(ValueError, match="already registered"):
            register_corpus("test.registry.rejected", aliases=("test.registry.alias",))(EmptyCorpus)
        assert "test.registry.rejected" not in available_corpora()
