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

from ensembl.vep.corpus import Corpus
from ensembl.vep.corpus.registry import (
    available_corpora,
    build_corpus,
    register_corpus,
    corpus_spec_from_options,
    build_corpus_from_options,
)
from ensembl.vep.corpus.spec import CorpusSpec


class TestCorpusRegistry:
    """Canonical registered names."""

    def test_lists_builtin_corpora(self) -> None:
        """The installed construction strategies are discoverable."""
        assert {"synthetic.all", "synthetic.default", "synthetic.heavy"} <= set(available_corpora())

    def test_unknown_name_lists_available_corpora(self) -> None:
        """Lookup failures tell callers how to correct the requested name."""
        with pytest.raises(KeyError, match="Unknown corpus 'missing.corpus'.*synthetic.default"):
            build_corpus(CorpusSpec(kind="missing.corpus"))

    def test_registration_rejects_empty_names(self) -> None:
        """Names cannot be empty strings."""
        with pytest.raises(ValueError, match="Corpus registry names must not be empty"):
            register_corpus("", aliases=("synthetic.default",))
        with pytest.raises(ValueError, match="Corpus registry names must not be empty"):
            register_corpus("synthetic.default", aliases=("default", ""))

    def test_registration_rejects_repeated_names(self) -> None:
        """Names cannot be repeated."""

        def factory(_specification: CorpusSpec) -> Corpus:
            return build_corpus(CorpusSpec(kind="synthetic.default"))

        with pytest.raises(ValueError, match="Corpus registration repeats a name"):
            register_corpus("synthetic.default", aliases=("synthetic.default",))(factory)

    def test_registration_of_aliases(self) -> None:
        """Register a corpus factory under an alias."""

        def factory(_specification: CorpusSpec) -> Corpus:
            return build_corpus(CorpusSpec(kind="synthetic.default"))

        register_corpus("test", aliases=("test_alias",))(factory)

    def test_registration_is_atomic_when_an_alias_conflicts(self) -> None:
        """A failed registration does not leave a partially registered factory."""

        def factory(_specification: CorpusSpec) -> Corpus:
            raise AssertionError("the rejected factory must never become visible")

        with pytest.raises(ValueError, match="synthetic.default.*already registered"):
            register_corpus("test.rejected", aliases=("synthetic.default",))(factory)
        assert "test.rejected" not in available_corpora()

    @pytest.mark.parametrize("kind", ["default", "heavy", "both"])
    def test_legacy_aliases_are_not_registered(self, kind: str) -> None:
        """Only canonical dotted names resolve through the registry."""
        with pytest.raises(KeyError, match=f"Unknown corpus '{kind}'"):
            build_corpus(CorpusSpec(kind=kind))

    @pytest.mark.parametrize(
        ("kind", "sequence_regions"),
        [("synthetic.default", 1), ("synthetic.heavy", 1), ("synthetic.all", 2)],
    )
    def test_registered_synthetic_names_build(self, kind: str, sequence_regions: int) -> None:
        """Every canonical synthetic name resolves through the same registry path."""
        corpus = build_corpus(CorpusSpec(kind=kind, fuzz=0, step=600))
        assert len(corpus.sequence_regions) == sequence_regions

    def test_corpus_spec_from_options(self) -> None:
        """Test creating CorpusSpec from an options dict"""
        _ = corpus_spec_from_options({"kind": "loaded.files", "corpus_dir": "test"})

    def test_build_corpus_from_options(self) -> None:
        """Test creating Corpus from an options dict"""
        _ = build_corpus_from_options({"kind": "synthetic.default", "corpus_dir": "test"})
