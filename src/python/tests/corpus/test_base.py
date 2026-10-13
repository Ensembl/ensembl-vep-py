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
"""Tests for the corpus contract."""

from __future__ import annotations

from ensembl.vep.corpus.spec import CorpusSpec
from .test_registry import EmptyCorpus


class TestCorpusContext:
    """The common resource-management surface."""

    def test_context_manager_returns_the_same_corpus(self) -> None:
        """Consumers can uniformly scope loaded and generated corpora."""
        corpus = EmptyCorpus(CorpusSpec(kind="test.empty"))
        with corpus as managed:
            assert managed is corpus

    def test_description_identifies_the_registered_kind(self) -> None:
        """Every implementation inherits a useful default description."""
        corpus = EmptyCorpus(CorpusSpec(kind="test.empty"))
        assert corpus.description == "corpus (test.empty)"
