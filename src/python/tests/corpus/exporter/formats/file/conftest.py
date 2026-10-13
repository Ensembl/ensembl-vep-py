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
"""Shared fixtures for file-format corpus exporter tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from ensembl.vep.corpus import Corpus, build_corpus
from ensembl.vep.corpus.exporter import get_exporter
from ensembl.vep.corpus.exporter.manifest import CorpusManifest
from ensembl.vep.corpus.spec import CorpusSpec


@pytest.fixture(name="source", scope="module")
def fixture_source() -> Corpus:
    """Build both small synthetic corpora once for file-format assertions."""
    return build_corpus(CorpusSpec(kind="synthetic.all", fuzz=5, step=600))


@pytest.fixture(name="saved", scope="module")
def fixture_saved(source: Corpus, tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, CorpusManifest]:
    """Export the shared corpus once and return its directory and manifest."""
    directory = tmp_path_factory.mktemp("corpus")
    return directory, get_exporter("file").export(source, directory)
