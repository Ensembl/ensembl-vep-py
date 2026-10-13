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
"""Corpus construction and management.

A corpus contains variants and their ``(variant, feature)`` overlaps, however they were obtained.
Keeping variants explicit preserves inputs with no feature overlap instead of forcing consumers to
reconstruct variants from derived overlap records. The synthetic generator lives here rather than
under ``annotators/`` because it depends only on ``common`` and ``models`` and backs the
session-scoped fixtures of the whole test suite, not one subsystem.
"""

from ensembl.vep.corpus.base import Corpus, NamedOverlap, NamedVariant
from ensembl.vep.corpus.golden import GoldenCorpus, load_golden_corpus
from ensembl.vep.corpus.loaded import CorpusResource, LoadedCorpus
from ensembl.vep.corpus.registry import (
    CorpusOptionsError,
    available_corpora,
    build_corpus,
    build_corpus_from_options,
    register_corpus,
)
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.corpus.synthetic import SyntheticCorpus

__all__ = [
    "Corpus",
    "CorpusOptionsError",
    "CorpusResource",
    "CorpusSpec",
    "GoldenCorpus",
    "LoadedCorpus",
    "NamedOverlap",
    "NamedVariant",
    "SyntheticCorpus",
    "available_corpora",
    "build_corpus",
    "build_corpus_from_options",
    "load_golden_corpus",
    "register_corpus",
]
