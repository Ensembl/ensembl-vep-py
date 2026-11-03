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

"""Shared fixtures for cache tests."""

# pylint: disable=redefined-outer-name

from __future__ import annotations

import pytest

from ensembl.vep.corpus.synthetic import SyntheticCorpus
from ensembl.vep.corpus.synthetic.default import build_genome, default_corpus
from ensembl.vep.corpus.synthetic.heavy import heavy_corpus
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.sequence import Slice

_TEST_FUZZ_COUNT = 300


@pytest.fixture(scope="session")
def synthetic_genome() -> tuple[Slice, tuple[Transcript, ...]]:
    """Create the synthetic chromosome slice and transcripts."""
    return build_genome()


@pytest.fixture(scope="session")
def synthetic_transcripts_by_stable_id(
    synthetic_genome: tuple[Slice, tuple[Transcript, ...]],
) -> dict[str, Transcript]:
    """Create synthetic transcripts keyed by stable ID."""
    _, transcripts = synthetic_genome
    return {transcript.stable_id: transcript for transcript in transcripts}


@pytest.fixture(scope="session")
def small_default_corpus() -> SyntheticCorpus:
    """Create a small synthetic annotation corpus."""
    return default_corpus(fuzz_count=_TEST_FUZZ_COUNT)


@pytest.fixture(scope="session")
def full_heavy_locus_corpus() -> SyntheticCorpus:
    """Create a synthetic corpus with shared coding regions."""
    return heavy_corpus()
