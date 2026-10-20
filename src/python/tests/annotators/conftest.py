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

# pylint: disable=redefined-outer-name
"""Shared fixtures for the annotators test suite."""

from __future__ import annotations

import pytest

from ensembl.vep.corpus.overlap import build_transcript_overlaps
from ensembl.vep.corpus.synthetic import SyntheticCorpus
from ensembl.vep.corpus.synthetic.components.variation import boundary_variant_alleles
from ensembl.vep.corpus.synthetic.default import build_genome, default_corpus
from ensembl.vep.corpus.synthetic.heavy import heavy_corpus
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.sequence import Slice

# A small, deterministic fuzz count keeps the test suite fast while still exercising a broad mix
# of random variants alongside the adversarial boundary set.
_TEST_FUZZ_COUNT = 300

# The "whole" corpus, matching the fuzz variant count socore's own reference invocation uses
# (``python3 -m socore.conformance 4000``) -- used by the cross-annotator conformance gate, which
# needs to be thorough rather than merely fast.
_FULL_FUZZ_COUNT = 4000


@pytest.fixture(scope="session")
def synthetic_genome() -> tuple[Slice, tuple[Transcript, ...]]:
    """The synthetic chromosome slice and its four transcripts, built once per test session."""
    return build_genome()


@pytest.fixture(scope="session")
def synthetic_transcripts_by_stable_id(
    synthetic_genome: tuple[Slice, tuple[Transcript, ...]],
) -> dict[str, Transcript]:
    """The four synthetic transcripts, keyed by stable_id, for direct lookup in tests."""
    _slice, transcripts = synthetic_genome
    return {transcript.stable_id: transcript for transcript in transcripts}


@pytest.fixture(scope="session")
def boundary_overlaps(
    synthetic_genome: tuple[Slice, tuple[Transcript, ...]],
) -> dict[str, TranscriptVariationAllele]:
    """Every adversarial boundary-variant overlap, keyed by its descriptive name.

    Built once per session (session-scoped) since the corpus itself is a pure function of a
    fixed seed and is comparatively expensive to build; individual tests only read from it.
    """
    synthetic_slice, transcripts = synthetic_genome
    named_variants = boundary_variant_alleles(synthetic_slice, transcripts)
    named_overlaps = build_transcript_overlaps(named_variants, transcripts)
    overlaps_by_name: dict[str, TranscriptVariationAllele] = {}
    for name, overlap in named_overlaps:
        # A boundary variant can legitimately pair with more than one transcript (e.g. a
        # downstream flank of one transcript sitting inside another); disambiguate by the
        # overlapped transcript's stable_id so every entry stays addressable by name.
        overlaps_by_name[f"{overlap.feature.stable_id}::{name}"] = overlap
    return overlaps_by_name


@pytest.fixture(scope="session")
def small_default_corpus() -> SyntheticCorpus:
    """The full synthetic corpus (boundary + a small fuzz set), built once per session."""
    return default_corpus(fuzz_count=_TEST_FUZZ_COUNT)


@pytest.fixture(scope="session")
def full_default_corpus() -> SyntheticCorpus:
    """The whole boundary+fuzz corpus (4000 fuzz variants), built once per session.

    Used by the cross-annotator conformance gate (``test_conformance.py``), which is meant to be
    thorough over the adversarial four-transcript corpus rather than merely fast.
    """
    return default_corpus(fuzz_count=_FULL_FUZZ_COUNT)


@pytest.fixture(scope="session")
def full_heavy_locus_corpus() -> SyntheticCorpus:
    """The heavy-locus corpus (2 strands x 18 exon-sharing transcripts), built once per session.

    The only corpus that exercises multi-variant-per-transcript zone reuse (and, in a later
    phase, cross-transcript codon sharing) -- see ``ensembl.vep.corpus.synthetic.
    heavy_corpus``.
    """
    return heavy_corpus()
