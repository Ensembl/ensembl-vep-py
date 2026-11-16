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
"""Shared fixtures for trace tests."""

# pylint: disable=redefined-outer-name

from __future__ import annotations

import pytest

from ensembl.vep.corpus.synthetic import (
    boundary_variant_alleles,
    build_synthetic_genome,
    build_transcript_variation_alleles,
)
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.sequence import Slice


@pytest.fixture(scope="session")
def synthetic_genome() -> tuple[Slice, list[Transcript]]:
    """Create the synthetic chromosome slice and transcripts."""
    return build_synthetic_genome()


@pytest.fixture(scope="session")
def boundary_overlaps(
    synthetic_genome: tuple[Slice, list[Transcript]],
) -> dict[str, TranscriptVariationAllele]:
    """Create boundary overlaps keyed by transcript and variant name."""
    synthetic_slice, transcripts = synthetic_genome
    variants = boundary_variant_alleles(synthetic_slice, transcripts)
    overlaps = build_transcript_variation_alleles(variants, transcripts)
    return {f"{overlap.feature.stable_id}::{name}": overlap for name, overlap in overlaps}
