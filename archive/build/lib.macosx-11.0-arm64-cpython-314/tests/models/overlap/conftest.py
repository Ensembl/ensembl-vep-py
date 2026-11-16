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
"""Shared fixtures for the overlap model tests."""

from __future__ import annotations

from collections.abc import Callable

import pytest

from ensembl.vep.common import CoordSystemLevel, ExonPhase, Strand
from ensembl.vep.models.genic import Exon, Transcript
from ensembl.vep.models.sequence import CoordSystem, Slice

SINGLE_EXON_CODING_SEQUENCE = "ATG" + "GGG" * 10 + "TAA"
"""The CDS of ``single_exon_transcript``, genomic positions 4 to 39 on the forward strand."""

SINGLE_EXON_CDS_START = 4
"""Genomic position of the first coding base of ``single_exon_transcript``."""

SingleExonTranscript = Callable[[ExonPhase | None], Transcript]
"""Builds the single-exon transcript, its exon carrying the given phase."""


def _slice(sequence: str) -> Slice:
    """Return a slice carrying ``sequence`` from position 1."""
    genomic_slice = Slice(
        coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
        seq_region_name="chrSingle",
        start=1,
        end=len(sequence),
        strand=Strand.FORWARD,
    )
    with genomic_slice.unsafe() as unsafe_slice:
        setattr(unsafe_slice, "_seq", sequence)
    return genomic_slice


def _build_single_exon_transcript(phase: ExonPhase | None) -> Transcript:
    """Return a forward transcript of one exon at 1-40, coding from 4 to 39."""
    genomic_slice = _slice("NNN" + SINGLE_EXON_CODING_SEQUENCE + "N")
    transcript = Transcript(
        start=1, end=40, strand=Strand.FORWARD, stable_id="ENST_SINGLE", slice=genomic_slice
    )
    exon = Exon(start=1, end=40, strand=Strand.FORWARD, stable_id="ENST_SINGLE_e0", slice=genomic_slice)
    exon.phase = phase
    transcript.spliced_exons = [exon]
    transcript.seed_attrs(
        first_coding_exon_index=0,
        last_coding_exon_index=0,
        exon_coding_start=SINGLE_EXON_CDS_START,
        exon_coding_end=39,
    )
    return transcript


@pytest.fixture(name="single_exon_transcript")
def fixture_single_exon_transcript() -> SingleExonTranscript:
    """A builder for the single-exon coding transcript, given its exon's phase."""
    return _build_single_exon_transcript
