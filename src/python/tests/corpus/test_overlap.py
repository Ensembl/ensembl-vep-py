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
"""Tests for generic overlap construction."""

from __future__ import annotations

import pytest

from ensembl.vep.corpus.overlap import build_transcript_overlaps
from ensembl.vep.corpus.synthetic.default import build_genome
from ensembl.vep.corpus.synthetic.components.variation import fuzz_variant_alleles


class TestTranscriptOverlaps:
    """The generic transcript overlap builder."""

    def test_rejects_a_negative_window(self) -> None:
        """A matching distance cannot be negative."""
        sequence_region, transcripts = build_genome()
        variants = fuzz_variant_alleles(sequence_region, variant_count=1)
        with pytest.raises(ValueError, match="window"):
            build_transcript_overlaps(variants, transcripts, window=-1)
