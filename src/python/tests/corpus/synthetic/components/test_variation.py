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
"""Tests for synthetic variant generation."""

from __future__ import annotations

import hashlib
import json

import pytest

from ensembl.vep.corpus.synthetic.components.variation import boundary_variant_alleles, fuzz_variant_alleles
from ensembl.vep.corpus.synthetic.default import build_genome


class TestVariantComponents:
    """Boundary and fuzz variant streams."""

    def test_boundary_stream_has_a_stable_seeded_signature(self) -> None:
        """Names, coordinates, alleles and order remain a regression fixture."""
        sequence_region, transcripts = build_genome()
        signature = [
            (name, allele.start, allele.end, allele.reference_allele, allele.alternate_allele)
            for name, allele in boundary_variant_alleles(sequence_region, transcripts)
        ]
        encoded = json.dumps(signature, separators=(",", ":")).encode()
        assert len(signature) == 655
        assert (
            hashlib.sha256(encoded).hexdigest()
            == "7b0662d7adee7623ce148a37ea08c7af6f300693c0b696cbc6e95960573c0ed8"
        )

    def test_fuzz_count_must_not_be_negative(self) -> None:
        """The low-level generator rejects impossible request sizes itself."""
        sequence_region, _transcripts = build_genome()
        with pytest.raises(ValueError, match="variant_count"):
            fuzz_variant_alleles(sequence_region, variant_count=-1)
