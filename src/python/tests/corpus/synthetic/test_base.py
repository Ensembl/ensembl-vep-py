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
"""Tests for the synthetic corpus model."""

from __future__ import annotations

import pytest

from ensembl.vep.corpus.synthetic import SyntheticCorpus
from ensembl.vep.corpus.synthetic.default import default_corpus


class TestSyntheticCorpus:
    """Materialised synthetic corpus membership."""

    def test_named_variants_are_required(self) -> None:
        """Overlap-only construction fails instead of exporting an empty variant collection."""
        corpus = default_corpus(fuzz_count=0)
        with pytest.raises(TypeError, match="named_variants"):
            SyntheticCorpus(  # type: ignore[call-arg]
                synthetic_slice=corpus.synthetic_slice,
                transcripts=corpus.transcripts,
                named_overlaps=corpus.named_overlaps,
            )

    def test_rejects_non_tuple_membership(self) -> None:
        """Require builders to materialise membership before construction."""
        corpus = default_corpus(fuzz_count=0)
        with pytest.raises(TypeError, match="tuple"):
            SyntheticCorpus(
                synthetic_slice=corpus.synthetic_slice,
                named_variants=list(corpus.named_variants),  # type: ignore[arg-type]
                transcripts=corpus.transcripts,
                named_overlaps=corpus.named_overlaps,
            )
