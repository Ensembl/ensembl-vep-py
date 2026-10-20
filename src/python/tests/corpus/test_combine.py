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
"""Tests for corpus composition."""

from __future__ import annotations

import pytest

from ensembl.vep.corpus.combine import combine_corpora
from ensembl.vep.corpus.registry import build_corpus
from ensembl.vep.corpus.spec import CorpusSpec


class TestCombinedCorpus:
    """Several corpora presented through one contract."""

    def test_rejects_an_empty_component_collection(self) -> None:
        """A combined corpus must have actual data to expose."""
        with pytest.raises(ValueError, match="at least one"):
            combine_corpora(CorpusSpec(kind="synthetic.all"), ())

    def test_rejects_non_tuple_components(self) -> None:
        """Require callers to materialise component membership."""
        default = build_corpus(CorpusSpec(kind="synthetic.default", fuzz=0))
        with pytest.raises(TypeError, match="tuple"):
            combine_corpora(
                CorpusSpec(kind="synthetic.all", fuzz=0),
                [default],  # type: ignore[arg-type]
            )

    def test_flattens_components_in_order(self) -> None:
        """The combined view preserves stable member ordering."""
        default = build_corpus(CorpusSpec(kind="synthetic.default", fuzz=0))
        heavy = build_corpus(CorpusSpec(kind="synthetic.heavy", step=600))
        combined = combine_corpora(CorpusSpec(kind="synthetic.all", fuzz=0, step=600), (default, heavy))
        assert combined.sequence_regions == (*default.sequence_regions, *heavy.sequence_regions)
        assert combined.named_variants == (*default.named_variants, *heavy.named_variants)
