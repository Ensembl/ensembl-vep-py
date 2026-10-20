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
"""Tests for the default synthetic corpus recipe."""

from __future__ import annotations

from ensembl.vep.corpus.synthetic.default import default_corpus


class TestDefaultCorpus:
    """Default recipe stability and sizing."""

    def test_overlap_membership_can_be_read_repeatedly(self) -> None:
        """Materialised overlap order is stable across repeated reads."""
        corpus = default_corpus(fuzz_count=0)
        first = [name for name, _overlap in corpus.named_overlaps]
        second = [name for name, _overlap in corpus.named_overlaps]
        assert first == second

    def test_fuzz_count_is_honoured(self) -> None:
        """The recipe adds exactly the requested number of fuzz variants."""
        small = default_corpus(fuzz_count=0)
        larger = default_corpus(fuzz_count=50)
        assert len(larger.named_variants) == len(small.named_variants) + 50
        assert len(larger.named_overlaps) > len(small.named_overlaps)
