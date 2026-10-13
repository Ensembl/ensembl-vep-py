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
"""Tests for the golden corpus model."""

from __future__ import annotations

import pytest

from ensembl.vep.corpus.golden import GoldenCorpus, load_golden_corpus
from ensembl.vep.corpus.spec import CorpusSpec


class _TrackingResource:
    """Minimal packaged-resource context used to verify ownership."""

    exits: int = 0

    def __enter__(self) -> _TrackingResource:
        """Return the context-managed resource."""
        return self

    def __exit__(self, *_exc_info: object) -> None:
        """Count cleanup calls."""
        self.exits += 1


class TestGoldenCorpus:
    """Wrapped corpus delegation and resource ownership."""

    def test_delegates_properties_and_closes_resources_once(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Delegate corpus views and close both owned resources once."""
        specification = CorpusSpec(kind="golden.file", parameters={"dataset": "default"})
        loaded = load_golden_corpus(specification)
        close_calls = 0

        def close_underlying() -> None:
            """Record delegated cleanup."""
            nonlocal close_calls
            close_calls += 1

        monkeypatch.setattr(loaded.corpus, "close", close_underlying)
        resource = _TrackingResource()
        golden = GoldenCorpus(
            spec=specification,
            corpus=loaded.corpus,
            expected_terms=loaded.expected_terms,
            resource_context=resource,
        )

        assert golden.sequence_regions is loaded.corpus.sequence_regions
        assert golden.named_variants is loaded.corpus.named_variants
        assert golden.transcripts is loaded.corpus.transcripts
        assert golden.named_overlaps is loaded.corpus.named_overlaps
        golden.close()
        golden.close()

        assert close_calls == 1
        assert resource.exits == 1
        loaded.close()
