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
"""Tests for materialised corpora and their owned resources."""

from __future__ import annotations

import pytest

from ensembl.vep.corpus.loaded import LoadedCorpus
from ensembl.vep.corpus.loaded.base import collect_sequence_regions
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.corpus.synthetic.default import default_corpus


class Resource:
    """Record how often a loaded resource is closed."""

    def __init__(self, close_error: Exception | None = None) -> None:
        self.close_count = 0
        self.close_error = close_error

    def close(self) -> None:
        """Record one release."""
        self.close_count += 1
        if self.close_error is not None:
            raise self.close_error


def _materialise(resource: Resource) -> LoadedCorpus:
    """Copy a synthetic corpus into the loaded-corpus contract."""
    source = default_corpus(fuzz_count=0)
    return LoadedCorpus(
        spec=CorpusSpec(kind="test.loaded"),
        sequence_regions=source.sequence_regions,
        named_variants=source.named_variants,
        transcripts=source.transcripts,
        named_overlaps=source.named_overlaps,
        resources=(resource,),
    )


class TestLoadedCorpus:
    """Materialised membership and resource ownership."""

    def test_preserves_materialised_membership(self) -> None:
        """Consumers see the same stable collections through the common contract."""
        corpus = _materialise(Resource())
        assert corpus.spec.kind == "test.loaded"
        assert len(corpus.sequence_regions) == 1
        assert len(corpus.transcripts) == 4
        assert corpus.description == "loaded corpus (test.loaded)"

    def test_closes_owned_resources_once(self) -> None:
        """Repeated cleanup is safe for callers and context managers."""
        resource = Resource()
        corpus = _materialise(resource)
        corpus.close()
        corpus.close()
        assert resource.close_count == 1

    def test_context_manager_releases_resources(self) -> None:
        """Loaded backends are released when a consumer leaves its scope."""
        resource = Resource()
        with _materialise(resource):
            assert resource.close_count == 0
        assert resource.close_count == 1

    def test_failed_construction_closes_resources(self) -> None:
        """Validation failure releases resources without hiding its original error."""
        resource = Resource(close_error=RuntimeError("resource cleanup failed"))
        with pytest.raises(ValueError, match="a corpus must contain at least one sequence region"):
            LoadedCorpus(
                spec=CorpusSpec(kind="test.loaded"),
                sequence_regions=(),
                named_variants=(),
                transcripts=(),
                named_overlaps=(),
                resources=(resource,),
            )
        assert resource.close_count == 1

    def test_collects_sequence_regions_once_in_encounter_order(self) -> None:
        """A loader can derive sequence-region membership from loaded features."""
        source = default_corpus(fuzz_count=0)
        features = (source.transcripts[0], source.named_variants[0][1], source.transcripts[1])
        assert collect_sequence_regions(features) == source.sequence_regions
