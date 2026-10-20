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
"""Tests for golden corpus resolution."""

from __future__ import annotations

from pathlib import Path

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the annotators and the facts the catalogue compiles against.
from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.corpus import Corpus, build_corpus
from ensembl.vep.corpus.golden import (
    GoldenCorpus,
    GoldenReport,
    load_golden_corpus,
    record_golden,
)
from ensembl.vep.corpus.spec import CorpusSpec

PACKAGED_GOLDEN_DIRECTORY = (
    Path(__file__).parents[3] / "ensembl" / "vep" / "corpus" / "golden" / "data" / "default"
)


class TestGoldenResolution:
    """Loaded golden data and its explicit generation fallback."""

    @staticmethod
    def _spec() -> CorpusSpec:
        """Return the packaged fixture's specification."""
        return CorpusSpec(kind="golden.file", parameters={"dataset": "default"})

    def test_packaged_fixture_loads_as_a_file_backed_corpus(self) -> None:
        """Load the checked-in bundle rather than regenerating it."""
        golden = build_corpus(self._spec())

        assert isinstance(golden, GoldenCorpus)
        assert isinstance(golden.corpus, Corpus)
        assert golden.sequence_regions
        assert golden.named_variants
        assert golden.transcripts
        assert golden.named_overlaps
        golden.close()

    @pytest.mark.parametrize("parameters", [{}, {"allow_generate": True}])
    def test_packaged_defaults_have_one_stable_golden_identity(self, parameters: dict[str, bool]) -> None:
        """Keep pinned identity independent of defaulted fallback policy."""
        with build_corpus(CorpusSpec(kind="golden.file", parameters=parameters)) as golden:
            assert isinstance(golden, GoldenCorpus)
            recorded = record_golden(
                golden,
                build_annotator(annotator=golden.expected_terms.annotator),
                golden.expected_terms.annotator,
            )
            report = GoldenReport.compare(golden.expected_terms, recorded)
        assert report.ok, "\n".join(report.lines())

    def test_missing_directory_generates_only_when_explicitly_allowed(self, tmp_path: Path) -> None:
        """Allow an absent optional bundle to use an explicit fallback."""
        missing_directory = tmp_path / "missing"
        fallback = CorpusSpec(kind="synthetic.default", parameters={"fuzz": 0})

        with pytest.raises(FileNotFoundError, match="missing"):
            load_golden_corpus(self._spec(), directory=missing_directory, fallback=fallback)

        golden = load_golden_corpus(
            self._spec(),
            directory=missing_directory,
            fallback=fallback,
            allow_generate=True,
        )
        assert golden.corpus.spec == fallback
        assert golden.expected_terms.terms
        golden.close()

    def test_malformed_existing_bundle_never_silently_generates(self, tmp_path: Path) -> None:
        """Reject a partial bundle even when fallback generation is enabled."""
        malformed_directory = tmp_path / "malformed"
        malformed_directory.mkdir()
        fallback = CorpusSpec(kind="synthetic.default", parameters={"fuzz": 0})

        with pytest.raises(FileNotFoundError, match="corpus.fa"):
            load_golden_corpus(
                self._spec(),
                directory=malformed_directory,
                fallback=fallback,
                allow_generate=True,
            )

    def test_explicit_bundle_cannot_replace_requested_identity(self) -> None:
        """Reject pinned terms identifying another requested dataset."""
        specification = CorpusSpec(kind="golden.file", parameters={"dataset": "another"})
        with pytest.raises(ValueError, match="not requested"):
            load_golden_corpus(specification, directory=PACKAGED_GOLDEN_DIRECTORY)
