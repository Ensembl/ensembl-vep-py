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
from shutil import copytree

import pytest

from ensembl.vep.corpus import Corpus, build_corpus
from ensembl.vep.corpus.golden import GoldenCorpus, GoldenTerms, load_golden_corpus
from ensembl.vep.corpus.spec import CorpusSpec

PACKAGED_GOLDEN_DIRECTORY = (
    Path(__file__).parents[3] / "ensembl" / "vep" / "corpus" / "golden" / "data" / "default"
)


class TestGoldenResolution:
    """Loaded golden data and its explicit fallback boundary."""

    @staticmethod
    def _spec() -> CorpusSpec:
        """Return the packaged fixture's specification."""
        return CorpusSpec(kind="golden.file", parameters={"dataset": "default"})

    @staticmethod
    def _copy_bundle(tmp_path: Path) -> Path:
        """Copy the packaged bundle so a test can safely alter its terms."""
        directory = tmp_path / "golden"
        copytree(PACKAGED_GOLDEN_DIRECTORY, directory)
        return directory

    @staticmethod
    def _write_terms(directory: Path, terms: dict[str, tuple[str, ...]]) -> None:
        """Replace a bundle's terms while retaining its pinned provenance."""
        terms_path = directory / "terms.json"
        original = GoldenTerms.load(terms_path)
        GoldenTerms(
            spec=original.spec,
            annotator=original.annotator,
            catalogue=original.catalogue,
            terms=terms,
        ).write(terms_path)

    def test_packaged_fixture_loads_as_a_file_backed_corpus(self) -> None:
        """Load the checked-in bundle rather than regenerating it."""
        golden = build_corpus(self._spec())

        assert isinstance(golden, GoldenCorpus)
        assert isinstance(golden.corpus, Corpus)
        assert golden.sequence_regions
        assert golden.named_variants
        assert golden.transcripts
        assert golden.named_overlaps
        assert golden.expected_terms.terms
        golden.close()

    def test_default_parameters_have_one_stable_golden_identity(self) -> None:
        """Keep pinned identity independent of defaulted loading policy."""
        with build_corpus(CorpusSpec(kind="golden.file")) as golden:
            assert isinstance(golden, GoldenCorpus)
            assert golden.spec == self._spec()

    def test_rejects_non_string_directory_parameter(self) -> None:
        """Do not treat an invalid directory parameter as the packaged default."""
        specification = CorpusSpec(kind="golden.file", parameters={"directory": 123})

        with pytest.raises(TypeError, match="directory must be a string"):
            build_corpus(specification)

    def test_missing_directory_is_not_silently_regenerated(self, tmp_path: Path) -> None:
        """Keep missing-data policy outside the reusable corpus loader."""

        with pytest.raises(FileNotFoundError, match="missing"):
            load_golden_corpus(self._spec(), directory=tmp_path / "missing")

    def test_malformed_existing_bundle_never_generates(self, tmp_path: Path) -> None:
        """Reject a partial bundle rather than masking corruption."""

        with pytest.raises(FileNotFoundError, match="corpus.fa"):
            load_golden_corpus(
                self._spec(),
                directory=tmp_path,
            )

    def test_explicit_bundle_cannot_replace_requested_identity(self) -> None:
        """Reject pinned terms identifying another requested dataset."""
        specification = CorpusSpec(kind="golden.file", parameters={"dataset": "another"})
        with pytest.raises(ValueError, match="not requested"):
            load_golden_corpus(specification, directory=PACKAGED_GOLDEN_DIRECTORY)

    def test_public_loader_ignores_location_when_comparing_identity(self) -> None:
        """Treat the directory as loading policy rather than pinned dataset identity."""
        specification = CorpusSpec(
            kind="golden.file",
            parameters={"dataset": "default", "directory": str(PACKAGED_GOLDEN_DIRECTORY)},
        )

        with load_golden_corpus(specification, directory=PACKAGED_GOLDEN_DIRECTORY) as golden:
            assert golden.spec == self._spec()

    def test_rejects_bundle_missing_terms_for_loaded_overlaps(self, tmp_path: Path) -> None:
        """Every loaded overlap needs an expected-term entry."""
        directory = self._copy_bundle(tmp_path)
        original = GoldenTerms.load(directory / "terms.json")
        missing_key = next(iter(original.terms))
        self._write_terms(
            directory,
            {key: terms for key, terms in original.terms.items() if key != missing_key},
        )

        with pytest.raises(ValueError, match="missing terms for"):
            load_golden_corpus(self._spec(), directory=directory)

    def test_rejects_bundle_with_terms_for_unknown_overlaps(self, tmp_path: Path) -> None:
        """Expected terms must not describe overlaps outside the loaded corpus."""
        directory = self._copy_bundle(tmp_path)
        original = GoldenTerms.load(directory / "terms.json")
        self._write_terms(
            directory,
            {**original.terms, "unknown_transcript::unknown_variant": ()},
        )

        with pytest.raises(ValueError, match="terms for unknown overlaps"):
            load_golden_corpus(self._spec(), directory=directory)
