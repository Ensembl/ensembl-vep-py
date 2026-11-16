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
"""Tests for golden terms."""

from __future__ import annotations

from pathlib import Path

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the annotators and the facts the catalogue compiles against.
from ensembl.vep.annotators.consequence.rules.catalogue import catalogue_digest, catalogue_path
from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.corpus import Corpus, build_corpus
from ensembl.vep.corpus.golden import (
    FORMAT_VERSION,
    GoldenCorpus,
    GoldenReport,
    GoldenTerms,
    load_golden_corpus,
    overlap_key,
    record_golden,
)
from ensembl.vep.corpus.overlaps import SyntheticOverlapSource
from ensembl.vep.corpus.spec import CorpusSpec

GOLDEN_FILE = Path(__file__).parent / "data" / "golden_linear.json"
PACKAGED_GOLDEN_DIRECTORY = (
    Path(__file__).parents[2] / "ensembl" / "vep" / "corpus" / "golden" / "data" / "default"
)
"""The checked-in golden terms the regression gate reads. Re-record with
``vep corpus golden --golden_file <this> --update_golden`` plus the spec in its header."""


@pytest.fixture(name="recorded", scope="module")
def fixture_recorded() -> GoldenTerms:
    """Linear's terms over a small corpus, recorded once."""
    source = SyntheticOverlapSource(kind="both", fuzz=5, step=600)
    return record_golden(source.corpus, build_annotator(annotator="consequence.linear"), "consequence.linear")


class TestRecording:
    """What a recording holds."""

    def test_every_overlap_is_filed_under_its_feature_and_name(self, recorded: GoldenTerms) -> None:
        """One entry per overlap, keyed feature::name."""
        source = SyntheticOverlapSource(kind="both", fuzz=5, step=600)
        expected_keys = {overlap_key(name, overlap) for name, overlap in source.named_overlaps()}
        assert set(recorded.terms) == expected_keys
        assert len(recorded.terms) == len(list(source.named_overlaps()))

    def test_it_names_its_corpus_annotator_and_catalogue(self, recorded: GoldenTerms) -> None:
        """Everything a later check needs to explain a mismatch."""
        assert recorded.spec == CorpusSpec(kind="both", fuzz=5, step=600)
        assert recorded.annotator == "consequence.linear"
        assert recorded.catalogue == catalogue_digest()

    def test_the_terms_are_the_annotators_own(self, recorded: GoldenTerms) -> None:
        """A key's terms are what the annotator emits for that overlap, in order."""
        assert recorded.terms["ENST0001::ENST0001_cds3_del3"] == ("inframe_deletion",)


class TestFile:
    """The file round trip."""

    def test_it_round_trips(self, recorded: GoldenTerms, tmp_path: Path) -> None:
        """Written then loaded is equal."""
        path = tmp_path / "golden.json"
        recorded.write(path)
        assert GoldenTerms.load(path) == recorded

    def test_one_overlap_per_line(self, recorded: GoldenTerms) -> None:
        """So a regeneration diffs as exactly the overlaps that moved."""
        text = recorded.to_json()
        term_lines = [line for line in text.splitlines() if "::" in line]
        assert len(term_lines) == len(recorded.terms)

    def test_another_format_version_is_rejected(self, recorded: GoldenTerms) -> None:
        """An old file is refused, not misread."""
        text = recorded.to_json().replace(f'"format_version": {FORMAT_VERSION}', '"format_version": 999')
        with pytest.raises(ValueError, match="999"):
            GoldenTerms.from_json(text)


class TestReport:
    """Telling the three kinds of mismatch apart."""

    def test_an_identical_recording_is_ok(self, recorded: GoldenTerms) -> None:
        """Nothing to report."""
        report = GoldenReport.compare(recorded, recorded)
        assert report.ok
        assert report.lines() == ["recorded terms match the golden file"]

    def test_a_moved_term_is_a_mismatch(self, recorded: GoldenTerms) -> None:
        """The regression case names the overlap and both term lists."""
        key = "ENST0001::ENST0001_cds3_del3"
        moved = GoldenTerms(
            spec=recorded.spec,
            annotator=recorded.annotator,
            catalogue=recorded.catalogue,
            terms={**recorded.terms, key: ("frameshift_variant",)},
        )
        report = GoldenReport.compare(recorded, moved)
        assert not report.ok
        assert report.mismatches == ((key, ("inframe_deletion",), ("frameshift_variant",)),)
        assert any(key in line for line in report.lines())

    def test_a_changed_catalogue_is_reported_even_when_terms_agree(self, recorded: GoldenTerms) -> None:
        """A stale label is still a problem: the file no longer says what produced it."""
        relabelled = GoldenTerms(
            spec=recorded.spec, annotator=recorded.annotator, catalogue="0" * 16, terms=recorded.terms
        )
        report = GoldenReport.compare(relabelled, recorded)
        assert report.catalogue_changed
        assert not report.ok
        assert "--update_golden" in report.lines()[0]

    def test_a_different_annotator_is_informational(self, recorded: GoldenTerms) -> None:
        """Checking staged against linear's golden terms is a legitimate check."""
        other = GoldenTerms(
            spec=recorded.spec,
            annotator="consequence.staged",
            catalogue=recorded.catalogue,
            terms=recorded.terms,
        )
        report = GoldenReport.compare(recorded, other)
        assert report.annotator_changed
        assert report.ok

    def test_missing_and_unexpected_overlaps_are_listed(self, recorded: GoldenTerms) -> None:
        """A corpus that grew or shrank shows which overlaps."""
        keys = sorted(recorded.terms)
        trimmed = GoldenTerms(
            spec=recorded.spec,
            annotator=recorded.annotator,
            catalogue=recorded.catalogue,
            terms={key: recorded.terms[key] for key in keys[1:]} | {"ENSTX::new": ()},
        )
        report = GoldenReport.compare(recorded, trimmed)
        assert report.missing == (keys[0],)
        assert report.unexpected == ("ENSTX::new",)


class TestCatalogueDigest:
    """The catalogue half of a golden file's provenance."""

    def test_the_digest_follows_the_file_bytes(self, tmp_path: Path) -> None:
        """Same bytes, same digest; one byte more, a different one."""
        copy = tmp_path / "catalogue.json"
        copy.write_bytes(catalogue_path().read_bytes())
        assert catalogue_digest(copy) == catalogue_digest()
        copy.write_bytes(catalogue_path().read_bytes() + b" ")
        assert catalogue_digest(copy) != catalogue_digest()


class TestCheckedInGolden:
    """The regression gate: linear still emits what the checked-in file says."""

    def test_linear_matches_the_checked_in_golden_terms(self) -> None:
        """A failure here is either a regression or a catalogue change; the report says which."""
        golden = GoldenTerms.load(GOLDEN_FILE)
        source = SyntheticOverlapSource.from_spec(golden.spec)
        recorded = record_golden(source.corpus, build_annotator(annotator=golden.annotator), golden.annotator)
        report = GoldenReport.compare(golden, recorded)
        assert report.ok, "\n".join(report.lines())


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
    """The loaded golden corpus and its explicit generation fallback."""

    @staticmethod
    def _spec() -> CorpusSpec:
        """Return the packaged fixture's specification."""
        return CorpusSpec(kind="golden.file", parameters={"dataset": "default"})

    def test_packaged_fixture_loads_as_a_file_backed_corpus(self) -> None:
        """The checked-in bundle is loaded through ``LoadedCorpus``, not regenerated."""
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
        """Defaulted datasets and operational fallback policy do not change pinned identity."""
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
        """An absent optional bundle may use the requested synthetic fallback."""
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
        """A partial bundle is an error even when a fallback was explicitly enabled."""
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

    def test_explicit_bundle_cannot_silently_replace_requested_identity(self) -> None:
        """Pinned terms remain authoritative without disguising one dataset as another."""
        specification = CorpusSpec(kind="golden.file", parameters={"dataset": "another"})
        with pytest.raises(ValueError, match="not requested"):
            load_golden_corpus(specification, directory=PACKAGED_GOLDEN_DIRECTORY)

    def test_delegates_properties_and_closes_resources_once(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Golden wrappers preserve the corpus interface and own both cleanup responsibilities."""
        loaded = load_golden_corpus(self._spec())
        close_calls = 0

        def close_underlying() -> None:
            """Record delegated cleanup."""
            nonlocal close_calls
            close_calls += 1

        monkeypatch.setattr(loaded.corpus, "close", close_underlying)
        resource = _TrackingResource()
        golden = GoldenCorpus(
            spec=self._spec(),
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
