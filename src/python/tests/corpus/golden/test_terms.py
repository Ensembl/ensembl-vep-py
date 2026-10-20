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
from ensembl.vep.corpus import build_corpus
from ensembl.vep.corpus.golden import (
    FORMAT_VERSION,
    GoldenReport,
    GoldenTerms,
    overlap_key,
    record_golden,
)
from ensembl.vep.corpus.spec import CorpusSpec

GOLDEN_FILE = Path(__file__).parents[1] / "data" / "golden_linear.json"
"""The checked-in golden terms the regression gate reads. Re-record with
``vep corpus golden --golden_file <this> --update_golden`` plus the spec in its header."""


@pytest.fixture(name="recorded", scope="module")
def fixture_recorded() -> GoldenTerms:
    """Linear's terms over a small corpus, recorded once."""
    corpus = build_corpus(CorpusSpec(kind="synthetic.all", fuzz=5, step=600))
    return record_golden(corpus, build_annotator(annotator="consequence.linear"), "consequence.linear")


class TestRecording:
    """What a recording holds."""

    def test_every_overlap_is_filed_under_its_feature_and_name(self, recorded: GoldenTerms) -> None:
        """One entry per overlap, keyed feature::name."""
        corpus = build_corpus(CorpusSpec(kind="synthetic.all", fuzz=5, step=600))
        expected_keys = {overlap_key(name, overlap) for name, overlap in corpus.named_overlaps}
        assert set(recorded.terms) == expected_keys
        assert len(recorded.terms) == len(corpus.named_overlaps)

    def test_it_names_its_corpus_annotator_and_catalogue(self, recorded: GoldenTerms) -> None:
        """Everything a later check needs to explain a mismatch."""
        assert recorded.spec == CorpusSpec(kind="synthetic.all", fuzz=5, step=600)
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
        corpus = build_corpus(golden.spec)
        recorded = record_golden(corpus, build_annotator(annotator=golden.annotator), golden.annotator)
        report = GoldenReport.compare(golden, recorded)
        assert report.ok, "\n".join(report.lines())
