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
"""Tests for golden term serialisation and comparison."""

from __future__ import annotations

from pathlib import Path

import pytest

from ensembl.vep.corpus.golden import FORMAT_VERSION, GoldenReport, GoldenTerms
from ensembl.vep.corpus.spec import CorpusSpec


@pytest.fixture(name="recorded")
def fixture_recorded() -> GoldenTerms:
    """Return a small representative recording."""
    return GoldenTerms(
        spec=CorpusSpec(kind="golden.file", parameters={"dataset": "test"}),
        annotator="consequence.linear",
        catalogue="catalogue-1",
        terms={"TX1::variant-1": ("missense_variant",), "TX1::variant-2": ("intron_variant",)},
    )


class TestGoldenTerms:
    """Stable file representation."""

    def test_round_trips_with_one_overlap_per_line(self, recorded: GoldenTerms, tmp_path: Path) -> None:
        """Write deterministic, reviewable JSON and load it unchanged."""
        path = tmp_path / "terms.json"
        recorded.write(path)

        assert GoldenTerms.load(path) == recorded
        assert len([line for line in path.read_text(encoding="utf-8").splitlines() if "::" in line]) == 2

    def test_rejects_another_format_version(self, recorded: GoldenTerms) -> None:
        """Refuse an old file rather than misreading it."""
        text = recorded.to_json().replace(f'"format_version": {FORMAT_VERSION}', '"format_version": 999')
        with pytest.raises(ValueError, match="999"):
            GoldenTerms.from_json(text)


class TestGoldenReport:
    """Separate provenance and term mismatches."""

    def test_identical_terms_are_ok(self, recorded: GoldenTerms) -> None:
        """Report an exact match succinctly."""
        report = GoldenReport.compare(recorded, recorded)
        assert report.ok
        assert report.lines() == ["recorded terms match the golden file"]

    def test_changed_terms_and_memberships_are_reported(self, recorded: GoldenTerms) -> None:
        """Name moved, missing, and unexpected overlap terms."""
        changed = GoldenTerms(
            spec=recorded.spec,
            annotator=recorded.annotator,
            catalogue=recorded.catalogue,
            terms={"TX1::variant-1": ("frameshift_variant",), "TX2::variant-3": ()},
        )
        report = GoldenReport.compare(recorded, changed)

        assert not report.ok
        assert report.missing == ("TX1::variant-2",)
        assert report.unexpected == ("TX2::variant-3",)
        assert report.mismatches == (("TX1::variant-1", ("missense_variant",), ("frameshift_variant",)),)

    def test_another_annotator_is_informational(self, recorded: GoldenTerms) -> None:
        """Permit checking another strategy against the same pinned expectations."""
        changed = GoldenTerms(
            spec=recorded.spec,
            annotator="consequence.staged",
            catalogue=recorded.catalogue,
            terms=recorded.terms,
        )
        report = GoldenReport.compare(recorded, changed)
        assert report.annotator_changed
        assert report.ok
