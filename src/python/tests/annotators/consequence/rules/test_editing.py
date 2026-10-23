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
"""Tests for catalogue editing."""

from __future__ import annotations

import copy

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the facts the when-atom checks validate against.
from ensembl.vep.annotators.consequence.rules.catalogue import load_catalogue
from ensembl.vep.annotators.consequence.rules.editing import (
    accession_index,
    apply_add,
    apply_remove,
    build_entry,
    derive_subsumes,
    insert_at_rank,
    next_rank,
    remove_term,
)


@pytest.fixture(name="payload")
def fixture_payload() -> dict:
    """A private copy of the packaged catalogue, safe to mutate."""
    return copy.deepcopy(load_catalogue())


class TestAdd:
    """Adding a term."""

    def test_a_term_is_inserted(self, payload: dict) -> None:
        """The entry count grows by exactly one."""
        before = len(payload["consequences"])
        result = apply_add(payload, term="sample", accession="SO:0009999", impact="MODIFIER", when='"is_snv"')
        assert result is not None
        updated, _ = result
        assert len(updated["consequences"]) == before + 1

    def test_the_entry_carries_what_was_asked_for(self, payload: dict) -> None:
        """Rank, tier and impact all come from the caller."""
        result = apply_add(
            payload, term="sample", accession="SO:0009999", impact="MODIFIER", when='"is_snv"', rank=20
        )
        assert result is not None
        updated, _ = result
        entry = next(e for e in updated["consequences"] if e["term"] == "sample")
        assert (entry["accession"], entry["rank"], entry["impact"]) == ("SO:0009999", 20, "MODIFIER")

    def test_a_missing_accession_is_rejected(self, payload: dict) -> None:
        """An entry with no accession cannot be traced back to the ontology."""
        with pytest.raises(ValueError, match="accession"):
            apply_add(payload, term="sample", accession=None, impact="MODIFIER")

    def test_a_malformed_accession_is_rejected(self, payload: dict) -> None:
        """Catches a typo before it reaches the catalogue."""
        with pytest.raises(ValueError, match="SO:NNNNNNN"):
            apply_add(payload, term="sample", accession="SO:1", impact="MODIFIER")

    def test_a_duplicate_term_is_rejected(self, payload: dict) -> None:
        """Terms are the catalogue's key."""
        with pytest.raises(ValueError, match="already has a term"):
            apply_add(payload, term="intron_variant", accession="SO:0009999", impact="MODIFIER")

    def test_a_duplicate_accession_is_rejected(self, payload: dict) -> None:
        """Two terms sharing an accession would make the SO mapping ambiguous."""
        existing = payload["consequences"][0]["accession"]
        with pytest.raises(ValueError, match="already has an entry"):
            apply_add(payload, term="sample", accession=existing, impact="MODIFIER")

    def test_an_unknown_subsumes_target_is_rejected(self, payload: dict) -> None:
        """Suppressing a term that does not exist would silently do nothing."""
        with pytest.raises(ValueError, match="subsumes target"):
            apply_add(payload, term="sample", accession="SO:0009999", impact="MODIFIER", subsumes=("nope",))

    def test_an_entry_with_no_when_is_refused_by_default(self, payload: dict) -> None:
        """A wrong formula is the one error conformance cannot catch, so none is ever guessed."""
        assert apply_add(payload, term="sample", accession="SO:0009999", impact="MODIFIER") is None

    def test_an_entry_with_no_when_can_be_allowed_explicitly(self, payload: dict) -> None:
        """The placeholder is written only when the caller acknowledges it."""
        result = apply_add(
            payload, term="sample", accession="SO:0009999", impact="MODIFIER", allow_missing_when=True
        )
        assert result is not None


class TestRemove:
    """Removing a term."""

    def test_the_term_is_gone(self, payload: dict) -> None:
        """The entry count shrinks by exactly one."""
        before = len(payload["consequences"])
        updated, _messages = apply_remove(payload, term="feature_elongation")
        assert len(updated["consequences"]) == before - 1

    def test_it_is_cleaned_out_of_other_entries(self, payload: dict) -> None:
        """A dangling subsumes target would suppress nothing and read as if it did."""
        updated, _messages = apply_remove(payload, term="feature_elongation")
        for entry in updated["consequences"]:
            assert "feature_elongation" not in entry.get("subsumes", [])

    def test_an_unknown_term_is_rejected(self, payload: dict) -> None:
        """Silently succeeding would hide a typo."""
        with pytest.raises(ValueError, match="no catalogue entry"):
            apply_remove(payload, term="nope")


class TestRanks:
    """Rank arithmetic."""

    def test_the_next_rank_follows_the_highest(self, payload: dict) -> None:
        """Appending must not collide with an existing rank."""
        ranks = [entry["rank"] for entry in payload["consequences"]]
        assert next_rank(payload["consequences"]) == max(ranks) + 1

    def test_inserting_shifts_later_ranks_down(self, payload: dict) -> None:
        """Ranks are a total order, so an insertion has to make room."""
        entry = build_entry(
            term="sample",
            accession="SO:0009999",
            rank=0,
            tier=3,
            impact="MODIFIER",
            so_supported_parents=(),
            when="is_snv",
            subsumes=(),
            notes=None,
        )
        updated = insert_at_rank(payload["consequences"], entry, 5)
        ranks = [item["rank"] for item in updated]
        assert len(ranks) == len(set(ranks))
        assert next(item["rank"] for item in updated if item["term"] == "sample") == 5

    def test_appending_when_no_rank_is_given(self, payload: dict) -> None:
        """The common case: a new term goes last."""
        entry = build_entry(
            term="sample",
            accession="SO:0009999",
            rank=0,
            tier=3,
            impact="MODIFIER",
            so_supported_parents=(),
            when="is_snv",
            subsumes=(),
            notes=None,
        )
        updated = insert_at_rank(payload["consequences"], entry, None)
        assert updated[-1]["term"] == "sample"


class TestSubsumesDerivation:
    """Deriving suppression from the SO hierarchy."""

    def test_a_parent_accession_yields_its_term(self, payload: dict) -> None:
        """If the new term is_a an existing one, it is the more specific and suppresses it."""
        parent = next(e for e in payload["consequences"] if e.get("when"))
        derived, _notes = derive_subsumes([parent["accession"]], accession_index(payload["consequences"]))
        assert parent["term"] in derived

    def test_an_unknown_accession_derives_nothing(self, payload: dict) -> None:
        """A parent outside the catalogue has no term to suppress."""
        derived, notes = derive_subsumes(["SO:0000001"], accession_index(payload["consequences"]))
        assert not derived
        assert notes


class TestRemoveTerm:
    """The low-level removal helper."""

    def test_it_reports_where_it_cleaned_up(self, payload: dict) -> None:
        """So the caller can say what else changed."""
        _remaining, cleaned_from = remove_term(payload["consequences"], "feature_elongation")
        assert isinstance(cleaned_from, list)
