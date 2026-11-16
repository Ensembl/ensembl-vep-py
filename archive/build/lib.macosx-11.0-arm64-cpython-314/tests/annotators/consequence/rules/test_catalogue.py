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
"""Tests for the checked-in overlap consequence catalogue (schema v4)."""

from __future__ import annotations

from collections.abc import Iterator


import pytest

from ensembl.vep.annotators.consequence.rules import catalogue as consequences_module
from ensembl.vep.annotators.consequence.rules.catalogue import (
    DEFAULT_TERM,
    EXCLUDED_TERMS,
    load_consequence_definitions,
)
from ensembl.vep.annotators.consequence.rules.expressions import expr_from_json
from ensembl.vep.annotators.consequence.facts.base import REGISTRY


class TestCatalogueContents:
    """What the shipped catalogue holds, and the shape the loader expects of it."""

    def test_catalogue_terms_count(self) -> None:
        """Catalogue terms count."""
        assert len(load_consequence_definitions()) == 41

    def test_catalogue_terms_are_unique(self) -> None:
        """Catalogue terms are unique."""
        definitions = load_consequence_definitions()
        terms = [definition.term for definition in definitions]
        assert len(terms) == len(set(terms))

    def test_default_term_is_on_the_exclusion_list(self) -> None:
        """intergenic_variant is applied by the caller, not matched by a rule (see rulesets.prune)."""
        assert DEFAULT_TERM in EXCLUDED_TERMS

    def test_catalogue_ranks_are_unique(self) -> None:
        """Catalogue ranks are unique."""
        definitions = load_consequence_definitions()
        ranks = [definition.rank for definition in definitions]
        assert len(ranks) == len(set(ranks))


class TestWhenFormulas:
    """The ``when`` of each entry: which terms carry one, and what it may reference."""

    def test_only_the_excluded_terms_have_no_when(self) -> None:
        """Only the excluded terms have no when."""
        definitions = load_consequence_definitions()
        terms_without_when = {definition.term for definition in definitions if definition.when is None}
        assert terms_without_when == set(EXCLUDED_TERMS)

    def test_every_when_compiles_to_an_expression_over_known_atoms(self) -> None:
        """Every when compiles to an expression over known atoms."""
        definitions = load_consequence_definitions()
        covered = [definition for definition in definitions if definition.when is not None]
        # Every catalogue term except the two on EXCLUDED_TERMS -- see the module docstring of
        # consequence/rules/catalogue.py for what those two are and why each is excluded.
        assert len(covered) == 39
        for definition in covered:
            expression = expr_from_json(definition.when)
            unknown_atoms = expression.atoms() - set(REGISTRY)
            assert not unknown_atoms, f"{definition.term} references unknown facts: {unknown_atoms}"

    def test_subsumed_terms_referenced_exist_in_the_catalogue(self) -> None:
        """Subsumed terms referenced should exist in the catalogue."""
        definitions = load_consequence_definitions()
        known_terms = {definition.term for definition in definitions}
        for definition in definitions:
            if definition.when is None:
                continue
            unknown_subsumed = set(definition.subsumes) - known_terms
            assert not unknown_subsumed, f"{definition.term} subsumes unknown terms: {unknown_subsumed}"


class TestSchemaValidation:
    """Exercise load_consequence_definitions' error paths via a monkeypatched payload.

    load_consequence_definitions is memoised with functools.cache, so the cache is explicitly
    cleared before and after monkeypatching to avoid either reusing a stale real-catalogue
    result or leaking a broken cached result into later tests.
    """

    @pytest.fixture(autouse=True)
    def _reset_cache_around_each_test(self) -> Iterator[None]:
        """Reset cache around each test."""
        load_consequence_definitions.cache_clear()
        yield
        load_consequence_definitions.cache_clear()

    def test_rejects_unsupported_schema_version(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Rejects unsupported schema version."""
        monkeypatch.setattr(
            consequences_module.json, "loads", lambda _text: {"schema_version": 999, "consequences": []}
        )
        with pytest.raises(ValueError, match="Unsupported overlap consequence catalogue schema version"):
            load_consequence_definitions()

    def test_rejects_duplicate_terms(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Rejects duplicate terms."""
        duplicate_payload = {
            "schema_version": 5,
            "consequences": [
                {"term": "stop_gained", "accession": "SO:0001587", "rank": 1},
                {"term": "stop_gained", "accession": "SO:0001587", "rank": 2},
            ],
        }
        monkeypatch.setattr(consequences_module.json, "loads", lambda _text: duplicate_payload)
        with pytest.raises(ValueError, match="Duplicate consequence catalogue term"):
            load_consequence_definitions()

    def test_real_catalogue_reloads_cleanly_after_cache_clear(self) -> None:
        """Real catalogue reloads cleanly after cache clear."""
        assert len(load_consequence_definitions()) == 41
