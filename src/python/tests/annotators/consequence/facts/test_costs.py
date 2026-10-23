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
"""Tests for measuring fact costs."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the facts being profiled.
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.common import CostTier
from ensembl.vep.annotators.consequence.facts.base import REGISTRY
from ensembl.vep.annotators.consequence.facts.costs import (
    declared_vs_measured,
    fact_costs_path,
    normalised_cost_table,
    profile_self_nanoseconds,
    write_cost_table,
)
from ensembl.vep.corpus.synthetic.default import default_corpus


@pytest.fixture(name="measured", scope="module")
def fixture_measured() -> dict[str, float]:
    """Self-times profiled over a small corpus."""
    corpus = default_corpus(fuzz_count=0)
    overlaps = [overlap for _name, overlap in corpus.named_overlaps][:120]
    return profile_self_nanoseconds(FactContext.for_overlaps(overlaps))


class TestProfiler:
    """Timing every fact's own work."""

    def test_every_registered_fact_is_measured(self, measured: dict[str, float]) -> None:
        """A fact absent from the table would have no measured cost to compare."""
        assert set(measured) == set(REGISTRY)

    def test_self_times_are_positive(self, measured: dict[str, float]) -> None:
        """Self-time excludes nested facts, but is never negative."""
        assert all(value >= 0 for value in measured.values())

    def test_the_registry_is_restored(self, measured: dict[str, float]) -> None:
        """The profiler swaps every compute for a timing shim and must put them all back."""
        assert measured
        for name, registered_fact in REGISTRY.items():
            assert registered_fact.compute.__name__ != "timed_compute", name

    def test_a_failure_still_restores_the_registry(self) -> None:
        """The restore is in a finally block, so a bad context cannot leave shims installed."""
        with pytest.raises(Exception):  # pylint: disable=broad-exception-caught
            profile_self_nanoseconds([FactContext(overlap=object())])  # type: ignore[arg-type]
        for registered_fact in REGISTRY.values():
            assert registered_fact.compute.__name__ != "timed_compute"


class TestNormalisation:
    """Turning nanoseconds into the same units the declared tiers use."""

    def test_costs_are_whole_numbers(self, measured: dict[str, float]) -> None:
        """A cost is summed and compared as an integer everywhere else."""
        assert all(isinstance(value, int) for value in normalised_cost_table(measured).values())

    def test_costs_are_never_negative(self, measured: dict[str, float]) -> None:
        """Rounding a very fast fact must floor at zero, not go below it."""
        assert all(value >= 0 for value in normalised_cost_table(measured).values())

    def test_every_measured_fact_gets_a_cost(self, measured: dict[str, float]) -> None:
        """The table is what apply_fact_cost_overrides consumes."""
        assert set(normalised_cost_table(measured)) == set(measured)

    def test_normalising_without_a_bounds_fact_is_rejected(self) -> None:
        """BOUNDS = 1 is the unit; with none measured there is nothing to normalise against."""
        free_only = {name: 10.0 for name, fact in REGISTRY.items() if fact.cost == CostTier.FREE}
        with pytest.raises(RuntimeError, match="BOUNDS"):
            normalised_cost_table(free_only)


class TestComparison:
    """Declared tier against what was measured."""

    def test_one_row_per_measured_fact(self, measured: dict[str, float]) -> None:
        """The report shows every fact, not only the surprising ones."""
        assert len(declared_vs_measured(measured)) == len(measured)

    def test_rows_are_ordered_by_declared_cost(self, measured: dict[str, float]) -> None:
        """So the table reads cheapest first, as the tiers are meant to."""
        costs = [REGISTRY[row[0]].cost for row in declared_vs_measured(measured)]
        assert costs == sorted(costs)

    def test_the_declared_tier_is_named(self, measured: dict[str, float]) -> None:
        """A tier number alone would not say which tier it is."""
        assert declared_vs_measured(measured)[0][1] in {tier.name for tier in CostTier}


class TestWriting:
    """Persisting a measured table."""

    def test_the_default_destination_is_the_table_the_loader_reads(self) -> None:
        """These disagreed before: the writer wrote a path nothing ever loaded."""
        assert fact_costs_path().parent.name == "data"
        assert fact_costs_path().parent.parent.name == "consequence"

    def test_a_table_round_trips(self, tmp_path: Path) -> None:
        """What is written is what apply_fact_cost_overrides will read back."""
        destination = write_cost_table({"is_snv": 1, "overlap_cds": 20}, tmp_path / "costs.json")
        assert json.loads(destination.read_text(encoding="utf-8")) == {"is_snv": 1, "overlap_cds": 20}

    def test_keys_are_sorted(self, tmp_path: Path) -> None:
        """So a re-measurement diffs cleanly against the checked-in table."""
        destination = write_cost_table({"z": 1, "a": 2}, tmp_path / "costs.json")
        assert destination.read_text(encoding="utf-8").index('"a"') < destination.read_text(
            encoding="utf-8"
        ).index('"z"')
