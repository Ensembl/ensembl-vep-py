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
"""Tests for results files and baseline comparison."""

from __future__ import annotations

from pathlib import Path

import pytest

from ensembl.vep.bench.baseline import (
    FORMAT_VERSION,
    BenchResults,
    RunInfo,
    compare_results,
    git_state,
    machine_info,
)


def _results(ms: float = 10.0, cost: float = 5.0, parameters: str = "abc") -> BenchResults:
    """Return a one-row result set."""
    row = {
        "corpus": "default",
        "configuration": "consequence.linear",
        "ms_per_1k": ms,
        "cost_per_context": cost,
    }
    return BenchResults(info=RunInfo(verb="bench.matrix", corpus_parameters=parameters), rows=(row,))


class TestFile:
    """Writing and reading results."""

    def test_it_round_trips(self, tmp_path: Path) -> None:
        """Written then loaded is equal."""
        path = tmp_path / "results.json"
        _results().write(path)
        assert BenchResults.load(path) == _results()

    def test_another_format_version_is_rejected(self) -> None:
        """An old file is refused, not misread."""
        text = _results().to_json().replace(f'"format_version": {FORMAT_VERSION}', '"format_version": 999')
        with pytest.raises(ValueError, match="999"):
            BenchResults.from_json(text)


class TestCompare:
    """Judging a run against a baseline."""

    def test_the_same_numbers_are_fine(self) -> None:
        """Nothing regressed, nothing changed."""
        comparison = compare_results(_results(), _results(), threshold=0.05)
        assert comparison.ok
        assert not comparison.regressions

    def test_a_slowdown_past_the_threshold_regresses(self) -> None:
        """Twenty percent slower against a five percent threshold."""
        comparison = compare_results(_results(ms=10.0), _results(ms=12.0), threshold=0.05)
        assert not comparison.ok
        (regression,) = comparison.regressions
        assert regression.metric == "ms_per_1k"
        assert regression.change == pytest.approx(0.2)

    def test_a_slowdown_within_the_threshold_does_not(self) -> None:
        """Noise is not a regression."""
        assert compare_results(_results(ms=10.0), _results(ms=10.4), threshold=0.05).ok

    def test_a_cost_change_is_reported_but_not_a_regression(self) -> None:
        """Cost is deterministic, so any change is behaviour, not timing."""
        comparison = compare_results(_results(cost=5.0), _results(cost=6.0), threshold=0.05)
        assert comparison.ok
        assert comparison.cost_changes

    def test_a_different_corpus_makes_the_comparison_meaningless(self) -> None:
        """Refused, rather than reported as a regression or an improvement."""
        comparison = compare_results(_results(parameters="abc"), _results(parameters="xyz"), threshold=0.05)
        assert comparison.corpus_changed
        assert not comparison.ok

    def test_rows_without_a_baseline_are_listed(self) -> None:
        """A new configuration has nothing to be compared against."""
        current = BenchResults(
            info=RunInfo(verb="bench.matrix", corpus_parameters="abc"),
            rows=({"corpus": "default", "configuration": "new", "ms_per_1k": 1.0},),
        )
        assert compare_results(_results(), current, threshold=0.05).unmatched == ("default / new",)


class TestContext:
    """What a results file records about where it came from."""

    def test_machine_and_git_state(self) -> None:
        """Enough to explain a different number."""
        assert {"system", "machine", "python", "cpus"} <= set(machine_info())
        assert set(git_state()) == {"revision", "dirty"}
