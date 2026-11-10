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
"""Tests for the bench runners."""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the annotators the benchmarks run.
from ensembl.vep.bench.suites import Configuration
from ensembl.vep.bench.workspace import Workspace
from ensembl.vep.cli import main
from ensembl.vep.config import Config
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.bench import endtoend as endtoend_runner
from ensembl.vep.runner.bench.endtoend import Runner

SMALL = [
    "--kind",
    "synthetic.default",
    "--fuzz",
    "0",
    "--repeats",
    "1",
    "--warmup",
    "0",
]
PAIR = ["--strategies", "consequence.linear,consequence.staged", "--caches", "none"]


class TestMatrix:
    """``vep bench matrix``."""

    def test_it_times_the_requested_configurations(self, capsys: pytest.CaptureFixture[str]) -> None:
        """Each strategy x cache set gets a row, with the speed-up column."""
        arguments = ["bench", "matrix", *SMALL, "--strategies", "consequence.linear,consequence.staged"]
        assert main([*arguments, "--caches", "none,zonemap"]) == ExitCode.OK
        out = capsys.readouterr().out
        assert "consequence.staged +zonemap" in out
        assert "vs linear" in out

    def test_a_run_against_its_own_results_is_fine(self, tmp_path: Path) -> None:
        """Written results serve as a baseline; a generous threshold absorbs timing noise."""
        results = tmp_path / "results.json"
        assert (
            main(["bench", "matrix", *SMALL, *PAIR, "--no_gate", "--results_file", str(results)])
            == ExitCode.OK
        )
        assert json.loads(results.read_text(encoding="utf-8"))["info"]["verb"] == "bench.matrix"
        rerun = ["bench", "matrix", *SMALL, *PAIR, "--no_gate", "--baseline_file", str(results)]
        assert main([*rerun, "--regression_threshold", "100"]) == ExitCode.OK

    def test_a_regression_fails(self, tmp_path: Path) -> None:
        """A baseline a thousand times faster makes this run a regression."""
        results = tmp_path / "results.json"
        main(["bench", "matrix", *SMALL, *PAIR, "--no_gate", "--results_file", str(results)])
        payload = json.loads(results.read_text(encoding="utf-8"))
        for row in payload["rows"]:
            row["ms_per_1k"] /= 1000
        results.write_text(json.dumps(payload), encoding="utf-8")
        assert main(["bench", "matrix", *SMALL, *PAIR, "--no_gate", "--baseline_file", str(results)]) == (
            ExitCode.FAILURE
        )

    def test_a_baseline_over_another_corpus_is_refused(self, tmp_path: Path) -> None:
        """Comparing across corpora would mean nothing."""
        results = tmp_path / "results.json"
        main(["bench", "matrix", *SMALL, *PAIR, "--no_gate", "--results_file", str(results)])
        other = [
            "bench",
            "matrix",
            "--kind",
            "synthetic.default",
            "--fuzz",
            "3",
            "--repeats",
            "1",
            "--warmup",
            "0",
            *PAIR,
        ]
        assert main([*other, "--no_gate", "--baseline_file", str(results)]) == ExitCode.USAGE

    def test_an_unknown_strategy_or_cache_is_a_usage_error(self) -> None:
        """Named on the command line, so the caller's mistake."""
        assert main(["bench", "matrix", *SMALL, "--strategies", "consequence.nope"]) == ExitCode.USAGE
        assert main(["bench", "matrix", *SMALL, "--caches", "bogus"]) == ExitCode.USAGE


class TestSweep:
    """``vep bench sweep``."""

    def test_it_sweeps_the_batch_strategies(self, capsys: pytest.CaptureFixture[str]) -> None:
        """One column per strategy, one row per value."""
        assert main(["bench", "sweep", *SMALL, "--sweep_values", "32,128"]) == ExitCode.OK
        assert "consequence.columnar-staged ms/1k ctx" in capsys.readouterr().out

    def test_an_unknown_strategy_is_a_usage_error(self) -> None:
        """Named on the command line, so the caller's mistake."""
        assert main(["bench", "sweep", *SMALL, "--strategies", "consequence.nope"]) == ExitCode.USAGE


class TestOracle:
    """``vep bench oracle``."""

    def test_without_perl_it_still_compares_our_strategies(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """A missing Perl side is a note, not a failure: the oracle never gates."""
        arguments = ["bench", "oracle", "--repeats", "1", "--oracle_number", "2"]
        assert main([*arguments, "--perl_modules", str(tmp_path / "absent")]) == ExitCode.OK
        out = capsys.readouterr().out
        assert "coding_exact" in out
        assert "unavailable" in out


class TestEndToEnd:
    """``vep bench endtoend``."""

    def test_an_unknown_strategy_is_a_usage_error(self) -> None:
        """Configurations are resolved before any corpus is written."""
        arguments = ["bench", "endtoend", *SMALL, "--strategies", "consequence.nope"]
        assert main(arguments) == ExitCode.USAGE

    def test_missing_indexing_tools_are_a_usage_error(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Without htslib the corpus cannot be indexed, and an unindexed run would measure
        a different code path, so the verb stops instead."""
        monkeypatch.setattr(shutil, "which", lambda name: None)
        arguments = ["bench", "endtoend", *SMALL, *PAIR, "--work_dir", str(tmp_path)]

        assert main(arguments) == ExitCode.USAGE
        assert "htslib" in caplog.text

    def test_the_reference_is_linear_wherever_it_appears(self) -> None:
        """The gate compares against the reference strategy, not against whatever ran first."""
        runs = [
            (Configuration(strategy="consequence.columnar"), 1),
            (Configuration(strategy="consequence.linear"), 1),
        ]

        assert Runner._reference_key(runs) == ("consequence.linear", 1)  # pylint: disable=protected-access

    def test_without_linear_the_first_run_is_the_reference(self) -> None:
        """Timing a subset that leaves out the reference still gates against something."""
        runs = [(Configuration(strategy="consequence.columnar", caches=frozenset({"zonemap"})), 2)]

        assert Runner._reference_key(runs) == (  # pylint: disable=protected-access
            "consequence.columnar +zonemap",
            2,
        )

    @pytest.mark.integration
    def test_the_gate_compares_the_outputs_of_every_configuration(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Each configuration is run once and checked against the reference before timing."""
        arguments = ["bench", "endtoend", *SMALL, *PAIR, "--work_dir", str(tmp_path)]

        assert main(arguments) == ExitCode.OK
        assert "output gate passed" in caplog.text

    @pytest.mark.integration
    def test_it_times_each_configuration(self, capsys: pytest.CaptureFixture[str], tmp_path: Path) -> None:
        """Every configuration gets a row, measured as a whole run of the command line."""
        arguments = [
            "bench", "endtoend", *SMALL, "--strategies", "consequence.columnar",
            "--caches", "none", "--no_gate", "--work_dir", str(tmp_path),
        ]  # fmt: skip

        assert main(arguments) == ExitCode.OK
        out = capsys.readouterr().out
        assert "ms/1k variants" in out
        assert "consequence.columnar" in out

    @pytest.mark.integration
    def test_cleaning_leaves_the_directory_as_it_was(self, tmp_path: Path) -> None:
        """What the run created is removed; what was already there is not."""
        kept = tmp_path / "kept.txt"
        kept.write_text("mine", encoding="utf-8")
        arguments = [
            "bench", "endtoend", *SMALL, "--strategies", "consequence.columnar",
            "--caches", "none", "--no_gate", "--work_dir", str(tmp_path),
        ]  # fmt: skip

        assert main(arguments) == ExitCode.OK
        assert list(tmp_path.iterdir()) == [kept]

    @pytest.mark.integration
    def test_no_clean_keeps_the_corpus_and_the_outputs(self, tmp_path: Path) -> None:
        """A kept workspace can be re-run against, or inspected after a surprising result."""
        arguments = [
            "bench", "endtoend", *SMALL, "--strategies", "consequence.columnar",
            "--caches", "none", "--no_gate", "--no-clean", "--work_dir", str(tmp_path),
        ]  # fmt: skip

        assert main(arguments) == ExitCode.OK
        assert (tmp_path / "corpus.vcf").exists()
        assert (tmp_path / "corpus.gff3.bgz.tbi").exists()


class TestEndToEndFailures:
    """``vep bench endtoend`` when something goes wrong."""

    @staticmethod
    def _arguments(tmp_path: Path) -> list[str]:
        """Return a two-configuration run writing into a directory of its own."""
        return ["bench", "endtoend", *SMALL, *PAIR, "--work_dir", str(tmp_path)]

    def test_a_temporary_directory_is_used_when_no_directory_was_given(self) -> None:
        """Nothing is written into the working directory unless one was asked for."""
        directory = Runner(Config())._directory()  # pylint: disable=protected-access

        assert directory.parent == Path(tempfile.gettempdir())
        assert directory.name.startswith("vep-bench-endtoend-")

    def test_failed_indexing_is_reported(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """A bgzip or tabix that fails stops the benchmark rather than being annotated around."""

        def fail(_self: Workspace) -> None:
            raise subprocess.CalledProcessError(1, ["bgzip"])

        monkeypatch.setattr(Workspace, "index", fail)

        assert main(self._arguments(tmp_path)) == ExitCode.FAILURE
        assert "indexing the corpus failed" in caplog.text

    def test_a_disagreeing_configuration_is_not_timed(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """The gate reports the mismatching variants and the run stops before any timing."""
        monkeypatch.setattr(Workspace, "index", lambda _self: None)
        monkeypatch.setattr(endtoend_runner, "run_once", lambda command, log_path: None)
        monkeypatch.setattr(
            endtoend_runner, "output_mismatches", lambda reference, candidate: ["1:1 A/G: differs"]
        )
        monkeypatch.setattr(
            endtoend_runner,
            "time_runs",
            lambda *args, **kwargs: pytest.fail("timed a configuration that failed the gate"),
        )

        assert main(self._arguments(tmp_path)) == ExitCode.FAILURE
        assert "output gate failed" in caplog.text
        assert "1:1 A/G: differs" in caplog.text

    def test_a_failing_run_is_reported_with_its_log(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """The run's own stderr is what explains the failure."""

        def fail(*_args: object, **_kwargs: object) -> None:
            raise subprocess.CalledProcessError(2, ["vep"], stderr="it did not start")

        monkeypatch.setattr(Workspace, "index", lambda _self: None)
        monkeypatch.setattr(endtoend_runner, "time_runs", fail)

        assert main([*self._arguments(tmp_path), "--no_gate"]) == ExitCode.FAILURE
        assert "it did not start" in caplog.text
