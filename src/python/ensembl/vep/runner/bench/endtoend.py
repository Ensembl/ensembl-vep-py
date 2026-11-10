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
"""End-to-end command-line benchmarks over exported corpora."""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import ClassVar

from ensembl.vep.bench.baseline import BenchResults
from ensembl.vep.bench.endtoend import (
    EndToEndRow,
    output_mismatches,
    run_command,
    run_once,
    summarise,
    time_runs,
)
from ensembl.vep.bench.harness import REFERENCE_STRATEGY, consequence_strategies
from ensembl.vep.bench.suites import Configuration, build_configurations
from ensembl.vep.bench.workspace import MissingTool, Workspace
from ensembl.vep.config import Config
from ensembl.vep.corpus.registry import build_corpus
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.render.table import TableFormat, render_table
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.bench import common
from ensembl.vep.runner.registry import register_runner

__all__ = ["Runner"]

logger = logging.getLogger(__name__)

HEADERS = ("Configuration", "processes", "s", "ms/1k variants", "variants/s", "peak RSS MB", "vs reference")
"""Columns of the results table. Peak RSS is the largest single process, not the tree's sum."""

_WORK_DIR_PREFIX = "vep-bench-endtoend-"


@register_runner("bench.endtoend")
class Runner:
    """End-to-end benchmark runner for full CLI configurations.

    The gate runs first and compares outputs rather than in-memory terms: a configuration whose
    annotated VCF disagrees with the reference's is not timed. Each run is a separate process given
    this working copy on ``PYTHONPATH``, so what is measured is the tree the benchmark was started
    from.

    Attributes:
        name (ClassVar[str]): The registry name.
        sections (ClassVar[tuple[str, ...]]): The configuration sections read.
        description (ClassVar[str]): One-line description for the subcommand list.
    """

    name: ClassVar[str] = "bench.endtoend"
    sections: ClassVar[tuple[str, ...]] = ("logging", "corpus", "bench", "bench_endtoend")
    description: ClassVar[str] = "time whole command-line runs over a corpus on disk"

    def __init__(self, config: Config) -> None:
        """Hold the configuration."""
        self.config = config

    def run(self) -> ExitCode:
        """Prepare the corpus, gate, measure, report, and clean up after itself."""
        endtoend = self.config.bench_endtoend
        try:
            strategies = common.resolve_strategies(endtoend.strategies, consequence_strategies())
            configurations = build_configurations(strategies, endtoend.caches, endtoend.batch_size)
            common.table_format(self.config)
        except ValueError as problem:
            logger.error("%s", problem)
            return ExitCode.USAGE

        workspace = Workspace(directory=self._directory(), clean=endtoend.clean)
        try:
            return self._measure(workspace, configurations)
        finally:
            workspace.close()

    def _directory(self) -> Path:
        """Return where the corpus, its indexes and the run outputs go."""
        chosen = self.config.bench_endtoend.work_dir or self.config.corpus.corpus_dir
        if chosen is not None:
            return chosen
        return Path(tempfile.gettempdir()) / f"{_WORK_DIR_PREFIX}{os.getpid()}"

    def _measure(self, workspace: Workspace, configurations: list[Configuration]) -> ExitCode:
        """Run every configuration over the prepared corpus and report the timings."""
        spec = CorpusSpec.from_options(self.config.corpus)
        workspace.open()
        try:
            with build_corpus(spec) as corpus:
                workspace.materialise(corpus)
            workspace.index()
        except MissingTool as missing:
            logger.error("%s", missing)
            return ExitCode.USAGE
        except subprocess.CalledProcessError as failure:
            logger.error("indexing the corpus failed: %s", failure)
            return ExitCode.FAILURE

        variants = workspace.variant_count()
        runs = [
            (configuration, processes)
            for configuration in configurations
            for processes in self.config.bench_endtoend.processes
        ]

        try:
            if not self.config.bench.no_gate and not self._gate(workspace, runs):
                return ExitCode.FAILURE
            rows = self._time(workspace, runs, corpus=self.config.corpus.kind, variants=variants)
        except subprocess.CalledProcessError as failure:
            logger.error("a run failed (%s):\n%s", failure.returncode, failure.stderr)
            return ExitCode.FAILURE

        print(f"\n{self.config.corpus.kind}: {variants} variants, {len(runs)} runs")
        print(render_table(HEADERS, self._rows(rows), TableFormat(common.table_format(self.config))))

        results = BenchResults(
            info=common.run_info(self.config, self.name), rows=tuple(row.to_json() for row in rows)
        )
        return common.record_results(self.config, results)

    def _gate(self, workspace: Workspace, runs: list[tuple[Configuration, int]]) -> bool:
        """Run each configuration once and check its output against the reference's.

        Args:
            workspace (Workspace): The prepared corpus.
            runs (list[tuple[Configuration, int]]): Each configuration with its worker count.

        Returns:
            bool: True when every configuration produced the reference's consequence terms.
        """
        outputs: dict[tuple[str, int], Path] = {}
        for configuration, processes in runs:
            outputs[(configuration.label, processes)] = self._run_for_gate(
                workspace, configuration, processes
            )

        reference_key = self._reference_key(runs)
        reference_output = outputs[reference_key]
        for key, output in outputs.items():
            if key == reference_key:
                continue
            mismatches = output_mismatches(reference_output, output)
            if mismatches:
                logger.error(
                    "output gate failed for %s on %d process(es), against %s:",
                    key[0],
                    key[1],
                    reference_key[0],
                )
                for mismatch in mismatches:
                    logger.error("  %s", mismatch)
                return False
        logger.info("output gate passed: %d configurations agree with %s", len(outputs), reference_key[0])
        return True

    @staticmethod
    def _reference_key(runs: list[tuple[Configuration, int]]) -> tuple[str, int]:
        """Return the run every other run's output is compared against."""
        for configuration, processes in runs:
            if configuration.strategy == REFERENCE_STRATEGY:
                return (configuration.label, processes)
        configuration, processes = runs[0]
        return (configuration.label, processes)

    def _run_for_gate(self, workspace: Workspace, configuration: Configuration, processes: int) -> Path:
        """Run one configuration once, keeping its output for the gate to compare."""
        output = self._output_path(workspace, configuration, processes)
        log_path = output.with_suffix(".log")
        workspace.track(output)
        workspace.track(log_path)
        command = run_command(configuration, processes=processes, workspace=workspace, output=output)
        run_once(command, log_path=log_path)
        return output

    def _time(
        self,
        workspace: Workspace,
        runs: list[tuple[Configuration, int]],
        *,
        corpus: str,
        variants: int,
    ) -> list[EndToEndRow]:
        """Time every run and return one row each."""
        bench = self.config.bench
        rows: list[EndToEndRow] = []
        for configuration, processes in runs:
            output = self._output_path(workspace, configuration, processes)
            log_path = output.with_suffix(".log")
            workspace.track(output)
            workspace.track(log_path)
            command = run_command(configuration, processes=processes, workspace=workspace, output=output)
            outcomes = time_runs(command, log_path=log_path, repeats=bench.repeats, warmup=bench.warmup)
            rows.append(
                summarise(
                    outcomes,
                    corpus=corpus,
                    configuration=configuration.label,
                    processes=processes,
                    variants=variants,
                )
            )
        return rows

    @staticmethod
    def _output_path(workspace: Workspace, configuration: Configuration, processes: int) -> Path:
        """Return where one configuration's annotated VCF goes."""
        stem = configuration.label.replace(" ", "").replace("+", "_").replace(".", "-")
        return workspace.directory / f"run_{stem}_p{processes}.vcf"

    @staticmethod
    def _rows(rows: list[EndToEndRow]) -> list[tuple[object, ...]]:
        """Return the rows formatted, with the speed-up against the first run."""
        reference = rows[0].seconds if rows else 0.0
        formatted: list[tuple[object, ...]] = []
        for row in rows:
            speedup = reference / row.seconds if row.seconds else float("inf")
            formatted.append(
                (
                    row.configuration,
                    row.processes,
                    f"{row.seconds:.2f}",
                    f"{row.ms_per_1k:.1f}",
                    f"{row.variants_per_second:.1f}",
                    f"{row.peak_rss_bytes / 1024 / 1024:.0f}",
                    f"{speedup:.2f}x",
                )
            )
        return formatted
