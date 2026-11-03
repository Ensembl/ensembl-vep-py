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
"""Benchmarks across values of one numeric annotation setting."""

from __future__ import annotations

import logging
from typing import ClassVar

from ensembl.vep.bench.baseline import BenchResults
from ensembl.vep.bench.sweep import DEFAULT_BATCH_STRATEGIES, Sweep, SweepPoint, sweep_strategy
from ensembl.vep.config import Config
from ensembl.vep.corpus.registry import CorpusOptionsError
from ensembl.vep.render.table import TableFormat, render_table
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.bench import common
from ensembl.vep.runner.registry import register_runner

__all__ = ["Runner"]

logger = logging.getLogger(__name__)


@register_runner("bench.sweep")
class Runner:
    """Benchmark runner sweeping one setting across strategies.

    Attributes:
        name (ClassVar[str]): The registry name.
        sections (ClassVar[tuple[str, ...]]): The configuration sections read.
        description (ClassVar[str]): One-line description for the subcommand list.
    """

    name: ClassVar[str] = "bench.sweep"
    sections: ClassVar[tuple[str, ...]] = ("logging", "corpus", "bench", "bench_sweep")
    description: ClassVar[str] = "sweep one numeric annotation setting and time each value"

    def __init__(self, config: Config) -> None:
        """Hold the configuration."""
        self.config = config

    def run(self) -> ExitCode:
        """Sweep, report, and record."""
        options = self.config.bench_sweep
        try:
            strategies = common.resolve_strategies(options.strategies, DEFAULT_BATCH_STRATEGIES)
            sweep = Sweep(knob=options.sweep_knob, values=tuple(options.sweep_values))
            table_format = TableFormat(common.table_format(self.config))
        except ValueError as problem:
            logger.error("%s", problem)
            return ExitCode.USAGE
        if not strategies:
            logger.error("no strategy to sweep is registered here")
            return ExitCode.USAGE

        failed = False
        points: list[SweepPoint] = []
        try:
            with common.corpora(self.config) as labelled:
                for label, overlaps in labelled:
                    corpus_points: list[SweepPoint] = []
                    for strategy in strategies:
                        measured, failures = sweep_strategy(
                            label, strategy, overlaps, sweep, common.policy(self.config)
                        )
                        for failure in failures:
                            logger.error("%s", failure)
                        failed = failed or bool(failures)
                        corpus_points.extend(measured)
                    print(f"\n{label}: {len(overlaps)} overlaps, sweeping {sweep.knob}")
                    print(
                        render_table(
                            ("value", *(f"{name} ms/1k ctx" for name in strategies)),
                            _rows(sweep, strategies, corpus_points),
                            table_format,
                        )
                    )
                    points.extend(corpus_points)
        except CorpusOptionsError as unusable:
            logger.error("%s", unusable)
            return ExitCode.USAGE

        if failed:
            return ExitCode.FAILURE
        results = BenchResults(
            info=common.run_info(self.config, self.name), rows=tuple(point.to_json() for point in points)
        )
        return common.record_results(self.config, results)


def _rows(sweep: Sweep, strategies: list[str], points: list[SweepPoint]) -> list[tuple[object, ...]]:
    """Return one row per value, one column per strategy."""
    by_key = {(point.strategy, point.value): point.ms_per_1k for point in points}
    rows: list[tuple[object, ...]] = []
    for value in sweep.values:
        cells: list[object] = [value]
        for strategy in strategies:
            measured = by_key.get((strategy, value))
            cells.append(f"{measured:.3f}" if measured is not None else "-")
        rows.append(tuple(cells))
    return rows
