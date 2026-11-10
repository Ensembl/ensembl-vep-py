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
"""Strategy and cache benchmark matrix."""

from __future__ import annotations

import logging
from typing import ClassVar

from ensembl.vep.bench.baseline import BenchResults
from ensembl.vep.bench.harness import consequence_strategies, native_kernels_available
from ensembl.vep.bench.suites import MatrixRow, build_configurations, run_matrix
from ensembl.vep.config import Config
from ensembl.vep.corpus.registry import CorpusOptionsError
from ensembl.vep.render.table import TableFormat, render_table
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.bench import common
from ensembl.vep.runner.registry import register_runner

__all__ = ["Runner"]

logger = logging.getLogger(__name__)

HEADERS = ("Configuration", "ms/1k ctx", "batch ms/1k ctx", "facts/ctx", "cost/ctx", "vs linear")
EXTRAS_HEADERS = ("Configuration", "Instrumentation")


@register_runner("bench.matrix")
class Runner:
    """Benchmark runner for strategy and cache combinations.

    ``cost/ctx`` and wall clock rank strategies differently on purpose: the cost is the number that
    survives a port to another language, while wall clock here is dominated by interpreter dispatch.

    Attributes:
        name (ClassVar[str]): The registry name.
        sections (ClassVar[tuple[str, ...]]): The configuration sections read.
        description (ClassVar[str]): One-line description for the subcommand list.
    """

    name: ClassVar[str] = "bench.matrix"
    sections: ClassVar[tuple[str, ...]] = ("logging", "corpus", "bench", "bench_matrix")
    description: ClassVar[str] = "time every strategy and cache configuration over a corpus"

    def __init__(self, config: Config) -> None:
        """Hold the configuration."""
        self.config = config

    def run(self) -> ExitCode:
        """Gate, measure, report, and record."""
        matrix = self.config.bench_matrix
        try:
            strategies = common.resolve_strategies(matrix.strategies, consequence_strategies())
            configurations = build_configurations(strategies, matrix.caches, matrix.batch_size)
            table_format = TableFormat(common.table_format(self.config))
        except ValueError as problem:
            logger.error("%s", problem)
            return ExitCode.USAGE
        if not native_kernels_available():
            logger.info(
                "consequence.columnar-kernels is not registered here: the native extension is not built"
            )

        try:
            with common.corpora(self.config) as labelled:
                if not self.config.bench.no_gate and not common.gate(labelled):
                    return ExitCode.FAILURE

                rows: list[MatrixRow] = []
                for label, overlaps in labelled:
                    corpus_rows = run_matrix(label, overlaps, configurations, common.policy(self.config))
                    print(f"\n{label}: {len(overlaps)} overlaps, {len(configurations)} configurations")
                    print(render_table(HEADERS, self._rows(corpus_rows), table_format))
                    extras = [
                        (row.configuration, _format_extras(row.extras)) for row in corpus_rows if row.extras
                    ]
                    if extras:
                        print(render_table(EXTRAS_HEADERS, extras, table_format))
                    rows.extend(corpus_rows)
        except CorpusOptionsError as unusable:
            logger.error("%s", unusable)
            return ExitCode.USAGE

        results = BenchResults(
            info=common.run_info(self.config, self.name), rows=tuple(row.to_json() for row in rows)
        )
        return common.record_results(self.config, results)

    @staticmethod
    def _rows(corpus_rows: list[MatrixRow]) -> list[tuple[object, ...]]:
        """Return a corpus's rows formatted, with the speed-up against the first configuration."""
        reference = corpus_rows[0].ms_per_1k if corpus_rows else 0.0
        formatted: list[tuple[object, ...]] = []
        for row in corpus_rows:
            speedup = reference / row.ms_per_1k if row.ms_per_1k else float("inf")
            formatted.append(
                (
                    row.configuration,
                    f"{row.ms_per_1k:.3f}",
                    f"{row.batch_ms_per_1k:.3f}",
                    f"{row.facts_per_context:.2f}",
                    f"{row.cost_per_context:.2f}",
                    f"{speedup:.2f}x",
                )
            )
        return formatted


def _format_extras(extras: dict[str, object]) -> str:
    """Return instrumentation as one line of key=value pairs."""
    parts = []
    for key, value in extras.items():
        parts.append(f"{key}={value:.4f}" if isinstance(value, float) else f"{key}={value}")
    return ", ".join(parts)
