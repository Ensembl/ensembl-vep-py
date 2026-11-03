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
"""Case-by-case comparison with the Perl ``VariationEffect.pm`` oracle."""

from __future__ import annotations

import logging
import subprocess
from typing import ClassVar

from ensembl.vep.annotators.consequence.rules.rulesets import build_fact_ruleset
from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.bench.baseline import BenchResults
from ensembl.vep.bench.harness import KERNELS_STRATEGY, REFERENCE_STRATEGY
from ensembl.vep.bench.oracle import (
    covered_cases,
    covered_through_subsumption,
    oracle_cases,
    perl_seconds,
    perl_terms,
    perl_unavailable_reason,
    python_seconds,
    python_terms,
)
from ensembl.vep.config import Config
from ensembl.vep.render.table import TableFormat, render_table
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.bench import common
from ensembl.vep.runner.registry import register_runner

__all__ = ["Runner"]

logger = logging.getLogger(__name__)

NO_TERMS = "(none)"


@register_runner("bench.oracle")
class Runner:
    """Oracle benchmark runner for Python and Perl strategies.

    A disagreement is reported, never gated on: the oracle is the external reference, and the Perl
    script covers only four categories, so the question asked is whether our terms cover its terms.

    Attributes:
        name (ClassVar[str]): The registry name.
        sections (ClassVar[tuple[str, ...]]): The configuration sections read.
        description (ClassVar[str]): One-line description for the subcommand list.
    """

    name: ClassVar[str] = "bench.oracle"
    sections: ClassVar[tuple[str, ...]] = ("logging", "bench", "bench_oracle")
    description: ClassVar[str] = "compare our strategies with the Perl VariationEffect.pm oracle"

    def __init__(self, config: Config) -> None:
        """Hold the configuration."""
        self.config = config

    def run(self) -> ExitCode:
        """Compare, time, report, and record."""
        options = self.config.bench_oracle
        try:
            strategies = common.resolve_strategies(options.strategies, (REFERENCE_STRATEGY, KERNELS_STRATEGY))
            table_format = TableFormat(common.table_format(self.config))
        except ValueError as problem:
            logger.error("%s", problem)
            return ExitCode.USAGE

        cases = oracle_cases()
        annotators = {name: build_annotator(annotator=name) for name in strategies}
        ours = {name: python_terms(annotator, cases) for name, annotator in annotators.items()}

        reason = perl_unavailable_reason(options.perl_modules, options.perl_oracle_script)
        perl: tuple[tuple[str, ...], ...] | None = None
        if reason is not None:
            logger.warning("Perl oracle unavailable, comparing our strategies only: %s", reason)
        else:
            assert options.perl_modules is not None
            try:
                perl = perl_terms(options.perl_modules, options.perl_oracle_script, len(cases))
            except (ValueError, OSError, subprocess.SubprocessError) as problem:
                logger.error("the Perl oracle failed: %s", problem)
                return ExitCode.FAILURE

        headers = ("case", "perl (oracle)", *strategies)
        rows = []
        for index, case in enumerate(cases):
            perl_cell = ",".join(perl[index]) or NO_TERMS if perl is not None else "unavailable"
            rows.append(
                (case.name, perl_cell, *(",".join(ours[name][index]) or NO_TERMS for name in strategies))
            )
        print(render_table(headers, rows, table_format))
        if perl is not None:
            subsumes = {rule.term: rule.subsumes for rule in build_fact_ruleset().rules}
            for name in strategies:
                logger.info(
                    "%s covers the Perl-reported terms in %d/%d cases directly, "
                    "%d/%d counting the terms its emitted ones suppress by subsumption",
                    name,
                    covered_cases(perl, ours[name]),
                    len(cases),
                    covered_through_subsumption(perl, ours[name], subsumes),
                    len(cases),
                )

        timings = self._timings(annotators, cases, perl is not None)
        print(
            render_table(
                ("subject", "ns/case (best)"),
                [(subject, f"{value:.1f}") for subject, value in timings],
                table_format,
            )
        )
        results = BenchResults(
            info=common.run_info(self.config, self.name),
            rows=tuple(
                {"corpus": "oracle", "configuration": subject, "ns_per_case": value}
                for subject, value in timings
            ),
        )
        return common.record_results(self.config, results)

    def _timings(self, annotators: dict, cases: tuple, with_perl: bool) -> list[tuple[str, float]]:
        """Return best nanoseconds per case for each strategy and, when available, for Perl."""
        options = self.config.bench_oracle
        repeats = self.config.bench.repeats
        number = options.oracle_number
        per_case = number * len(cases)
        timings = []
        for name, annotator in annotators.items():
            best = min(python_seconds(annotator, cases, number=number, repeats=repeats))
            timings.append((f"python {name}", best / per_case * 1e9))
        if with_perl and options.perl_modules is not None:
            best = min(
                perl_seconds(options.perl_modules, options.perl_oracle_script, number=number, repeats=repeats)
            )
            timings.append(("perl (includes interpreter start-up)", best / per_case * 1e9))
        return timings
