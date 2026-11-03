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
"""Consequence fact cost measurement."""

from __future__ import annotations

import logging
from typing import ClassVar

from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.facts.costs import (
    declared_vs_measured,
    normalised_cost_table,
    profile_self_nanoseconds,
    write_cost_table,
)
from ensembl.vep.config import Config
from ensembl.vep.config.report import resolve_report_format
from ensembl.vep.corpus.registry import CorpusOptionsError, build_corpus_from_options
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.render.table import TableFormat, render_table
from ensembl.vep.runner.registry import register_runner

__all__ = ["Runner"]

logger = logging.getLogger(__name__)

HEADERS = ("Fact", "Declared", "Measured ns", "Measured cost")


@register_runner("catalogue.costs")
class Runner:
    """Runner for profiling fact costs against declared tiers.

    Profiling replaces every fact's ``compute`` with a timing shim, so this is a deliberate
    measurement run and never something an annotation path does.
    """

    name: ClassVar[str] = "catalogue.costs"
    sections: ClassVar[tuple[str, ...]] = ("logging", "catalogue", "catalogue_costs", "corpus")
    description: ClassVar[str] = "measure each fact's real cost against its declared tier"

    def __init__(self, config: Config) -> None:
        """Hold the configuration."""
        self.config = config

    def run(self) -> ExitCode:
        """Profile, report, and optionally persist the measured table."""
        try:
            with build_corpus_from_options(self.config.corpus) as corpus:
                logger.info("profiling over %s", corpus.description)
                contexts = FactContext.for_overlaps([overlap for _name, overlap in corpus.named_overlaps])
                logger.info("%d contexts", len(contexts))
                measured = profile_self_nanoseconds(contexts)
        except CorpusOptionsError as unusable:
            logger.error("%s", unusable)
            return ExitCode.USAGE

        rows: list[tuple[object, ...]] = []
        for name, declared, nanoseconds, cost in declared_vs_measured(measured):
            rows.append((name, declared, f"{nanoseconds:.1f}", cost))
        try:
            table_format = resolve_report_format(
                self.config.catalogue_costs.output_format, graph=False, table_default="text", graph_formats=()
            )
        except ValueError as unusable:
            logger.error("%s", unusable)
            return ExitCode.USAGE
        print(render_table(HEADERS, rows, TableFormat(table_format)))

        if self.config.catalogue_costs.output_file is not None:
            destination = write_cost_table(
                normalised_cost_table(measured), self.config.catalogue_costs.output_file
            )
            logger.info("wrote measured cost table to %s", destination)
        return ExitCode.OK
