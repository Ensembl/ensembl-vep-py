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
"""Annotated-file summaries and corpus fact profiles."""

from __future__ import annotations

import logging
from typing import ClassVar

from attrs import asdict

from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.config import Config
from ensembl.vep.corpus.registry import CorpusOptionsError, build_corpus_from_options
from ensembl.vep.output.parser import parser_for
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.registry import register_runner
from ensembl.vep.stats import Reportable, RunSummaryCollector, render_reportable
from ensembl.vep.stats.profile import trace_run

__all__ = ["Runner"]

logger = logging.getLogger(__name__)


@register_runner("stats")
class Runner:
    """Runner for run summaries and corpus fact profiles.

    ``run`` produces the summary as a by-product; this command re-reads one from what a run wrote
    (any tool's CSQ/ANN, not only ours), or -- with ``--profile`` -- turns tracing on over a corpus
    and reports what each fact cost. The profile is opt-in because it is measurably slower.

    Attributes:
        name (ClassVar[str]): The registry name.
        sections (ClassVar[tuple[str, ...]]): The configuration sections read.
        description (ClassVar[str]): One-line description for the subcommand list.
    """

    name: ClassVar[str] = "stats"
    sections: ClassVar[tuple[str, ...]] = ("logging", "source", "corpus", "annotation", "stats")
    description: ClassVar[str] = "summarise an annotated file, or profile the facts over a corpus"

    def __init__(self, config: Config) -> None:
        """Hold the configuration."""
        self.config = config

    def run(self) -> ExitCode:
        """Build the requested report and write it."""
        stats = self.config.stats
        input_file = self.config.source.input_file
        if stats.profile and input_file is not None:
            logger.error("--profile measures a generated corpus; give it or --input_file, not both")
            return ExitCode.USAGE
        if not stats.profile and input_file is None:
            logger.error("give --input_file to summarise an annotated file, or --profile to profile a corpus")
            return ExitCode.USAGE

        try:
            report = self._profile() if stats.profile else self._summary()
        except CorpusOptionsError as problem:
            logger.error("%s", problem)
            return ExitCode.USAGE
        except (OSError, NotImplementedError, ValueError) as problem:
            logger.error("%s", problem)
            return ExitCode.FAILURE

        self._emit(render_reportable(report, stats.output_format))
        return ExitCode.OK

    def _summary(self) -> Reportable:
        """Re-read the annotated input into a run summary."""
        source = self.config.source
        assert source.input_file is not None
        reader = parser_for(source.input_file, source.format)
        collector = RunSummaryCollector(
            metadata={"input": str(source.input_file), "reader": type(reader).__name__}
        )
        for variation in reader.parse(source.input_file):
            collector.observe(variation)
        summary = collector.finish()
        logger.info("summarised %s: %s", source.input_file, summary.subtitle)
        return summary

    def _profile(self) -> Reportable:
        """Trace every overlap of the configured corpus and fold the derivations."""
        annotation = self.config.annotation
        corpus = self.config.corpus
        built_corpus = build_corpus_from_options(corpus)
        annotator = build_annotator(
            **asdict(annotation),
            upstream_distance=self.config.source.upstream_distance,
            downstream_distance=self.config.source.downstream_distance,
        )
        logger.info("profiling %s over %s", annotation.annotator, built_corpus.description)
        cost_table = "measured (fact_costs.json)" if annotation.fact_cost_override else "declared tiers"
        try:
            return trace_run(
                annotator,
                [overlap for _name, overlap in built_corpus.named_overlaps],
                annotator_name=annotation.annotator,
                cost_table=cost_table,
            )
        finally:
            built_corpus.close()

    def _emit(self, rendered: str) -> None:
        """Write the report to the configured destination."""
        destination = self.config.stats.output_file
        if destination is None:
            print(rendered)
            return
        destination.write_text(rendered, encoding="utf-8")
        logger.info("wrote %s", destination)
