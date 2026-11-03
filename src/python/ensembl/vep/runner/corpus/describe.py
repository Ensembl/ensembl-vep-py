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
"""Registered corpus descriptions."""

from __future__ import annotations

import logging
from typing import ClassVar

from attrs import asdict

from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.config import Config
from ensembl.vep.config.report import resolve_report_format
from ensembl.vep.corpus.annotate import annotated_variations
from ensembl.vep.corpus.base import Corpus
from ensembl.vep.corpus.describe import (
    CHROMOSOME_HEADERS,
    CLASS_HEADERS,
    OVERLAP_HEADERS,
    TRANSCRIPT_HEADERS,
    CorpusDescription,
)
from ensembl.vep.corpus.registry import CorpusOptionsError, build_corpus_from_options
from ensembl.vep.render.page import Section, render_page
from ensembl.vep.render.table import TableFormat, render_table
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.registry import register_runner
from ensembl.vep.stats.summary import RunSummary, RunSummaryCollector

__all__ = ["Runner"]

logger = logging.getLogger(__name__)


@register_runner("corpus.describe")
class Runner:
    """Runner for describing a registered corpus.

    With ``--annotate``, also run the configured annotator over the whole corpus and count the
    terms it emits -- what the corpus exercises, not only what it contains.

    Attributes:
        name (ClassVar[str]): The registry name.
        sections (ClassVar[tuple[str, ...]]): The configuration sections read.
        description (ClassVar[str]): One-line description for the subcommand list.
    """

    name: ClassVar[str] = "corpus.describe"
    sections: ClassVar[tuple[str, ...]] = ("logging", "corpus", "annotation", "corpus_describe")
    description: ClassVar[str] = "describe a registered corpus and list its selectable overlaps"

    def __init__(self, config: Config) -> None:
        """Hold the configuration."""
        self.config = config

    def run(self) -> ExitCode:
        """Describe the configured corpus in the requested format."""
        corpus = self.config.corpus
        try:
            built_corpus = build_corpus_from_options(corpus)
        except CorpusOptionsError as unusable:
            logger.error("%s", unusable)
            return ExitCode.USAGE
        try:
            return self._describe(built_corpus)
        finally:
            built_corpus.close()

    def _describe(self, built_corpus: Corpus) -> ExitCode:
        """Describe one successfully built corpus."""
        describe = self.config.corpus_describe
        description = CorpusDescription.build(
            built_corpus,
            select=describe.select,
            batch_size=self.config.annotation.batch_size,
        )
        if describe.select is not None and not description.overlaps:
            logger.error("no overlap is named %r in %s", describe.select, built_corpus.description)
            return ExitCode.USAGE

        logger.info(
            "%s: %d transcripts, %d overlaps, %.1f distinct transcripts per chunk of %d",
            built_corpus.description,
            len(description.transcripts),
            len(description.overlaps),
            description.transcripts_per_chunk,
            description.batch_size,
        )
        summary = self._summary(built_corpus) if describe.annotate else None
        try:
            table_format = resolve_report_format(
                describe.output_format, graph=False, table_default="markdown", graph_formats=()
            )
        except ValueError as unusable:
            logger.error("%s", unusable)
            return ExitCode.USAGE
        self._emit(self._render(description, table_format, summary))
        return ExitCode.OK

    def _summary(self, corpus: Corpus) -> RunSummary:
        """Annotate every variant of the corpus with the configured annotator and summarise the terms."""
        annotation = self.config.annotation
        annotator = build_annotator(
            **asdict(annotation),
            upstream_distance=self.config.source.upstream_distance,
            downstream_distance=self.config.source.downstream_distance,
        )
        collector = RunSummaryCollector(metadata={"annotator": annotation.annotator})
        for variation in annotated_variations(corpus, annotator):
            collector.observe(variation)
        summary = collector.finish()
        logger.info("annotated %s with %s", summary.subtitle, annotation.annotator)
        return summary

    @staticmethod
    def _render(description: CorpusDescription, table_format: str, summary: RunSummary | None = None) -> str:
        """Render the description's tables, and the term counts when there are any, as one page or text."""
        tables: list[tuple[str, str, tuple[str, ...], list[tuple[object, ...]]]] = [
            ("chromosomes", "Chromosomes", CHROMOSOME_HEADERS, description.chromosome_rows()),
            ("transcripts", "Transcripts", TRANSCRIPT_HEADERS, description.transcript_rows()),
            ("classes", "Variant classes", CLASS_HEADERS, description.class_rows()),
        ]
        if summary is not None:
            tables.extend(
                [
                    (
                        "most_severe",
                        "Most severe consequences",
                        ("Term", "Variants"),
                        list(summary.most_severe),
                    ),
                    ("consequences", "Consequences", ("Term", "Occurrences"), list(summary.consequences)),
                    (
                        "coding",
                        "Coding consequences",
                        ("Term", "Occurrences"),
                        list(summary.coding_consequences),
                    ),
                ]
            )
        tables.append(("overlaps", "Overlaps", OVERLAP_HEADERS, description.overlap_rows()))
        chosen = TableFormat(table_format)
        subtitle = (
            f"{description.transcripts_per_chunk:.1f} distinct transcripts per chunk "
            f"of {description.batch_size} overlaps"
        )
        if chosen is TableFormat.HTML:
            sections = []
            for identifier, title, headers, rows in tables:
                sections.append(
                    Section(identifier=identifier, title=title, blocks=(render_table(headers, rows, chosen),))
                )
            return render_page(f"Corpus: {description.kind}", sections, subtitle)

        blocks = [f"# Corpus: {description.kind}", subtitle]
        for _identifier, title, headers, rows in tables:
            blocks.append(f"\n## {title}\n")
            blocks.append(render_table(headers, rows, chosen))
        return "\n".join(blocks)

    def _emit(self, rendered: str) -> None:
        """Write the rendering to the configured destination."""
        destination = self.config.corpus_describe.output_file
        if destination is None:
            print(rendered)
            return
        destination.write_text(rendered, encoding="utf-8")
        logger.info("wrote %s", destination)
