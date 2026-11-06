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
"""The default annotation runner.

Annotating wires three things together and then gets out of the way: the input source reads
records, a ``PayloadIterator`` groups them, and an ``Executor`` runs each payload against an
``AnnotationBufferWorker``. Results come back in payload order whatever the worker count, and
only detached annotations cross back, which is what keeps a run's live model objects inside the
process that built them.
"""

from __future__ import annotations

import logging
import sys
from collections.abc import Iterator
from contextlib import nullcontext
from functools import cached_property, partial
from pathlib import Path
from typing import ClassVar

from ensembl.vep.config import Config
from ensembl.vep.output.factory import AnnotatedVariant, AnnotatedVariants
from ensembl.vep.output.factory.write import available_output_writers, get_output_writer
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.registry import register_runner
from ensembl.vep.runner.run.output_transaction import OutputTransaction
from ensembl.vep.runner.run.progress import VariantProgress
from ensembl.vep.sources.build import (
    InputSource,
    build_chromosome_synonyms,
    build_input_source,
    detect_input_format,
)
from ensembl.vep.stats import RunSummaryCollector, render_reportable
from ensembl.vep.workflow.build import workflow_for_command
from ensembl.vep.workflow.executor import ExecutionSpec, Executor, PayloadIterator
from ensembl.vep.workflow.worker import AnnotationBufferResult, AnnotationBufferWorker

__all__ = ["Runner"]

logger = logging.getLogger(__name__)

_SUMMARY_SUFFIXES = {"html": ".html", "text": ".txt", "json": ".json"}


@register_runner("run.default", default_for_group=True)
class Runner:
    """Runner for annotating an input file.

    Attributes:
        name (ClassVar[str]): The registry name.
        sections (ClassVar[tuple[str, ...]]): The configuration sections read.
        description (ClassVar[str]): One-line description for the subcommand list.
    """

    name: ClassVar[str] = "run.default"
    sections: ClassVar[tuple[str, ...]] = (
        "logging",
        "source",
        "annotation",
        "execution",
        "run",
    )
    description: ClassVar[str] = "annotate variants against a set of annotation sources"

    def __init__(self, config: Config) -> None:
        """Initialize the runner."""
        self.config = config

    @cached_property
    def input_source(self) -> InputSource:
        """Return the source this process reads input records with, built once.

        The parent only reads: it streams records, counts them, and passes the input's own header
        through. Every model is built in the worker that annotates it, so this source is given no
        reference sequence -- reading records needs none.

        Built on first use rather than in ``__init__`` because a runner is constructed before
        anything has validated what it was given: a configuration naming no input file has to be
        reported by ``run()``, which turns it into a usage error, rather than raise while the
        runner is being built.

        Returns:
            InputSource: The input source the configured format reads with.

        Raises:
            ValueError: If no input file is configured.
        """
        source_config = self.config.source
        if source_config.input_file is None:
            raise ValueError("give --input_file to annotate")
        return build_input_source(
            source_config,
            sequence_source=None,
            synonyms=build_chromosome_synonyms(source_config),
        )

    def run(self) -> ExitCode:
        """Annotate the input and write its output and optional summary."""
        destination = self.config.run.output_file
        summary_destination = self._summary_destination()

        try:
            self._validate_run(destination, summary_destination)
            if destination is None:
                self._run_to_stdout(summary_destination)
            else:
                self._run_to_file(destination, summary_destination)
        except (OSError, ValueError) as error:
            logger.error("%s", error)
            return ExitCode.USAGE

        return ExitCode.OK

    def annotate(self, *, show_progress: bool = False) -> AnnotatedVariants:
        """Return a closable stream of annotated variants in input order.

        Args:
            show_progress (bool): Show interactive progress when configuration and stderr allow it.

        Returns:
            AnnotatedVariants: The ordered stream and the summary collector it fills.

        Raises:
            ValueError: If no input file is configured.
        """
        input_file = self.config.source.input_file
        if input_file is None:
            raise ValueError("give --input_file to annotate")

        specification = ExecutionSpec.from_config(self.config.execution)
        executor: Executor[tuple[str, ...], AnnotationBufferResult] = Executor(
            specification=specification,
            worker_factory=partial(AnnotationBufferWorker, self.config),
            payload_iterator=PayloadIterator(
                records=self.input_source.iter_all_records(),
                batch_size=specification.buffer_size,
                region_of=self.input_source.region_of,
            ),
        )
        collector = RunSummaryCollector(
            metadata={
                "input": str(input_file),
                "annotator": self.config.annotation.annotator,
            }
        )
        variants = self._iter_annotated_variants(
            executor.execute(),
            collector=collector,
            show_progress=self._progress_wanted(show_progress),
        )
        return AnnotatedVariants(variants=variants, collector=collector)

    def _progress_wanted(self, show_progress: bool) -> bool:
        """Return whether interactive progress was asked for and is worth rendering."""
        return (
            show_progress
            and not self.config.run.no_progress
            and not self.config.logging.quiet
            and sys.stderr.isatty()
        )

    def _validate_run(
        self,
        destination: Path | None,
        summary_destination: Path | None,
    ) -> None:
        """Validate source, format, and every destination before doing work."""
        source_config = self.config.source
        run_config = self.config.run
        # Before any worker starts, so an unimplemented command is reported once here rather than
        # by every worker process once it has opened its sources.
        workflow_for_command(self.config.command)
        if source_config.input_file is None:
            raise ValueError("give --input_file to annotate; 'vep corpus describe' explores a corpus instead")
        # Resolved here as well, so an input format with no reader is reported once before any
        # worker starts rather than by every worker that opens the input.
        detect_input_format(source_config.input_file, source_config.format)
        available_formats = available_output_writers()
        if run_config.output_format not in available_formats:
            supported_formats = ", ".join(available_formats)
            raise ValueError(
                f"no writer for {run_config.output_format!r} yet; supported: {supported_formats}"
            )
        if destination is not None and destination.is_dir():
            # Reached by --output_file . , which is a mistake rather than a request to stream.
            raise ValueError(f"{destination} is a directory; give a file to write, or STDOUT to stream")
        if destination is not None and summary_destination == destination:
            raise ValueError("stats_file must differ from output_file")
        if run_config.force_overwrite:
            return
        for planned_destination in (destination, summary_destination):
            if planned_destination is not None and planned_destination.exists():
                raise ValueError(f"{planned_destination} exists; pass --force_overwrite to replace it")

    def _run_to_stdout(self, summary_destination: Path | None) -> None:
        """Stream annotations to stdout and write or log their summary."""
        _variant_count, collector = self._run_buffered(None)
        if not self.config.run.no_stats:
            self._write_collector_summary(collector, summary_destination)
            if summary_destination is not None:
                logger.info("wrote %s", summary_destination)

    def _run_to_file(self, destination: Path, summary_destination: Path | None) -> None:
        """Commit annotation and summary files after both render successfully."""
        transaction = OutputTransaction(replace_existing=self.config.run.force_overwrite)
        try:
            temporary_output = transaction.stage(destination)
            variant_count, collector = self._run_buffered(temporary_output)
            if summary_destination is not None:
                temporary_summary = transaction.stage(summary_destination)
                self._write_collector_summary(collector, temporary_summary)
            transaction.commit()
        finally:
            transaction.cleanup()

        logger.info("wrote %s (%d variants)", destination, variant_count)
        if summary_destination is not None:
            logger.info("wrote %s", summary_destination)

    def _run_buffered(self, destination: Path | None) -> tuple[int, RunSummaryCollector]:
        """Consume the Python annotation stream through the selected writer."""
        input_file = self.config.source.input_file
        assert input_file is not None
        writer = get_output_writer(
            self.config.run.output_format,
            command_line=f"vep run --input_file {input_file}",
        )
        output = nullcontext(sys.stdout) if destination is None else destination.open("w", encoding="utf-8")
        annotations = self.annotate(show_progress=True)
        with output as handle:
            with annotations:
                for header in writer.headers(self.input_source.header_lines()):
                    handle.write(f"{header}\n")
                for variant in annotations:
                    handle.write(f"{writer.line(variant)}\n")
        return annotations.variant_count, annotations.collector

    def _write_collector_summary(
        self,
        collector: RunSummaryCollector,
        destination: Path | None,
    ) -> None:
        """Write a streamed summary assembled from worker deltas."""
        summary = collector.finish()
        if destination is None:
            logger.info("run summary: %s", summary.subtitle)
            return
        rendered = render_reportable(summary, self.config.run.stats_format)
        destination.write_text(rendered, encoding="utf-8")

    def _iter_annotated_variants(
        self,
        results: Iterator[AnnotationBufferResult],
        *,
        collector: RunSummaryCollector,
        show_progress: bool,
    ) -> Iterator[AnnotatedVariant]:
        """Merge worker summaries while yielding each detached variant.

        Args:
            results (Iterator[AnnotationBufferResult]): What the executor yields, in input order.
            collector (RunSummaryCollector): The run summary each result's delta is merged into.
            show_progress (bool): Whether to draw the interactive progress line.

        Yields:
            Iterator[AnnotatedVariant]: Each annotated variant, in input order.
        """
        # Never `scan=True`: a count the input does not declare would cost a whole extra pass over
        # it, which progress reporting is not worth. Without one the bar is drawn unsized.
        variant_total = self.input_source.record_count() if show_progress else None
        try:
            with VariantProgress(total=variant_total, enabled=show_progress) as progress:
                for result in results:
                    collector.merge_ordered_delta(result.summary_delta)
                    for variant in result.variants:
                        yield variant
                        progress.advance(1)
        finally:
            close = getattr(results, "close", None)
            if close is not None:
                close()

    def _summary_destination(self) -> Path | None:
        """Return the explicitly named or output-derived summary destination."""
        if self.config.run.no_stats:
            return None
        return self.config.run.stats_file or self._derived_summary_path()

    def _derived_summary_path(self) -> Path | None:
        """Return the summary path derived from output, or ``None`` for stdout."""
        output_file = self.config.run.output_file
        if output_file is None:
            return None
        suffix = _SUMMARY_SUFFIXES.get(self.config.run.stats_format, ".txt")
        return output_file.with_name(f"{output_file.name}_summary{suffix}")
