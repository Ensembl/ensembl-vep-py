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
"""The default annotation runner."""

from __future__ import annotations

import logging
import sys
from contextlib import nullcontext
from pathlib import Path
from typing import ClassVar

from ensembl.vep.buffer.reader import input_header_lines
from ensembl.vep.config import Config
from ensembl.vep.output.factory import AnnotatedVariants
from ensembl.vep.output.factory.write import available_output_writers, get_output_writer
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.run.output_transaction import OutputTransaction
from ensembl.vep.runner.registry import register_runner
from ensembl.vep.runner.run.stream import annotate_buffers
from ensembl.vep.stats import RunSummaryCollector, render_reportable

__all__ = ["Runner"]

logger = logging.getLogger(__name__)

_SUMMARY_SUFFIXES = {"html": ".html", "text": ".txt", "json": ".json"}


@register_runner("run.default", default_for_group=True)
class Runner:
    """Annotate an input file.

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
        """Return ordered annotations for Python API consumers.

        Args:
            show_progress (bool): Show interactive progress when configuration and stderr allow it.
        """
        return annotate_buffers(self.config, show_progress=show_progress)

    def _validate_run(
        self,
        destination: Path | None,
        summary_destination: Path | None,
    ) -> None:
        """Validate source, format, and every destination before doing work."""
        source_config = self.config.source
        run_config = self.config.run
        if source_config.input_file is None:
            raise ValueError("give --input_file to annotate; 'vep corpus describe' explores a corpus instead")
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
                for header in writer.headers(input_header_lines(Path(input_file))):
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
