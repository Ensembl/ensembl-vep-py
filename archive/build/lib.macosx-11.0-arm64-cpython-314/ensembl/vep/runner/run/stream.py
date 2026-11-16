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
"""Stream ordered annotated variants from buffered workers."""

from __future__ import annotations

import sys
from collections.abc import Iterator
from functools import partial
from pathlib import Path

from ensembl.vep.buffer.reader import count_vcf_records, iter_vcf_buffers
from ensembl.vep.config import Config
from ensembl.vep.output.factory.base import AnnotatedVariant, AnnotatedVariants
from ensembl.vep.runner.run.progress import VariantProgress
from ensembl.vep.runner.run.worker import AnnotatedBufferResult, annotate_buffer, build_buffered_worker
from ensembl.vep.stats import RunSummaryCollector
from ensembl.vep.workflow.execution import ExecutionSpec, ordered_map

__all__ = ["annotate_buffers"]


def annotate_buffers(config: Config, *, show_progress: bool = False) -> AnnotatedVariants:
    """Return a closable stream of annotated variants in input order."""
    input_file = config.source.input_file
    if input_file is None:
        raise ValueError("give --input_file to annotate")

    input_path = Path(input_file)
    specification = ExecutionSpec.from_config(config.execution)
    payloads = iter_vcf_buffers(
        input_path,
        batch_size=specification.buffer_size,
    )
    results = ordered_map(
        payloads,
        specification=specification,
        worker_factory=partial(build_buffered_worker, config),
        execute=annotate_buffer,
    )
    collector = RunSummaryCollector(
        metadata={
            "input": str(input_file),
            "annotator": config.annotation.annotator,
        }
    )
    variants = _iter_annotated_variants(
        results,
        collector=collector,
        input_path=input_path,
        show_progress=(
            show_progress and not config.run.no_progress and not config.logging.quiet and sys.stderr.isatty()
        ),
    )
    return AnnotatedVariants(variants=variants, collector=collector)


def _iter_annotated_variants(
    results: Iterator[AnnotatedBufferResult],
    *,
    collector: RunSummaryCollector,
    input_path: Path,
    show_progress: bool,
) -> Iterator[AnnotatedVariant]:
    """Merge worker state while yielding each detached variant."""
    variant_total = count_vcf_records(input_path) if show_progress else 0
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
