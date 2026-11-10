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
"""The annotation worker: a payload of unparsed input records in, detached annotations out.

One worker per process, built after the process starts, holding what every track shares: the
reference sequence, the synonym map, the input source turning lines into variations, and the
output factory detaching what a payload produced. Which annotations are made, and from which
sources, belongs to the workflow's tracks rather than here.

The input source is what reads: the worker hands it one payload's records and is given variations
back, so which format the input is in, and what a record of it is, stays with the source.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from itertools import batched
from typing import Any

from attrs import frozen

from ensembl.vep.buffer import FixedSizeInputBuffer
from ensembl.vep.config import Config
from ensembl.vep.models.variation import Variation, VariationAllele
from ensembl.vep.output.factory.base import AnnotatedVariant
from ensembl.vep.output.factory.base import AnnotatedVariantFactory
from ensembl.vep.sources.build import (
    build_chromosome_synonyms,
    build_input_source,
    build_sequence_source,
)
from ensembl.vep.stats import RunSummaryCollector, RunSummaryDelta
from ensembl.vep.workflow.base import Workflow
from ensembl.vep.workflow.build.registry import build_workflow
from ensembl.vep.workflow.tracks.base import TrackResources, TrackResult

__all__ = ["AnnotationBufferResult", "AnnotationBufferWorker"]

logger = logging.getLogger(__name__)


@frozen(kw_only=True, auto_attribs=True)
class AnnotationBufferResult:
    """Ordered worker result detached from live model objects.

    Attributes:
        variants (tuple[AnnotatedVariant, ...]): Detached annotated variants.
        summary_delta (RunSummaryDelta): Mergeable run-summary counts.
    """

    variants: tuple[AnnotatedVariant, ...]
    summary_delta: RunSummaryDelta


class AnnotationBufferWorker:
    """Worker-local sources and workflow state for payloads of input records.

    The ``Worker`` protocol's ``run`` for annotation, which is why nothing here is constructed
    until the process that annotates owns it: sources are opened and caches warmed per process.
    """

    def __init__(self, config: Config) -> None:
        """Build the shared sources and the configured workflow after the worker starts."""
        source_config = config.source
        assert source_config.input_file is not None
        self.config = config
        synonyms = build_chromosome_synonyms(source_config)
        self.sequence_source = build_sequence_source(source_config)
        self.input_source = build_input_source(
            source_config, sequence_source=self.sequence_source, synonyms=synonyms
        )
        self.workflow: Workflow = build_workflow(
            config,
            TrackResources(
                source_config=source_config,
                sequence_source=self.sequence_source,
                synonyms=synonyms,
            ),
        )
        self.snapshot_factory = AnnotatedVariantFactory()
        self._track_identifiers = frozenset(track.identifier for track in self.workflow.tracks)

    def run(self, records: tuple[str, ...]) -> AnnotationBufferResult:
        """Return detached annotations for one payload of input records."""
        variations = self.input_source.records_to_features(records)
        unknown_regions = frozenset(self.sequence_source.unknown_regions)
        alleles = self._annotatable_alleles(variations, unknown_regions)
        annotated_variants: list[AnnotatedVariant] = []
        collector = RunSummaryCollector()

        # Each track fetches and caches what the whole batch needs once. Windows retain only
        # their temporary overlap graph.
        self.workflow.begin_payload(self._payload_buffer(alleles))
        # Read here rather than in the parent: output is the only consumer of the columns, and
        # the parser reads the record itself.
        input_columns = (self.input_source.record_columns(record) for record in records)
        variation_columns = zip(variations, input_columns, strict=True)
        for group in batched(variation_columns, self.config.execution.annotation_window):
            group_variations = [variation for variation, _columns in group]
            group_alleles = self._annotatable_alleles(group_variations, unknown_regions)
            results: list[TrackResult[Any]] = []
            if group_alleles:
                results = self.workflow.annotate(
                    FixedSizeInputBuffer(objects=tuple(group_alleles)),
                    threads=self.config.execution.threads,
                )
            self._attach_results(group_variations, results)
            for variation, columns in group:
                region = variation.slice.seq_region_name if variation.slice is not None else None
                annotated_variants.append(
                    self.snapshot_factory.build_annotated_variant(
                        variation,
                        columns,
                        skipped=region in unknown_regions,
                    )
                )
                collector.observe(variation)
            self._attach_results(group_variations, ())
        return AnnotationBufferResult(
            variants=tuple(annotated_variants),
            summary_delta=collector.delta(),
        )

    @staticmethod
    def _payload_buffer(alleles: Sequence[VariationAllele]) -> FixedSizeInputBuffer[VariationAllele]:
        """Return the buffer priming every track for one whole batch.

        The buffer is filled through its own refill policy, which is what checks that the input is
        sorted and that one payload names one sequence region. Its cap is the payload itself:
        batching already capped it at ``buffer_size`` records, and a smaller cap here would hold
        alleles back from the prefetch the tracks prime from.

        Args:
            alleles (Sequence[VariationAllele]): Every annotatable allele the batch describes.

        Returns:
            FixedSizeInputBuffer[VariationAllele]: The filled buffer.
        """
        buffer: FixedSizeInputBuffer[VariationAllele] = FixedSizeInputBuffer(buffer_size=max(len(alleles), 1))
        buffer.next(iter(alleles))
        return buffer

    def _attach_results(
        self,
        variations: Sequence[Variation],
        results: Sequence[TrackResult[Any]],
    ) -> None:
        """Hand each track its own results to attach, and release what the others held.

        Every track is given a group, empty or not, so a window that found nothing for a track
        releases the overlaps the previous window attached for it. Attaching nothing at all is
        therefore how a rendered window's overlap graph is cleared.

        Raises:
            ValueError: If a result names a track this workflow does not have.
        """
        grouped: dict[str, list[TrackResult[Any]]] = {}
        for result in results:
            grouped.setdefault(result.track_identifier, []).append(result)
        unknown = sorted(set(grouped) - self._track_identifiers)
        if unknown:
            raise ValueError(f"workflow '{self.workflow.name}' has no track {', '.join(unknown)}")
        for track in self.workflow.tracks:
            track.attach(variations, grouped.get(track.identifier, ()))

    @staticmethod
    def _annotatable_alleles(
        variations: Sequence[Variation],
        unknown_regions: frozenset[str],
    ) -> list[VariationAllele]:
        """Return alleles whose sequence region is available to annotation sources."""
        alleles: list[VariationAllele] = []
        for variation in variations:
            region = variation.slice.seq_region_name if variation.slice is not None else None
            if region not in unknown_regions:
                alleles.extend(variation.variation_alleles)
        return alleles
