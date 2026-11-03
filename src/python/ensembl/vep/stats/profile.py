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
"""Fact-level profiling collected by instrumenting an annotation run.

Profiling is opt-in because counting replaces each fact's ``compute`` callback and tracing records
every overlap derivation. The process-global callback replacement is restored before returning and
is not thread-safe; computed values remain local to each ``FactContext``.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager

from attrs import evolve

from ensembl.vep.annotators.base import Annotator
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.facts.base import REGISTRY
from ensembl.vep.models.overlap import OverlapConsequence, TranscriptVariationAllele
from ensembl.vep.trace import annotation_tracing
from ensembl.vep.trace.annotation import AnnotationTrace
from ensembl.vep.trace.run import RunTrace

__all__ = ["count_fact_computations", "trace_run"]

logger = logging.getLogger(__name__)


@contextmanager
def count_fact_computations() -> Iterator[dict[str, int]]:
    """Count every fact ``compute`` call and its declared cost, process-wide, for the duration.

    Yields:
        dict[str, int]: A live ``{"computed": ..., "cost": ...}`` counter, updated as facts are
            computed while the block is active.
    """
    counters = {"computed": 0, "cost": 0}
    original_compute_functions = {name: registered_fact.compute for name, registered_fact in REGISTRY.items()}

    def make_counting_compute(
        fact_name: str, original_compute: Callable[[FactContext], object]
    ) -> Callable[[FactContext], object]:
        def counting_compute(fact_context: FactContext) -> object:
            counters["computed"] += 1
            counters["cost"] += REGISTRY[fact_name].cost
            return original_compute(fact_context)

        return counting_compute

    try:
        for fact_name, registered_fact in list(REGISTRY.items()):
            REGISTRY[fact_name] = evolve(
                registered_fact, compute=make_counting_compute(fact_name, registered_fact.compute)
            )
        yield counters
    finally:
        for fact_name, original_compute in original_compute_functions.items():
            REGISTRY[fact_name] = evolve(REGISTRY[fact_name], compute=original_compute)


def trace_run(
    annotator: Annotator[OverlapConsequence],
    overlaps: Sequence[TranscriptVariationAllele],
    *,
    annotator_name: str = "",
    cost_table: str = "",
) -> RunTrace:
    """Annotate every overlap with tracing on and fold the derivations into one run-scoped view.

    Args:
        annotator (Annotator[OverlapConsequence]): The annotator to run. A native one never
            traces, so its run folds to nothing.
        overlaps (Sequence[TranscriptVariationAllele]): The overlaps to annotate, in order.
        annotator_name (str): The annotator's registry name, for the report.
        cost_table (str): Which cost table the declared costs came from, for the report.

    Returns:
        RunTrace: The fold.
    """
    traces: list[AnnotationTrace] = []
    with annotation_tracing(quiet=True):
        for overlap in overlaps:
            annotator.annotate_feature_overlap(overlap)
            trace = getattr(annotator, "last_trace", None)
            if trace is not None:
                traces.append(trace)
    if not traces and overlaps:
        logger.warning("%s recorded no traces; a native annotator never does", annotator_name or "annotator")
    return RunTrace.fold(traces, annotator_name=annotator_name, cost_table=cost_table)
