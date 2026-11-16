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
"""Tests for format-neutral annotated variants."""

from __future__ import annotations

from collections.abc import Iterator

import pytest

from ensembl.vep.output.factory import AnnotatedVariant, AnnotatedVariants
from ensembl.vep.stats import RunSummaryCollector, RunSummaryDelta


def _annotated_variant(identifier: str) -> AnnotatedVariant:
    """Build a minimal detached variant."""
    return AnnotatedVariant(
        id=identifier,
        seq_region_name="1",
        start=10,
        end=10,
        strand=1,
        allele_string="A/G",
    )


class TestAnnotatedVariants:
    """Annotation stream lifecycle and summary access."""

    def test_completed_stream_exposes_its_count_and_summary(self) -> None:
        """Publish summary state only after every variant is consumed."""
        collector = RunSummaryCollector()
        collector.merge(RunSummaryDelta(variants=2))
        annotations = AnnotatedVariants(
            variants=(_annotated_variant("first"), _annotated_variant("second")),
            collector=collector,
        )

        assert list(annotations)
        assert annotations.completed
        assert annotations.variant_count == 2
        assert annotations.collector is collector
        assert annotations.summary.variants == 2

    def test_incomplete_stream_has_no_final_summary(self) -> None:
        """Reject a final-summary request while annotations remain unread."""
        annotations = AnnotatedVariants(
            variants=(_annotated_variant("first"),),
            collector=RunSummaryCollector(),
        )

        with pytest.raises(RuntimeError, match="after the annotation stream is exhausted"):
            _summary = annotations.summary
        with pytest.raises(RuntimeError, match="after the annotation stream is exhausted"):
            _collector = annotations.collector

    def test_context_manager_closes_a_partially_consumed_iterator(self) -> None:
        """Release the underlying iterator when an API consumer stops early."""
        closed = False

        def variants() -> Iterator[AnnotatedVariant]:
            nonlocal closed
            try:
                yield _annotated_variant("first")
                yield _annotated_variant("second")
            finally:
                closed = True

        annotations = AnnotatedVariants(
            variants=variants(),
            collector=RunSummaryCollector(),
        )

        with annotations:
            assert next(annotations).id == "first"

        assert closed
        assert not annotations.completed
        assert not list(annotations)
