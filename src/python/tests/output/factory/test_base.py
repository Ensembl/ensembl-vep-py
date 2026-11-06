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

from ensembl.vep.common.constants import UNKNOWN_VALUE
from ensembl.vep.common import CoordSystemLevel, Locus
from ensembl.vep.common.enums import Strand
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.models.variation import Variation
from ensembl.vep.output.factory import AnnotatedVariant, AnnotatedVariants, AnnotationRecord
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


class TestAnnotatedVariantBuilding:
    """Detaching one variation into a snapshot."""

    @staticmethod
    def _variation(uploaded_name: str | None = "1") -> Variation:
        """Build one minimal SNV, named as its input spelled it and with no slice of its own."""
        return Variation(
            start=10,
            end=10,
            strand=Strand.FORWARD,
            stable_id="rs1",
            allele_string="A/G",
            uploaded_locus=Locus(seq_region_name=uploaded_name, start=10, end=10),
        )

    def test_it_carries_the_records_it_was_given(self) -> None:
        """Which records a run writes is the output factory's decision, not the snapshot's."""
        annotations: tuple[AnnotationRecord, ...] = ({"Consequence": ["missense_variant"]},)

        snapshot = AnnotatedVariant.build(
            self._variation(), ("1", "10", "rs1", "A", "G"), annotations=annotations
        )

        assert snapshot.annotations == annotations
        assert snapshot.input_columns == ("1", "10", "rs1", "A", "G")
        assert snapshot.seq_region_name == "1"
        assert snapshot.variant_class == "SNV"

    def test_the_region_is_named_as_the_input_named_it(self) -> None:
        """Report the uploaded spelling, not the one the annotation sources were queried with."""
        variation = self._variation(uploaded_name="chr1")
        variation.slice = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            start=1,
            end=10**6,
            strand=Strand.FORWARD,
        )

        snapshot = AnnotatedVariant.build(variation)

        assert snapshot.seq_region_name == "chr1"

    def test_a_variation_that_names_no_region_is_unknown(self) -> None:
        """The uploaded locus is the only answer: a slice names the region for the sources."""
        variation = self._variation(uploaded_name=None)
        variation.slice = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name="1",
            start=1,
            end=10**6,
            strand=Strand.FORWARD,
        )

        assert AnnotatedVariant.build(variation).seq_region_name == UNKNOWN_VALUE

    def test_a_variation_with_no_uploaded_locus_is_unknown(self) -> None:
        """Nothing names the region, so output says so rather than guessing."""
        variation = Variation(start=10, end=10, strand=Strand.FORWARD, allele_string="A/G")

        assert AnnotatedVariant.build(variation).seq_region_name == UNKNOWN_VALUE

    def test_the_most_severe_term_falls_back_to_the_first_record(self) -> None:
        """A term the model reports only inside its records still reaches the snapshot."""
        annotations: tuple[AnnotationRecord, ...] = ({"Consequence": ["splice_donor_variant"]},)

        snapshot = AnnotatedVariant.build(self._variation(), annotations=annotations)

        assert snapshot.most_severe_consequence == "splice_donor_variant"

    def test_a_skipped_variant_carries_no_annotation(self) -> None:
        """A region no source could annotate is reported as skipped rather than as intergenic."""
        snapshot = AnnotatedVariant.build(self._variation(), skipped=True)

        assert snapshot.skipped
        assert not snapshot.annotations
        assert snapshot.most_severe_consequence is None


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
