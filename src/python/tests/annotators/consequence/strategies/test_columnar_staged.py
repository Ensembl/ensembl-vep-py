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

# pylint: disable=redefined-outer-name
"""Tests for the columnar staged consequence annotator."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager

import pytest
from attrs import evolve

from ensembl.vep.annotators.base import Annotator
from ensembl.vep.annotators.consequence.facts.base import REGISTRY
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.strategies.columnar_staged import (
    ColumnarStagedAnnotator,
    topological_requires_order,
)
from ensembl.vep.caches.codon_cache import CodonEditLRUCache
from ensembl.vep.caches.zone_map import ZoneMapCache
from ensembl.vep.corpus.synthetic import SyntheticCorpus
from ensembl.vep.models.overlap import TranscriptVariationAllele


class _NotATranscriptOverlap:
    """Object that is not a transcript variation allele."""


@contextmanager
def _count_fact_cost() -> Iterator[dict[str, int]]:
    """Temporarily count every fact ``compute`` invocation's declared cost, globally.

    Mirrors ``bench/consequence_annotators.py``'s ``_count_fact_computations`` technique: every
    registered fact's ``compute`` is the single choke point every ``get_fact`` call ultimately
    reaches, regardless of which annotator module holds its own imported reference to
    ``get_fact``, so wrapping it here counts cost without needing any annotator to expose its
    internal ``FactContext``. Restores every original ``compute`` in a ``finally`` block.

    Yields:
        dict[str, int]: A live ``{"cost": ...}`` counter, updated as facts are computed.
    """
    counters = {"cost": 0}
    original_compute_functions = {name: registered_fact.compute for name, registered_fact in REGISTRY.items()}

    def make_counting_compute(
        fact_name: str, original_compute: Callable[[FactContext], object]
    ) -> Callable[[FactContext], object]:
        def counting_compute(fact_context: FactContext) -> object:
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


def _cost_per_context(annotator: Annotator, overlaps: Sequence[TranscriptVariationAllele]) -> float:
    """Return the declared-cost-weighted fact computations per context for ``annotator``.

    Args:
        annotator (Annotator): The annotator to measure.
        overlaps (Sequence[TranscriptVariationAllele]): The overlaps to run it over.

    Returns:
        float: Total declared cost of every fact computed, divided by the number of overlaps.
    """
    context_count = max(1, len(overlaps))
    with _count_fact_cost() as counters:
        annotator.annotate_feature_overlaps(overlaps)
    return counters["cost"] / context_count


def _sample_overlaps_by_variant(corpus: SyntheticCorpus, *, step: int) -> list[TranscriptVariationAllele]:
    """Sample variants while retaining every transcript overlap for each selected variant.

    Args:
        corpus (SyntheticCorpus): The corpus to sample.
        step (int): Retain every ``step``-th variant in corpus order.

    Returns:
        list[TranscriptVariationAllele]: Every overlap belonging to a selected variant.
    """
    variant_names = tuple(dict.fromkeys(name for name, _overlap in corpus.named_overlaps))
    selected_names = frozenset(variant_names[::step])
    return [overlap for name, overlap in corpus.named_overlaps if name in selected_names]


class TestTopologicalRequiresOrder:
    """Guards on the DAG sort the bit-parallel ``requires`` sweep depends on."""

    def test_topological_requires_order_succeeds_against_the_live_registry(self) -> None:
        """topological_requires_order() does not raise (no cycle) against the current REGISTRY."""
        order = topological_requires_order()
        assert order  # the registry does declare at least one `requires` edge today

    def test_every_fact_precedes_none_of_its_own_requires(self) -> None:
        """Every fact in the topological order appears strictly after everything it requires."""
        order = topological_requires_order()
        position = {fact_name: index for index, (fact_name, _requires) in enumerate(order)}
        for fact_name, requires_names in order:
            for required_name in requires_names:
                if required_name in position:
                    assert (
                        position[required_name] < position[fact_name]
                    ), f"{fact_name!r} is ordered before its own requirement {required_name!r}"

    def test_constructing_columnar_staged_annotator_does_not_raise(self) -> None:
        """Building ColumnarStagedAnnotator() (which computes and validates the order) succeeds."""
        annotator = ColumnarStagedAnnotator()
        assert annotator.name == "consequence.columnar-staged"

    def test_constructing_columnar_staged_zonemap_annotator_does_not_raise(self) -> None:
        """Building ColumnarStagedAnnotator with a zone map succeeds."""
        annotator = ColumnarStagedAnnotator(zone_map=ZoneMapCache())
        assert annotator.name == "consequence.columnar-staged"

    def test_cycle_is_detected(self) -> None:
        """A cycle introduced into REGISTRY's requires graph raises ValueError, not corrupt output."""
        original_fact = REGISTRY["overlap_feature"]
        # overlap_feature has no `requires` today; give it one that points back at a fact
        # that (transitively) requires overlap_feature, e.g. "overlap_exon" -- fabricating a cycle
        # purely for this test, restored in the `finally` block.
        REGISTRY["overlap_feature"] = evolve(original_fact, requires=("overlap_exon",))
        try:
            with pytest.raises(ValueError, match="cycle"):
                topological_requires_order()
        finally:
            REGISTRY["overlap_feature"] = original_fact


class TestColumnarStagedAnnotator:
    """Tests for staged columnar batch handling and validation."""

    def test_empty_helpers_return_empty_results(self) -> None:
        """Test empty evaluation and tracing chunks are no-ops."""
        annotator = ColumnarStagedAnnotator()

        assert annotator._evaluate_staged_chunk([]) == []  # pylint: disable=protected-access
        assert not annotator._trace_staged_chunk([])  # pylint: disable=protected-access

    def test_non_transcript_overlap_is_rejected(self) -> None:
        """Test a batch containing a non-transcript overlap raises a type error."""
        annotator = ColumnarStagedAnnotator()

        with pytest.raises(TypeError, match="only annotates TranscriptVariationAllele"):
            annotator.annotate_feature_overlaps([_NotATranscriptOverlap()])  # type: ignore[list-item]


class TestOptionalCaches:
    """Tests for the optional strategy caches."""

    def test_zone_map_is_used_without_changing_terms(
        self, boundary_overlaps: dict[str, TranscriptVariationAllele]
    ) -> None:
        """Test the enabled zone map serves facts without changing the result."""
        overlap = boundary_overlaps["ENST0001::ENST0001_cds3_del1"]
        expected = ColumnarStagedAnnotator().annotate_feature_overlaps([overlap])
        zone_map = ZoneMapCache()
        annotator = ColumnarStagedAnnotator(zone_map=zone_map)

        actual = annotator.annotate_feature_overlaps([overlap])

        assert actual == expected
        assert annotator.annotator_context.zone_map_cache is zone_map
        assert zone_map.maps_built == 1
        assert zone_map.facts_served > 0

    def test_codon_cache_is_reused_without_changing_terms(
        self, boundary_overlaps: dict[str, TranscriptVariationAllele]
    ) -> None:
        """Test the enabled codon cache serves a repeated coding edit."""
        overlap = boundary_overlaps["ENST0001::ENST0001_cds3_del1"]
        codon_cache = CodonEditLRUCache()
        annotator = ColumnarStagedAnnotator(codon_cache=codon_cache)

        consequences = annotator.annotate_feature_overlaps([overlap, overlap])

        assert consequences[1] == consequences[0]
        assert annotator.annotator_context.codon_cache is codon_cache
        assert codon_cache.miss_count == 1
        assert codon_cache.hit_count == 1
