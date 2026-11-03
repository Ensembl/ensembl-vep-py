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
"""Tests for the consequence annotator strategy base class."""

from __future__ import annotations

from collections.abc import Sequence

import pytest

from ensembl.vep.annotators.consequence.context import AnnotatorContext
from ensembl.vep.annotators.consequence.strategies.base import ConsequenceAnnotator
from ensembl.vep.caches.codon_cache import CodonEditLRUCache
from ensembl.vep.caches.zone_map import ZoneMapCache
from ensembl.vep.common import Readthrough
from ensembl.vep.models.overlap import FeatureOverlap, OverlapConsequence, TranscriptVariationAllele


class _TestAnnotator(ConsequenceAnnotator[str]):
    """Concrete annotator used to exercise the shared strategy contract."""

    name = "test"

    def __init__(self, *args: object, **kwargs: object) -> None:
        """Initialize the annotator and its received-overlap record."""
        super().__init__("rules", *args, **kwargs)  # type: ignore[arg-type]
        self.received_overlaps: list[FeatureOverlap] = []

    def annotate_feature_overlaps(self, overlaps: Sequence[FeatureOverlap]) -> list[list[OverlapConsequence]]:
        """Record a batch and return one empty result per overlap."""
        self.received_overlaps.extend(overlaps)
        return [[] for _overlap in overlaps]


class _AbstractMethodAnnotator(ConsequenceAnnotator[str]):
    """Concrete wrapper that calls the abstract method implementation."""

    # pylint: disable-next=useless-parent-delegation
    def annotate_feature_overlaps(self, overlaps: Sequence[FeatureOverlap]) -> list[list[OverlapConsequence]]:
        """Call the base implementation."""
        return super().annotate_feature_overlaps(overlaps)  # type: ignore[safe-super]


class TestConsequenceAnnotator:
    """Tests for shared consequence annotator behaviour."""

    def test_abstract_base_cannot_be_instantiated(self) -> None:
        """Test an implementation must provide the batch annotation method."""
        with pytest.raises(TypeError):
            ConsequenceAnnotator("rules")  # type: ignore[abstract] # pylint: disable=abstract-class-instantiated

    def test_abstract_batch_implementation_raises_not_implemented(self) -> None:
        """Test the abstract batch method has no fallback implementation."""
        annotator = _AbstractMethodAnnotator("rules")

        with pytest.raises(NotImplementedError):
            annotator.annotate_feature_overlaps([])

    def test_default_context_is_private_and_untraced(self) -> None:
        """Test construction creates a fresh context and no trace."""
        first_annotator = _TestAnnotator()
        second_annotator = _TestAnnotator()

        assert first_annotator.annotator_context is not second_annotator.annotator_context
        assert first_annotator.ruleset == "rules"
        assert first_annotator.name == "test"
        assert first_annotator.last_trace is None

    def test_supplied_options_override_context_values(self) -> None:
        """Test explicit cache and readthrough options update the supplied context."""
        original_zone_map = ZoneMapCache()
        original_codon_cache = CodonEditLRUCache()
        context = AnnotatorContext(
            zone_map_cache=original_zone_map,
            codon_cache=original_codon_cache,
        )
        zone_map = ZoneMapCache()
        codon_cache = CodonEditLRUCache()

        annotator = _TestAnnotator(
            context,
            zone_map=zone_map,
            codon_cache=codon_cache,
            readthrough=Readthrough.UTR,
        )

        assert annotator.annotator_context is context
        assert context.zone_map_cache is zone_map
        assert context.codon_cache is codon_cache
        assert context.readthrough is Readthrough.UTR

    def test_cache_counters_are_read_from_the_context(self) -> None:
        """Test cache counter properties expose the shared context counters."""
        context = AnnotatorContext()
        context.cache_hits = 3
        context.cache_misses = 2
        annotator = _TestAnnotator(context)

        assert annotator.cache_hits == 3
        assert annotator.cache_misses == 2

    def test_single_overlap_delegates_to_a_one_item_batch(
        self, boundary_overlaps: dict[str, TranscriptVariationAllele]
    ) -> None:
        """Test the single-overlap adapter delegates to the batch implementation."""
        overlap = boundary_overlaps["ENST0001::ENST0001_cds3_del1"]
        annotator = _TestAnnotator()

        consequences = annotator.annotate_feature_overlap(overlap)

        assert consequences == []
        assert annotator.received_overlaps == [overlap]
