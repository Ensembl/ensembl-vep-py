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
"""Tests for annotator construction."""

from __future__ import annotations

from collections.abc import Sequence

import pytest

from ensembl.vep.annotators import factory as factory_module
from ensembl.vep.annotators import registry as registry_module
from ensembl.vep.annotators.consequence.context.annotator import AnnotatorContext
from ensembl.vep.annotators.consequence.facts import base as facts_base
from ensembl.vep.annotators.factory import build_annotator, build_annotator_context
from ensembl.vep.caches.codon_cache import CodonEditLRUCache
from ensembl.vep.caches.mapper_cache import MapperIndexLRUCache
from ensembl.vep.caches.zone_map import ZoneMapCache
from ensembl.vep.common import Readthrough
from ensembl.vep.models.overlap import FeatureOverlap


class _ContextOnlyAnnotator:
    """Annotator stand-in that only accepts a shared context."""

    name = "test.context_only"

    def __init__(self, annotator_context: AnnotatorContext | None = None) -> None:
        self.annotator_context = annotator_context

    def annotate_feature_overlap(self, overlap: FeatureOverlap) -> list[str]:  # noqa: D102
        """Annotate feature overlap."""
        del overlap
        return []

    def annotate_feature_overlaps(self, overlaps: Sequence[FeatureOverlap]) -> list[list[str]]:  # noqa: D102
        """Annotate a batch of feature overlaps."""
        return [self.annotate_feature_overlap(overlap) for overlap in overlaps]


class _ConfiguredAnnotator:
    """Annotator stand-in that accepts every keyword build_annotator may forward."""

    name = "test.configured"

    def __init__(
        self,
        annotator_context: AnnotatorContext | None = None,
        batch_size: int = 0,
        early_subsumption_skip: bool = False,
    ) -> None:
        self.annotator_context = annotator_context
        self.batch_size = batch_size
        self.early_subsumption_skip = early_subsumption_skip

    def annotate_feature_overlap(self, overlap: FeatureOverlap) -> list[str]:  # noqa: D102
        """Annotate feature overlap."""
        del overlap
        return ["configured"]

    def annotate_feature_overlaps(self, overlaps: Sequence[FeatureOverlap]) -> list[list[str]]:  # noqa: D102
        """Annotate a batch of feature overlaps."""
        return [self.annotate_feature_overlap(overlap) for overlap in overlaps]


class TestBuildAnnotatorContext:
    """The context builder wires caches and readthrough together."""

    def test_defaults_disable_every_cache(self) -> None:
        """The default context is empty apart from readthrough."""
        context = build_annotator_context()
        assert context.zone_map_cache is None
        assert context.codon_cache is None
        assert context.mapper_index_cache is None
        assert context.readthrough is Readthrough.NONE

    def test_enables_requested_caches_and_coerces_readthrough(self) -> None:
        """Requested caches are created and string readthrough names are accepted."""
        context = build_annotator_context(
            zone_map=True,
            zone_map_bytes=56,
            codon_cache=True,
            codon_cache_entries=12,
            mapper_cache=True,
            mapper_cache_bytes=34,
            readthrough="utr",
        )

        assert isinstance(context.zone_map_cache, ZoneMapCache)
        assert context.zone_map_cache.max_bytes == 56
        assert isinstance(context.codon_cache, CodonEditLRUCache)
        assert context.codon_cache.max_entries == 12
        assert isinstance(context.mapper_index_cache, MapperIndexLRUCache)
        assert context.mapper_index_cache.max_bytes == 34
        assert context.readthrough is Readthrough.UTR

    @pytest.mark.parametrize("enabled", [False, True])
    def test_fact_cost_override_applies_the_packaged_table_only_when_set(
        self, enabled: bool, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The packaged cost table is applied exactly when ``fact_cost_override`` is set."""
        calls: list[bool] = []
        monkeypatch.setattr(factory_module, "apply_packaged_fact_cost_overrides", lambda: calls.append(True))

        build_annotator_context(fact_cost_override=enabled)

        assert calls == ([True] if enabled else [])

    def test_fact_cost_override_changes_registered_costs_and_not_the_default(self) -> None:
        """Building without the option leaves declared costs; with it the packaged table is live."""
        declared = {name: registered.cost for name, registered in facts_base.REGISTRY.items()}
        packaged = facts_base.load_fact_cost_overrides_from_package()
        try:
            build_annotator_context()
            assert {name: registered.cost for name, registered in facts_base.REGISTRY.items()} == declared

            build_annotator_context(fact_cost_override=True)
            for name, cost in packaged.items():
                if name in facts_base.REGISTRY:
                    assert facts_base.REGISTRY[name].cost == cost
        finally:
            facts_base.apply_fact_cost_overrides(declared)


class TestBuildAnnotator:
    """The annotator builder filters kwargs by factory signature."""

    def test_builds_an_annotator_with_a_context_only_signature(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Unsupported kwargs are filtered out before instantiation."""
        monkeypatch.setitem(
            registry_module._ANNOTATOR_REGISTRY,  # pylint: disable=protected-access
            "test.context_only",
            _ContextOnlyAnnotator,
        )

        annotator = build_annotator(
            annotator="test.context_only",
            zone_map=True,
            zone_map_bytes=78,
            codon_cache=True,
            mapper_cache=True,
            readthrough="beyond_transcript",
        )

        assert isinstance(annotator, _ContextOnlyAnnotator)
        context = annotator.annotator_context
        assert context is not None
        assert isinstance(context.zone_map_cache, ZoneMapCache)
        assert context.zone_map_cache.max_bytes == 78
        assert isinstance(context.codon_cache, CodonEditLRUCache)
        assert isinstance(context.mapper_index_cache, MapperIndexLRUCache)
        assert context.readthrough is Readthrough.BEYOND_TRANSCRIPT

    def test_builds_with_an_explicit_context_and_forwards_supported_keywords(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """An explicit context is reused and supported kwargs reach the factory."""
        monkeypatch.setitem(
            registry_module._ANNOTATOR_REGISTRY,  # pylint: disable=protected-access
            "test.configured",
            _ConfiguredAnnotator,
        )  # pylint: disable=protected-access
        context = AnnotatorContext(readthrough=Readthrough.UTR)

        annotator = build_annotator(
            annotator="test.configured",
            annotator_context=context,
            batch_size=128,
            early_subsumption_skip=True,
        )

        assert isinstance(annotator, _ConfiguredAnnotator)
        assert annotator.annotator_context is context
        assert annotator.batch_size == 128
        assert annotator.early_subsumption_skip is True

    def test_build_annotator_forwards_fact_cost_override(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """The option reaches the context builder when no shared context is supplied."""
        calls: list[bool] = []
        monkeypatch.setattr(factory_module, "apply_packaged_fact_cost_overrides", lambda: calls.append(True))
        monkeypatch.setitem(
            registry_module._ANNOTATOR_REGISTRY,  # pylint: disable=protected-access
            "test.context_only",
            _ContextOnlyAnnotator,
        )

        build_annotator(annotator="test.context_only", fact_cost_override=True)

        assert calls == [True]

    def test_unknown_name_raises_key_error(self) -> None:
        """Unknown annotator names bubble up as KeyError."""
        with pytest.raises(KeyError, match="Unknown annotator 'not.a.real.annotator'"):
            build_annotator(annotator="not.a.real.annotator")
