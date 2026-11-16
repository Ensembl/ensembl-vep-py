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
"""Unit tests for ensembl.vep.annotators.registry."""

from __future__ import annotations

from collections.abc import Sequence
from typing import cast

import pytest

from ensembl.vep.annotators import registry as registry_module
from ensembl.vep.annotators.registry import available_annotators, get_annotator, register_annotator
from ensembl.vep.annotators.base import Annotator
from ensembl.vep.models.overlap import FeatureOverlap


class _StringAnnotator(Annotator[str]):
    """Minimal annotator used to test registration and lookup."""

    name = "test.string"

    def annotate_feature_overlap(self, overlap: FeatureOverlap) -> list[str]:  # noqa: D102
        """Annotate feature overlap."""
        del overlap
        return ["alpha"]

    def annotate_feature_overlaps(self, overlaps: Sequence[FeatureOverlap]) -> list[list[str]]:  # noqa: D102
        """Annotate a batch of feature overlaps."""
        return [self.annotate_feature_overlap(overlap) for overlap in overlaps]


class TestAvailableAnnotators:
    """The registry always exposes the built-in annotators."""

    def test_includes_the_built_ins(self) -> None:
        """Built-ins are registered on import."""
        assert "consequence.linear" in available_annotators()
        assert "haplo.linear" in available_annotators()


class TestLookup:
    """Lookup behavior for registered and unregistered annotators."""

    def test_unknown_name_raises_key_error_with_the_available_list(self) -> None:
        """Unknown names raise KeyError with the registry summary in the message."""
        with pytest.raises(KeyError, match="Unknown annotator 'not.a.real.annotator'"):
            get_annotator("not.a.real.annotator")

    def test_registered_string_annotator_can_be_resolved(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """The registry can store annotators that do not produce consequences."""
        monkeypatch.setitem(registry_module._ANNOTATOR_REGISTRY, "test.string", _StringAnnotator)

        annotator = cast(_StringAnnotator, get_annotator("test.string"))
        overlap = cast(FeatureOverlap, object())

        assert isinstance(annotator, _StringAnnotator)
        assert annotator.annotate_feature_overlap(overlap) == ["alpha"]
        assert annotator.annotate_feature_overlaps([overlap, overlap]) == [["alpha"], ["alpha"]]


class TestRegistration:
    """Registration works as a plain call or as a decorator."""

    def test_register_annotator_as_a_plain_call(self) -> None:
        """A factory can be registered directly."""
        register_annotator("test.string", _StringAnnotator)
        try:
            assert "test.string" in available_annotators()
            assert isinstance(get_annotator("test.string"), _StringAnnotator)
        finally:
            registry_module._ANNOTATOR_REGISTRY.pop("test.string", None)  # pylint: disable=protected-access

    def test_register_annotator_as_a_decorator(self) -> None:
        """A factory can also be registered as a decorator."""

        @register_annotator("test.decorated")  # type: ignore[arg-type]
        class _DecoratedStringAnnotator(_StringAnnotator):
            name = "test.decorated"

        try:
            assert "test.decorated" in available_annotators()
            assert isinstance(get_annotator("test.decorated"), _DecoratedStringAnnotator)
        finally:
            registry_module._ANNOTATOR_REGISTRY.pop(
                "test.decorated", None
            )  # pylint: disable=protected-access
