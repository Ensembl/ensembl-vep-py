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
"""Tests for the annotator protocol."""

from __future__ import annotations

from collections.abc import Sequence
from typing import cast

import pytest

from ensembl.vep.annotators.base import Annotator
from ensembl.vep.models.overlap import FeatureOverlap


class _StructuralAnnotator:
    """A concrete object that satisfies the annotator protocol structurally."""

    name = "test.structural"

    def annotate_feature_overlap(self, overlap: FeatureOverlap) -> list[str]:  # noqa: D102
        """Annotate feature overlap."""
        del overlap
        return ["term"]

    def annotate_feature_overlaps(self, overlaps: Sequence[FeatureOverlap]) -> list[list[str]]:  # noqa: D102
        """Annotate a batch of feature overlaps."""
        return [self.annotate_feature_overlap(overlap) for overlap in overlaps]


class _ProtocolAnnotator(Annotator[str]):
    """Concrete subclass that inherits the protocol defaults."""

    name = "test.protocol"

    def annotate_feature_overlap(self, overlap: FeatureOverlap) -> list[str]:  # noqa: D102
        """Defer to the protocol default."""
        return super().annotate_feature_overlap(overlap)  # type: ignore[safe-super]

    def annotate_feature_overlaps(self, overlaps: Sequence[FeatureOverlap]) -> list[list[str]]:  # noqa: D102
        """Defer to the protocol default."""
        return super().annotate_feature_overlaps(overlaps)  # type: ignore[safe-super]


class TestAnnotatorProtocol:
    """The base annotator contract behaves like a runtime-checkable protocol."""

    def test_matching_objects_satisfy_the_protocol(self) -> None:
        """A structurally compatible object passes ``isinstance``."""
        annotator = _StructuralAnnotator()
        assert isinstance(annotator, Annotator)
        overlap = cast(FeatureOverlap, object())
        assert annotator.annotate_feature_overlap(overlap) == ["term"]
        assert annotator.annotate_feature_overlaps([overlap, overlap]) == [["term"], ["term"]]

    def test_the_default_methods_raise_not_implemented_error(self) -> None:
        """A concrete subclass that inherits the protocol methods raises when called."""
        annotator = _ProtocolAnnotator()
        overlap = cast(FeatureOverlap, object())

        with pytest.raises(NotImplementedError):
            annotator.annotate_feature_overlap(overlap)

        with pytest.raises(NotImplementedError):
            annotator.annotate_feature_overlaps([overlap])
