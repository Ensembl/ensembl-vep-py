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
"""Contracts shared by overlap finders."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from typing import Protocol, TypeVar, runtime_checkable

from ensembl.vep.models import Feature
from ensembl.vep.models.overlap import FeatureOverlap
from ensembl.vep.models.variation import VariationAllele
from ensembl.vep.caches.interval_index import Interval
from ensembl.vep.matcher.base import IntervalLike

__all__ = [
    "CandidateLookup",
    "CandidateIntervals",
    "FeatureMatcher",
    "MatchSubject",
    "OverlapFactory",
    "OverlapFinder",
    "VariationAlleleBuffer",
]

FeatureT = TypeVar("FeatureT", bound=Feature)
"""The feature type returned by a candidate lookup."""

CandidateFeatureT_co = TypeVar("CandidateFeatureT_co", bound=Feature, covariant=True)
"""The covariant feature type returned by a candidate lookup."""

OverlapT = TypeVar("OverlapT", bound=FeatureOverlap)
"""The overlap model returned by a finder."""

FoundOverlapT_co = TypeVar("FoundOverlapT_co", bound=FeatureOverlap, covariant=True)
"""The covariant overlap model yielded by a finder."""


@runtime_checkable
class VariationAlleleBuffer(Protocol):
    """Expose variation alleles held by an input buffer."""

    @property
    def objects(self) -> Iterable[VariationAllele]:
        """Return the variation alleles held by the buffer."""
        raise NotImplementedError

    def get_objects_by_intervals(self, intervals: Iterable[Interval]) -> list[VariationAllele]:
        """Return alleles that are candidates for any query interval."""
        raise NotImplementedError


class CandidateLookup(Protocol[CandidateFeatureT_co]):
    """Return source candidates for an input buffer."""

    def __call__(self, input_buffer: VariationAlleleBuffer) -> Iterable[CandidateFeatureT_co]:
        """Return candidate features for ``input_buffer``."""
        raise NotImplementedError


CandidateIntervals = Callable[[FeatureT], Iterable[Interval]]
"""Build buffer-index query intervals for one source feature."""


class FeatureMatcher(Protocol):
    """Determine whether a feature and variation allele are associated."""

    def matches(self, feature: IntervalLike, variation_allele: IntervalLike) -> bool:
        """Return whether ``feature`` matches ``variation_allele``."""
        raise NotImplementedError


MatchSubject = Callable[[FeatureT], IntervalLike]
"""Build the interval used for final association of one source feature."""


OverlapFactory = Callable[[FeatureT, VariationAllele], OverlapT]
"""Construct a typed overlap from a matched feature and variation allele."""


@runtime_checkable
class OverlapFinder(Protocol[FoundOverlapT_co]):
    """Find typed overlaps for one input buffer."""

    def find(self, input_buffer: VariationAlleleBuffer) -> Iterator[FoundOverlapT_co]:
        """Yield overlaps for ``input_buffer``."""
        raise NotImplementedError
