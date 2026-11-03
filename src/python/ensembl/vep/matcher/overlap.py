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
"""General interval-overlap matching."""

from __future__ import annotations

from attrs import define

from ensembl.vep.caches.interval_tree import Interval
from ensembl.vep.common.enums import Strand
from ensembl.vep.matcher.base import InsertionMatchPolicy, IntervalLike

__all__ = ["OverlapMatcher"]


@define(kw_only=True, auto_attribs=True, frozen=True)
class OverlapMatcher:
    """Matcher for containing or partially intersecting intervals.

    A normal interval uses inclusive base coordinates. An insertion is encoded as
    ``start=end+1`` and remains a zero-width interbase location during matching.

    Attributes:
        insertion_policy (InsertionMatchPolicy): Rule for associating an insertion with a feature.
    """

    insertion_policy: InsertionMatchPolicy = InsertionMatchPolicy.CONTAINED

    def matches(self, left: IntervalLike, right: IntervalLike) -> bool:
        """Return whether two intervals are associated."""

        left_interval = Interval(start=left.start, end=left.end)
        right_interval = Interval(start=right.start, end=right.end)

        if left_interval.is_insertion and right_interval.is_insertion:
            return left_interval.start == right_interval.start
        if left_interval.is_insertion:
            return self._insertion_matches(left_interval, right_interval, right)
        if right_interval.is_insertion:
            return self._insertion_matches(right_interval, left_interval, left)
        return left_interval.start <= right_interval.end and left_interval.end >= right_interval.start

    def _insertion_matches(
        self,
        insertion: Interval,
        feature_interval: Interval,
        feature: IntervalLike,
    ) -> bool:
        """Return whether an insertion point associates with a feature."""

        if self.insertion_policy is InsertionMatchPolicy.CONTAINED:
            return feature_interval.start < insertion.start <= feature_interval.end
        if self.insertion_policy is InsertionMatchPolicy.TOUCHING:
            return feature_interval.start <= insertion.start <= feature_interval.end + 1
        return self._three_prime_flank_matches(insertion, feature_interval, feature)

    @staticmethod
    def _three_prime_flank_matches(
        insertion: Interval,
        feature_interval: Interval,
        feature: IntervalLike,
    ) -> bool:
        """Return whether an insertion's feature-relative 3-prime flank is covered."""

        strand = getattr(feature, "strand", Strand.UNKNOWN)
        if strand is Strand.FORWARD:
            flank = insertion.start
        elif strand is Strand.REVERSE:
            flank = insertion.end
        else:
            return False
        return feature_interval.start <= flank <= feature_interval.end
