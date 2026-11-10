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
"""Reciprocal interval-overlap association."""

from __future__ import annotations

from attrs import define

from ensembl.vep.matcher.base import IntervalLike

__all__ = ["ReciprocalOverlapMatcher"]


@define(kw_only=True, auto_attribs=True, frozen=True)
class ReciprocalOverlapMatcher:
    """Matcher requiring both configured reference-span fractions.

    Zero-width insertions never match reciprocally because their reference-span fraction is
    undefined. Use ``OverlapMatcher`` for point association or ``ExactCoordinateMatcher`` for
    identical insertion positions.

    Attributes:
        left_minimum (float): Minimum fraction of the left reference span covered.
        right_minimum (float): Minimum fraction of the right reference span covered.
    """

    left_minimum: float = 0.5
    right_minimum: float = 0.5

    def __attrs_post_init__(self) -> None:
        for fraction in (self.left_minimum, self.right_minimum):
            if not 0.0 <= fraction <= 1.0:
                raise ValueError("Reciprocal overlap fractions must be between zero and one.")

    def matches(self, left: IntervalLike, right: IntervalLike) -> bool:
        """Return whether both intervals meet their configured overlap fraction."""
        left_length = left.end - left.start + 1
        right_length = right.end - right.start + 1
        if left_length <= 0 or right_length <= 0:
            return False

        overlap_length = max(0, min(left.end, right.end) - max(left.start, right.start) + 1)
        return (
            overlap_length / left_length >= self.left_minimum
            and overlap_length / right_length >= self.right_minimum
        )
