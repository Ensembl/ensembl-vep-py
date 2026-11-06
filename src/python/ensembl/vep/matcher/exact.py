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
"""Exact-coordinate association."""

from __future__ import annotations

from attrs import define

from ensembl.vep.matcher.base import IntervalLike

__all__ = ["ExactCoordinateMatcher"]


@define(kw_only=True, auto_attribs=True, frozen=True)
class ExactCoordinateMatcher:
    """Matcher for identical inclusive interval coordinates."""

    def matches(self, left: IntervalLike, right: IntervalLike) -> bool:
        """Return whether both coordinate pairs are identical."""
        return left.start == right.start and left.end == right.end
