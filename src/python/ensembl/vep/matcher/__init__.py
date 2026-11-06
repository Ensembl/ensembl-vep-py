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
"""Association matchers."""

from __future__ import annotations

from .allele import get_matched_variant_alleles
from .base import AlleleLike, InsertionMatchPolicy, IntervalLike
from .exact import ExactCoordinateMatcher
from .overlap import OverlapMatcher
from .reciprocal import ReciprocalOverlapMatcher

__all__ = [
    "AlleleLike",
    "ExactCoordinateMatcher",
    "InsertionMatchPolicy",
    "IntervalLike",
    "OverlapMatcher",
    "ReciprocalOverlapMatcher",
    "get_matched_variant_alleles",
]
