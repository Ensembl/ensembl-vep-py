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
"""Variant-allele association helpers."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TypeVar

from ensembl.vep.matcher.base import AlleleLike
from ensembl.vep.matcher.exact import ExactCoordinateMatcher

__all__ = ["get_matched_variant_alleles"]

AlleleT = TypeVar("AlleleT", bound="AlleleLike")
"""The concrete allele type retained by the matching helper."""


def get_matched_variant_alleles(
    variation_alleles: Iterable[AlleleT],
    query: AlleleLike,
    *,
    match_reference: bool = True,
    match_alternate: bool = True,
) -> list[AlleleT]:
    """Return alleles with matching coordinates and selected allele strings.

    Args:
        variation_alleles (Iterable[AlleleT]): Candidate alleles to inspect.
        query (AlleleLike): Coordinate and allele values to match.
        match_reference (bool): Require equal reference strings.
        match_alternate (bool): Require equal alternate strings.

    Returns:
        list[AlleleT]: Matching candidates in their supplied order.
    """
    coordinate_matcher = ExactCoordinateMatcher()
    matched: list[AlleleT] = []
    for variation_allele in variation_alleles:
        if not coordinate_matcher.matches(variation_allele, query):
            continue
        if match_reference and variation_allele.reference_allele != query.reference_allele:
            continue
        if match_alternate and variation_allele.alternate_allele != query.alternate_allele:
            continue
        matched.append(variation_allele)
    return matched
