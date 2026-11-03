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
"""Shared fixtures for buffer tests."""

from __future__ import annotations

import pytest

from ensembl.vep.common import CoordSystemLevel, Strand
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.models.variation import Variation


def _genomic_slice(seq_region_name: str) -> Slice:
    """Return a whole-chromosome slice for one sequence region."""
    return Slice(
        coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
        seq_region_name=seq_region_name,
        start=1,
        end=250_000_000,
        strand=Strand.FORWARD,
    )


def _variations(count: int, *, start: int = 1, seq_region_name: str = "1") -> list[Variation]:
    """Return consecutive single-base variations on one sequence region."""
    return [
        Variation(
            start=position,
            end=position,
            strand=Strand.FORWARD,
            slice=_genomic_slice(seq_region_name),
            allele_string="A/T",
        )
        for position in range(start, start + count)
    ]


@pytest.fixture
def variation_list() -> list[Variation]:
    """Return eight variations on region 1 followed by four on region 2."""
    return _variations(8) + _variations(4, start=10, seq_region_name="2")


@pytest.fixture
def variation_list_out_of_order() -> list[Variation]:
    """Return variations on region 1 whose starts jump backwards three times."""
    return [
        Variation(
            start=position,
            end=position,
            strand=Strand.FORWARD,
            slice=_genomic_slice("1"),
            allele_string="A/T",
        )
        for position in (100, 20, 120, 60, 5)
    ]
