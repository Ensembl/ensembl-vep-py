# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Shared fixtures for buffer unit tests."""

from __future__ import annotations
from typing import Iterator
from attrs import define
import pytest

from ensembl.vep.common import CoordSystemLevel, Strand
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.models.variation.variation import Variation
from ensembl.vep.sources.feature_source import FeatureSource


def genomic_slice(seq_region_name: str = "1") -> Slice:
    """Create a Slice instance for testing."""
    return Slice(
        coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
        seq_region_name=seq_region_name,
        start=1,
        end=1000,
        strand=Strand.FORWARD,
    )


def create_variation_objects(length: int, vf_per_chr: int = 4) -> list[Variation]:
    """Function to create a set of test variation objects"""
    num_slices = length // vf_per_chr + 1
    slices = [genomic_slice(str(seq_region)) for seq_region in range(1, num_slices + 1)]
    var_list: list = []
    for value in range(length):
        slice_obj = slices[value // vf_per_chr]
        var_list.append(
            Variation(start=value, end=value, strand=Strand.FORWARD, slice=slice_obj, allele_string="A/T")
        )
    return var_list


@pytest.fixture
def variation_list() -> list[Variation]:
    """Create a list of Variation objects for testing."""
    return create_variation_objects(11)


@define(kw_only=True, auto_attribs=True)
class DummyVariationSource(FeatureSource):
    """Dummy VariationSource for the Concrete Buffer"""

    features: list[Variation]

    def get_all_features(self) -> list[Variation]:
        """Return all features from this source."""

        return self.features

    def iter_all_features(self) -> Iterator[Variation]:
        """Lazily yield all features from this source."""

        return iter(self.features)

    def get_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> list[Variation]:
        """Return a subset of the features."""

        return self.features[0:3]

    def iter_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> Iterator[Variation]:
        """Lazily yield a subset of the features."""

        return iter(self.features[5:8])

    def close(self) -> None:
        """Empty features list."""

        self.features = []
