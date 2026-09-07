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

from typing import Iterator
from attrs import define
import pytest

from ensembl.vep.common import Strand
from ensembl.vep.models import Feature
from ensembl.vep.models.variation.variation import Variation
from ensembl.vep.sources.feature_source import FeatureSource


def create_feature_objects(length: int) -> list[Feature]:
    """Function to create a set of test features"""
    return([Feature(start=value, end=value, strand=Strand.FORWARD) for value in range(length)])

@pytest.fixture
def features_list():
    return create_feature_objects(11)

def create_variation_objects(length: int) -> list[Feature]:
    """Function to create a set of test variation objects"""
    return([Variation(start=value, end=value, strand=Strand.FORWARD, allele_string="A/T") for value in range(length)])

@pytest.fixture
def variation_list():
    return create_variation_objects(11)

@define(kw_only=True, auto_attribs=True)
class DummyFeatureSource(FeatureSource):
    """Dummy FeatureSource for the Concrete Buffer"""

    features: list[Feature]

    def get_all_features(self) -> list[Feature]:
        """Return all features from this source."""

        return self.features

    def iter_all_features(self) -> Iterator[Feature]:
        """Lazily yield all features from this source."""

        return iter(self.features)

    def get_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> list[Feature]:
        """Return a subset of the features."""

        return self.features[0:3]

    def iter_features_by_region(
        self,
        region: str | None = None,
        *,
        seq_region_name: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> Iterator[Feature]:
        """Lazily yield a subset of the features."""

        return iter(self.features[5:8])

    def close(self) -> None:
        """Empty features list."""

        self.features = []
