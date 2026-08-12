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
"""Vcf-back Variation feature source"""

from typing import Generic, TypeVar
from pathlib import Path

from attrs import define, field

from ensembl.vep.models.variation.variation import Variation, VariationAllele
from ensembl.vep.models.variation.structural_variation import StructuralVariation, StructuralVariationAllele
from ensembl.vep.parsers.vcf import VCFParser
from ensembl.vep.sources.file.file_feature_source import FileFeatureSource
from ensembl.vep.ensembl_vep_core import (  # noqa: F401,F403 pylint: disable=E0401,E0611
    VCFRecord
)

V = TypeVar("V", bound=Variation)

@define(kw_only=True, auto_attribs=True)
class VcfVariationSource(FileFeatureSource[VCFParser, V]):
    """Variation source backed by a Vcf file.

    Attributes:
        adaptor: The parser adaptor used to stream Vcf records.
    """

    file_path: str | Path | None = field(default=None)
    adaptor: VCFParser = field(init=False, factory=VCFParser)

    def _records_to_features(self, records: list) -> list[V]:
        """Convert VcfRecords into Variation objects."""
        
        pass

    @staticmethod
    def _records_to_variation(records: list) -> list[Variation]:
        """Convert VcfRecords into Variation objects."""
        
        pass

    @staticmethod
    def _records_to_strutural_variation(records: list) -> list[StructuralVariation]:
        """Convert VcfRecords into StructuralVariation objects."""
        
        pass