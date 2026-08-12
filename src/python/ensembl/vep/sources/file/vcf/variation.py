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
"""VCF-back Variation feature source"""

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

@define(kw_only=True, auto_attribs=True)
class VCFVariationSource(FileFeatureSource[VCFParser, Variation]):
    """Variation source backed by a VCF file.

    Attributes:
        adaptor: The parser adaptor used to stream VCF records.
    """

    file_path: str | Path | None = field(default=None)
    adaptor: VCFParser = field(init=False, factory=VCFParser)

    def _records_to_features(self, records: list[VCFRecord]) -> list[Variation]:
        """Convert VCFRecords into Variation or StructuralVariation objects.

        Args:
            records (list[VCFRecord]): VCF records to convert into variation features.

        Returns:
            list[Variation]: Variation and StructuralVariation objects created from the VCF records.
        """

        features: list[Variation] = []
        for record in records:
            if self._is_structural_variation_record(record):
                features.extend(self._records_to_structural_variation([record]))
            else:
                features.extend(self._records_to_variation([record]))

        return features

    @staticmethod
    def _records_to_variation(records: list) -> list[Variation]:
        """Convert VCFRecords into Variation objects."""
        
        pass

    @staticmethod
    def _records_to_strutural_variation(records: list) -> list[StructuralVariation]:
        """Convert VCFRecords into StructuralVariation objects."""
        
        pass

    @staticmethod
    def _stable_id_from_vcf_record(record: VCFRecord) -> str | None:
        """Return a stable identifier from the VCF ID column, when present.

        Args:
            record (VCFRecord): VCF record containing ID values.

        Returns:
            str | None: Semicolon-delimited stable identifier, or None when the record has no IDs.
        """

        return ";".join(record.id) or None

    def _end_from_vcf_record(self, record: VCFRecord, is_symbolic: bool = False) -> int:
        """Return the VCF record end coordinate.

        Args:
            record (VCFRecord): VCF record containing position, reference allele, and optional END info.

        Returns:
            int: End coordinate from the SVLEN/END info field, or the reference-derived end coordinate.
        """

        if is_symbolic:
            lengths = record.info.get("SVLEN")
            if lengths is not None:
                length = abs(max([int(allele_length) for allele_length in lengths.split(",")], key=abs))
                return record.pos + length - 1

            end = self._int_info_value(record, "END")
            if end is not None:
                return end
            return record.pos

        return record.pos + len(record.reference) - 1

    @staticmethod
    def _int_info_value(record: VCFRecord, key: str) -> int | None:
        """Return an INFO value parsed as an integer.

        Args:
            record (VCFRecord): VCF record containing INFO values.
            key (str): INFO field key to parse.

        Returns:
            int | None: Parsed integer value, or None when the key is absent or not an integer.
        """

        value = record.info.get(key)
        if value is None:
            return None

        try:
            return int(value)
        except (TypeError, ValueError):
            return None
