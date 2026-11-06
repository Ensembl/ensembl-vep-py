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
"""Ensembl VEP JSON-schema record generation."""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator
from typing import cast

from ensembl.vep.common.constants import EMPTY_ALLELE, UNKNOWN_VALUE
from ensembl.vep.output.factory.base import AnnotatedVariant, AnnotationRecord

__all__ = ["JSONStreamFactory"]

_FEATURE_GROUPS = {
    "Transcript": "transcript_consequences",
    "RegulatoryFeature": "regulatory_feature_consequences",
    "MotifFeature": "motif_feature_consequences",
    "Intergenic": "intergenic_consequences",
}
_CORE_FIELD_NAMES = {
    "Allele": "variant_allele",
    "Amino_acids": "amino_acids",
    "BIOTYPE": "biotype",
    "Codons": "codons",
    "Consequence": "consequence_terms",
    "DISTANCE": "distance",
    "EXON": "exon",
    "Existing_variation": "existing_variation",
    "Feature": "feature_id",
    "Gene": "gene_id",
    "HGNC_ID": "hgnc_id",
    "IMPACT": "impact",
    "INTRON": "intron",
    "OVERLAPBP": "bp_overlap",
    "OVERLAPPC": "percentage_overlap",
    "PICK": "pick",
    "REFSEQ": "refseq_transcript_ids",
    "REF_ALLELE": "ref_allele",
    "STRAND": "strand",
    "SYMBOL": "gene_symbol",
    "SYMBOL_SOURCE": "gene_symbol_source",
    "ENSP": "protein_id",
}
_POSITION_FIELD_NAMES = {
    "cDNA_position": "cdna",
    "CDS_position": "cds",
    "Protein_position": "protein",
}
_SKIPPED_FIELDS = frozenset({"Feature_type", "Location", "Uploaded_variation", "VARIANT_CLASS"})
_NUMBERIFY_EXEMPT = frozenset(
    {
        "seq_region_name",
        "id",
        "gene_id",
        "gene_symbol",
        "transcript_id",
        "SYMBOL",
    }
)
_PREDICTION_PATTERN = re.compile(
    r"(?:(?P<prediction>[A-Za-z_]+)(?:\((?P<score>[\d.]+)\))?|(?P<score_only>[\d.]+))"
)


class JSONStreamFactory:
    """Ensembl VEP JSON-schema record factory."""

    def records(self, variants: Iterable[AnnotatedVariant]) -> Iterator[dict[str, object]]:
        """Yield each variant as a JSON-compatible dictionary."""
        for variant in variants:
            yield self.record(variant)

    def record(self, variant: AnnotatedVariant) -> dict[str, object]:
        """Return one Ensembl VEP JSON-schema dictionary."""
        record: dict[str, object] = {
            "id": variant.id,
            "seq_region_name": variant.seq_region_name,
            "start": variant.start,
            "end": variant.end,
            "strand": variant.strand,
            "allele_string": variant.allele_string,
            "assembly_name": variant.assembly_name,
            "most_severe_consequence": variant.most_severe_consequence or UNKNOWN_VALUE,
        }
        if variant.input_columns:
            record["input"] = "\t".join(variant.input_columns)
        if variant.variant_class is not None:
            record["variant_class"] = variant.variant_class
        if variant.skipped:
            return {"input": record.get("input", "")}

        grouped_annotations: dict[str, list[dict[str, object]]] = {}
        for annotation in variant.annotations:
            feature_type = annotation.get("Feature_type")
            group_name = _FEATURE_GROUPS.get(str(feature_type))
            if feature_type is None:
                group_name = "intergenic_consequences"
            if group_name is None:
                continue
            grouped_annotations.setdefault(group_name, []).append(self._annotation(annotation))
        record.update(grouped_annotations)
        return cast(dict[str, object], self._numberify(record))

    def _annotation(self, annotation: AnnotationRecord) -> dict[str, object]:
        """Convert one output hash to its nested JSON representation."""
        converted: dict[str, object] = {}
        feature_type = annotation.get("Feature_type")
        for name, value in annotation.items():
            if name in _SKIPPED_FIELDS or self._empty(name, value):
                continue
            if name in _POSITION_FIELD_NAMES:
                self._add_position(converted, _POSITION_FIELD_NAMES[name], value)
                continue
            json_name = _CORE_FIELD_NAMES.get(name, name.lower())
            if name == "Feature":
                json_name = self._feature_name(feature_type)
            converted[json_name] = 1 if value == "YES" else value
        self._split_predictions(converted)
        self._structure_domains(converted)
        return converted

    @staticmethod
    def _split_predictions(record: dict[str, object]) -> None:
        """Split SIFT and PolyPhen values into prediction and score fields."""
        for tool in ("sift", "polyphen"):
            value = record.pop(tool, None)
            if value is None:
                continue
            match = _PREDICTION_PATTERN.fullmatch(str(value))
            if match is None:
                continue
            prediction = match.group("prediction")
            score = match.group("score") or match.group("score_only")
            if prediction:
                record[f"{tool}_prediction"] = prediction
            if score is not None:
                record[f"{tool}_score"] = score

    @staticmethod
    def _structure_domains(record: dict[str, object]) -> None:
        """Convert domain names into Ensembl VEP's database/name objects."""
        domains = record.get("domains")
        if not isinstance(domains, (list, tuple)):
            return
        structured = []
        for domain in domains:
            database, separator, name = str(domain).partition(":")
            if separator and database and name:
                structured.append({"db": database, "name": name})
        record["domains"] = structured

    @staticmethod
    def _feature_name(feature_type: object) -> str:
        """Return the JSON identifier name for one feature type."""
        if feature_type == "Transcript":
            return "transcript_id"
        if feature_type == "RegulatoryFeature":
            return "regulatory_feature_id"
        if feature_type == "MotifFeature":
            return "motif_feature_id"
        return "feature_id"

    @staticmethod
    def _empty(name: str, value: object) -> bool:
        """Return whether a field is absent from Perl VEP JSON output."""
        return value is None or value == "" or (name != "Allele" and value == EMPTY_ALLELE)

    @staticmethod
    def _add_position(record: dict[str, object], prefix: str, value: object) -> None:
        """Add one VEP position field as JSON start and end fields."""
        position = str(value)
        start, separator, end = position.partition("-")
        start_value = JSONStreamFactory._position_value(start)
        end_value = JSONStreamFactory._position_value(end if separator else start)
        if start_value is not None:
            record[f"{prefix}_start"] = start_value
        if end_value is not None:
            record[f"{prefix}_end"] = end_value

    @staticmethod
    def _position_value(value: str) -> int | str | None:
        """Return a numeric coordinate or a VEP position-with-total string."""
        if value.isdigit():
            return int(value)
        position, separator, total = value.partition("/")
        if separator and position.isdigit() and total.isdigit():
            return value
        return None

    @classmethod
    def _numberify(cls, value: object, *, key: str | None = None) -> object:
        """Recursively convert numeric strings as Ensembl VEP JSON does."""
        if isinstance(value, dict):
            return {name: cls._numberify(item, key=name) for name, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [cls._numberify(item) for item in value]
        if key in _NUMBERIFY_EXEMPT or not isinstance(value, str):
            return value
        try:
            return int(value)
        except ValueError:
            try:
                return float(value)
            except ValueError:
                return value
