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
"""Tests for the JSON annotation stream factory."""

from __future__ import annotations

from ensembl.vep.output.factory import AnnotatedVariant
from ensembl.vep.output.factory.stream import JSONStreamFactory


class TestJSONStreamFactory:
    """Per-variant JSON schema conversion."""

    def test_record_writes_locus_input_and_transcript_consequences(self) -> None:
        """Build an Ensembl VEP record from a transcript output hash."""
        variant = AnnotatedVariant(
            id="rs1",
            seq_region_name="1",
            start=10,
            end=10,
            strand=1,
            allele_string="A/G",
            assembly_name="GRCh38",
            input_columns=("1", "10", "rs1", "A", "G", ".", ".", "."),
            most_severe_consequence="missense_variant",
            variant_class="SNV",
            annotations=(
                {
                    "Allele": "G",
                    "Amino_acids": "K/R",
                    "BIOTYPE": "protein_coding",
                    "CDS_position": "4-6",
                    "Consequence": ["missense_variant"],
                    "Feature": "ENST000001",
                    "Feature_type": "Transcript",
                    "Gene": "ENSG000001",
                    "IMPACT": "MODERATE",
                    "Protein_position": "2",
                    "STRAND": 1,
                    "SYMBOL": "GENE1",
                    "Uploaded_variation": "rs1",
                    "VARIANT_CLASS": "SNV",
                },
            ),
        )

        record = JSONStreamFactory().record(variant)

        assert record == {
            "id": "rs1",
            "seq_region_name": "1",
            "start": 10,
            "end": 10,
            "strand": 1,
            "allele_string": "A/G",
            "input": "1\t10\trs1\tA\tG\t.\t.\t.",
            "assembly_name": "GRCh38",
            "most_severe_consequence": "missense_variant",
            "variant_class": "SNV",
            "transcript_consequences": [
                {
                    "variant_allele": "G",
                    "amino_acids": "K/R",
                    "biotype": "protein_coding",
                    "cds_start": 4,
                    "cds_end": 6,
                    "consequence_terms": ["missense_variant"],
                    "transcript_id": "ENST000001",
                    "gene_id": "ENSG000001",
                    "impact": "MODERATE",
                    "protein_start": 2,
                    "protein_end": 2,
                    "strand": 1,
                    "gene_symbol": "GENE1",
                }
            ],
        }

    def test_record_groups_each_supported_feature_type(self) -> None:
        """Keep each feature type in its VEP JSON consequence collection."""
        variant = AnnotatedVariant(
            id=None,
            seq_region_name="1",
            start=10,
            end=10,
            strand=1,
            allele_string="A/G",
            annotations=(
                {
                    "Allele": "G",
                    "Consequence": ["regulatory_region_variant"],
                    "Feature": "ENSR000001",
                    "Feature_type": "RegulatoryFeature",
                },
                {
                    "Allele": "G",
                    "Consequence": ["TF_binding_site_variant"],
                    "Feature": "ENSM000001",
                    "Feature_type": "MotifFeature",
                },
                {
                    "Allele": "G",
                    "Consequence": ["intergenic_variant"],
                    "Feature": "intergenic",
                    "Feature_type": "Intergenic",
                },
                {"Allele": "G", "Consequence": ["intergenic_variant"]},
                {"Allele": "G", "Consequence": ["other"], "Feature_type": "Other"},
            ),
        )

        record = JSONStreamFactory().record(variant)

        assert record["regulatory_feature_consequences"] == [
            {
                "variant_allele": "G",
                "consequence_terms": ["regulatory_region_variant"],
                "regulatory_feature_id": "ENSR000001",
            }
        ]
        assert record["motif_feature_consequences"] == [
            {
                "variant_allele": "G",
                "consequence_terms": ["TF_binding_site_variant"],
                "motif_feature_id": "ENSM000001",
            }
        ]
        assert record["intergenic_consequences"] == [
            {
                "variant_allele": "G",
                "consequence_terms": ["intergenic_variant"],
                "feature_id": "intergenic",
            },
            {"variant_allele": "G", "consequence_terms": ["intergenic_variant"]},
        ]

    def test_record_omits_empty_values_but_keeps_a_deletion_allele(self) -> None:
        """Discard VEP's empty values without losing an empty alternate allele."""
        variant = AnnotatedVariant(
            id=None,
            seq_region_name="1",
            start=10,
            end=11,
            strand=1,
            allele_string="AC/-",
            annotations=(
                {
                    "Allele": "-",
                    "Consequence": ["inframe_deletion"],
                    "Feature_type": "Transcript",
                    "Codons": "-",
                    "EXON": "",
                    "Existing_variation": None,
                    "custom_flag": "on",
                    "cDNA_position": "?-9",
                },
            ),
        )

        record = JSONStreamFactory().record(variant)

        assert record["transcript_consequences"] == [
            {
                "variant_allele": "-",
                "consequence_terms": ["inframe_deletion"],
                "custom_flag": "on",
                "cdna_end": 9,
            }
        ]

    def test_record_preserves_position_totals_and_converts_yes_flags(self) -> None:
        """Match VEP's scalar conversion for totals and enabled flags."""
        variant = AnnotatedVariant(
            id=None,
            seq_region_name="1",
            start=10,
            end=10,
            strand=1,
            allele_string="A/G",
            annotations=(
                {
                    "Allele": "G",
                    "Consequence": ["missense_variant"],
                    "Feature_type": "Transcript",
                    "cDNA_position": "12/100",
                    "CANONICAL": "YES",
                    "REF_ALLELE": "A",
                    "OVERLAPBP": 3,
                    "OVERLAPPC": 50.0,
                },
            ),
        )

        record = JSONStreamFactory().record(variant)

        assert record["transcript_consequences"] == [
            {
                "variant_allele": "G",
                "consequence_terms": ["missense_variant"],
                "cdna_start": "12/100",
                "cdna_end": "12/100",
                "canonical": 1,
                "ref_allele": "A",
                "bp_overlap": 3,
                "percentage_overlap": 50.0,
            }
        ]

    def test_record_structures_predictions_domains_and_numeric_values(self) -> None:
        """Apply the nested conversions used by Ensembl VEP JSON output."""
        variant = AnnotatedVariant(
            id="123",
            seq_region_name="1",
            start=10,
            end=10,
            strand=1,
            allele_string="A/G",
            annotations=(
                {
                    "Allele": "G",
                    "Consequence": ["missense_variant"],
                    "Feature": "000123",
                    "Feature_type": "Transcript",
                    "Gene": "000456",
                    "DOMAINS": ["Pfam:PF00001", "invalid"],
                    "SIFT": "deleterious(0.01)",
                    "PolyPhen": "0.875",
                    "custom_score": "2.5",
                },
            ),
        )

        record = JSONStreamFactory().record(variant)

        assert record["id"] == "123"
        assert record["seq_region_name"] == "1"
        assert record["transcript_consequences"] == [
            {
                "variant_allele": "G",
                "consequence_terms": ["missense_variant"],
                "transcript_id": "000123",
                "gene_id": "000456",
                "domains": [{"db": "Pfam", "name": "PF00001"}],
                "sift_prediction": "deleterious",
                "sift_score": 0.01,
                "polyphen_score": 0.875,
                "custom_score": 2.5,
            }
        ]

    def test_skipped_record_carries_only_its_input(self) -> None:
        """Do not expose partially annotated fields when annotation was skipped."""
        variant = AnnotatedVariant(
            id="rs1",
            seq_region_name="1",
            start=10,
            end=10,
            strand=1,
            allele_string="A/G",
            input_columns=("1", "10", "rs1", "A", "G"),
            skipped=True,
        )

        assert JSONStreamFactory().record(variant) == {"input": "1\t10\trs1\tA\tG"}

    def test_records_converts_an_iterable_lazily(self) -> None:
        """Yield one JSON record for each input variant."""
        variants = (
            AnnotatedVariant(
                id=None,
                seq_region_name="1",
                start=position,
                end=position,
                strand=1,
                allele_string="A/G",
            )
            for position in (10, 20)
        )

        records = JSONStreamFactory().records(variants)

        assert next(records)["start"] == 10
        assert next(records)["start"] == 20
