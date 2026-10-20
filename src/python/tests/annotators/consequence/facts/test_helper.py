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

# pylint: disable=redefined-outer-name
"""Tests for consequence fact helpers."""

from __future__ import annotations

from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.facts.helper import variant_class_so_term
from ensembl.vep.common import ClassSOTerm, Strand
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import (
    StructuralVariation,
    StructuralVariationAllele,
    Variation,
    VariationAllele,
)


def _create_fact_context(
    transcript: Transcript, variation_allele: VariationAllele
) -> FactContext[TranscriptVariationAllele]:
    """Create a fact context for a transcript and variation allele."""
    transcript_variation_allele = TranscriptVariationAllele(
        feature=transcript,
        overlapping_feature=variation_allele,
    )
    return FactContext(overlap=transcript_variation_allele)


class TestCodonWindowHelper:
    """``codon_window`` is shared by the codon helpers and caches."""

    def test_codon_window_none_when_no_footprint(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """``codon_window`` returns ``None`` without a CDS footprint."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        mapper = transcript.cds_mapper
        assert mapper is not None
        assert mapper.codon_window(None) is None

    def test_codon_window_spans_whole_codons(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """``codon_window`` expands a footprint to whole codons."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        mapper = transcript.cds_mapper
        assert mapper is not None
        window = mapper.codon_window((4, 4))
        assert window == (3, 6)
        assert (window[1] - window[0]) % 3 == 0


class TestVariantClassSOTerm:
    """``variant_class_so_term`` follows the owning variation."""

    def test_returns_class_so_term_from_variation(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """The class SO term comes from the owning variation."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        variation = Variation(
            start=transcript.start,
            end=transcript.start,
            strand=Strand.FORWARD,
            slice=transcript.slice,
            allele_string="A/T",
        )
        variation_allele = variation.variation_alleles[0]
        fact_context = _create_fact_context(transcript, variation_allele)

        assert variant_class_so_term(fact_context) is ClassSOTerm.SNV

    def test_returns_none_for_an_allele_without_a_variation(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """A standalone variation allele has no class SO term."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        variation_allele = VariationAllele(
            start=transcript.start,
            end=transcript.start,
            strand=Strand.FORWARD,
            slice=transcript.slice,
            reference_allele="A",
            alternate_allele="T",
        )
        fact_context = _create_fact_context(transcript, variation_allele)

        assert variant_class_so_term(fact_context) is None

    def test_returns_class_so_term_from_structural_variation(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """A structural allele reads the term from its structural variation."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        structural_variation = StructuralVariation(
            start=transcript.start,
            end=transcript.end,
            strand=Strand.FORWARD,
            slice=transcript.slice,
            allele_string="-/DUP",
        )
        structural_variation.seed_attrs(class_so_term=ClassSOTerm.COPY_NUMBER_GAIN)
        structural_variation_allele = StructuralVariationAllele(
            start=transcript.start,
            end=transcript.end,
            strand=Strand.FORWARD,
            slice=transcript.slice,
            reference_allele="-",
            alternate_allele="DUP",
            variation=structural_variation,
        )
        fact_context = _create_fact_context(transcript, structural_variation_allele)

        assert variant_class_so_term(fact_context) is ClassSOTerm.COPY_NUMBER_GAIN
