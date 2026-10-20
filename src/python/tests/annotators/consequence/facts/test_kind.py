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
"""Tests for consequence classification facts."""

from __future__ import annotations

import pytest

from ensembl.vep.annotators.consequence.facts.base import REGISTRY, closure_cost, get_fact
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.common import CostTier, Strand
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import FeatureOverlap, TranscriptVariationAllele
from ensembl.vep.models.regulatory import Motif, Regulatory
from ensembl.vep.models.sequence import Slice
from ensembl.vep.models.variation import StructuralVariationAllele, VariationAllele

_CLASS_FACTS = ("is_transcript", "is_motif", "is_regulatory", "is_structural", "is_short")


def _create_variation_allele(transcript_slice: Slice | None, *, structural: bool = False) -> VariationAllele:
    """Create a one-base structural or ordinary variation allele."""
    allele_class = StructuralVariationAllele if structural else VariationAllele
    return allele_class(
        start=1051,
        end=1051,
        strand=Strand.FORWARD,
        slice=transcript_slice,
        reference_allele="A",
        alternate_allele="T",
    )


def _create_transcript_variation_allele(
    transcript: Transcript, *, start: int, end: int, reference: str, alternate: str
) -> TranscriptVariationAllele:
    """Create a transcript variation allele at explicit genomic coordinates."""
    variation_allele = VariationAllele(
        start=start,
        end=end,
        strand=Strand.FORWARD,
        slice=transcript.slice,
        reference_allele=reference,
        alternate_allele=alternate,
    )
    return TranscriptVariationAllele(feature=transcript, overlapping_feature=variation_allele)


class TestClassMembershipFactsAreFree:
    """The facts read Python types only, so they must cost nothing to compute."""

    def test_every_class_fact_is_registered_at_the_free_tier(self) -> None:
        """Every class-membership fact is registered at the FREE tier."""
        for name in _CLASS_FACTS:
            assert REGISTRY[name].cost == CostTier.FREE

    def test_every_class_fact_has_a_free_closure(self) -> None:
        """No class fact pulls in a dependency, so its closure cost is FREE too.

        This is what makes them usable as cheap guards: FactRule.cheap_guard only picks up a
        conjunct whose *closure* cost is at most BOUNDS.
        """
        for name in _CLASS_FACTS:
            assert closure_cost(name) == CostTier.FREE

    def test_no_class_fact_declares_a_dependency_or_requirement(self) -> None:
        """They are leaves of the fact DAG."""
        for name in _CLASS_FACTS:
            assert REGISTRY[name].deps == ()
            assert REGISTRY[name].requires == ()


class TestClassMembershipFacts:
    """Test facts that check class membership."""

    def test_transcript_is_true(self, synthetic_transcripts_by_stable_id: dict[str, Transcript]) -> None:
        """is_transcript is true for a Transcript overlap."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        overlap = TranscriptVariationAllele(
            feature=transcript, overlapping_feature=_create_variation_allele(transcript.slice)
        )
        context = FactContext(overlap=overlap)
        assert get_fact(context, "is_transcript") is True
        assert get_fact(context, "is_regulatory") is False
        assert get_fact(context, "is_motif") is False

    def test_motif_overlap_is_motif_only(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """A Motif overlap is a motif and *not* a regulatory feature."""
        transcript_slice = synthetic_transcripts_by_stable_id["ENST0001"].slice
        motif = Motif(start=1000, end=1100, strand=Strand.FORWARD, slice=transcript_slice)
        overlap = FeatureOverlap(
            feature=motif, overlapping_feature=_create_variation_allele(transcript_slice)
        )
        context = FactContext(overlap=overlap)
        assert get_fact(context, "is_motif") is True
        assert get_fact(context, "is_regulatory") is False
        assert get_fact(context, "is_transcript") is False

    def test_regulatory_overlap_is_regulatory_only(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """A Regulatory overlap is a regulatory feature and not a motif."""
        transcript_slice = synthetic_transcripts_by_stable_id["ENST0001"].slice
        regulatory = Regulatory(start=1000, end=1100, strand=Strand.FORWARD, slice=transcript_slice)
        overlap = FeatureOverlap(
            feature=regulatory, overlapping_feature=_create_variation_allele(transcript_slice)
        )
        context = FactContext(overlap=overlap)
        assert get_fact(context, "is_regulatory") is True
        assert get_fact(context, "is_motif") is False
        assert get_fact(context, "is_transcript") is False

    def test_is_short_is_true(self, synthetic_transcripts_by_stable_id: dict[str, Transcript]) -> None:
        """is_short is true for a short variant overlap."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        overlap = TranscriptVariationAllele(
            feature=transcript, overlapping_feature=_create_variation_allele(transcript.slice)
        )
        context = FactContext(overlap=overlap)
        assert get_fact(context, "is_short") is True
        assert get_fact(context, "is_structural") is False

    def test_is_structural_is_true(self, synthetic_transcripts_by_stable_id: dict[str, Transcript]) -> None:
        """is_structural is true for a structural variant overlap."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        overlap = TranscriptVariationAllele(
            feature=transcript,
            overlapping_feature=_create_variation_allele(transcript.slice, structural=True),
        )
        context = FactContext(overlap=overlap)
        assert get_fact(context, "is_structural") is True
        assert get_fact(context, "is_short") is False


class TestBiotypeFacts:
    """Biotype facts read the transcript's biotype, falling back to the gene's."""

    @pytest.mark.parametrize(
        ("stable_id", "expected_protein_coding", "expected_non_coding", "expected_nmd"),
        [
            ("ENST0001", True, False, False),
            ("ENST0002", True, False, False),
            ("ENST0003", False, True, False),
            ("ENST0004", False, False, True),
        ],
    )
    def test_biotype_facts(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        stable_id: str,
        expected_protein_coding: bool,
        expected_non_coding: bool,
        expected_nmd: bool,
    ) -> None:
        """Test that transcript biotypes set the expected facts."""
        transcript = synthetic_transcripts_by_stable_id[stable_id]
        overlap = _create_transcript_variation_allele(
            transcript, start=transcript.start, end=transcript.start, reference="A", alternate="T"
        )
        context = FactContext(overlap=overlap)
        assert get_fact(context, "biotype_protein_coding") is expected_protein_coding
        assert get_fact(context, "biotype_non_coding") is expected_non_coding
        assert get_fact(context, "biotype_nmd") is expected_nmd


class TestVariantShapeFacts:
    """FREE facts derived purely from the allele strings."""

    @pytest.mark.parametrize(
        (
            "start",
            "end",
            "reference",
            "alternate",
            "expected_is_snv",
            "expected_is_insertion",
            "expected_is_deletion",
            "expected_is_short_insertion",
            "expected_is_short_deletion",
            "expected_is_structural_insertion",
            "expected_is_structural_deletion",
            "expected_cds_length_change_inframe",
        ),
        [
            (1051, 1051, "A", "T", True, False, False, None, None, False, False, True),
            (1051, 1050, "", "AC", False, True, False, True, False, False, False, False),
            (1051, 1050, "", "ACC", False, True, False, True, None, False, None, True),
            (1051, 1052, "AA", "", False, False, True, None, True, None, False, False),
            (1051, 1053, "AAA", "", False, False, True, None, True, None, False, True),
        ],
    )
    def test_shape_facts(
        self,
        synthetic_transcripts_by_stable_id: dict[str, Transcript],
        start: int,
        end: int,
        reference: str,
        alternate: str,
        expected_is_snv: bool,
        expected_is_insertion: bool,
        expected_is_deletion: bool,
        expected_is_short_insertion: bool | None,
        expected_is_short_deletion: bool | None,
        expected_is_structural_insertion: bool | None,
        expected_is_structural_deletion: bool | None,
        expected_cds_length_change_inframe: bool,
    ) -> None:
        """Variant-shape facts behave consistently across the supported allele patterns."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        overlap = _create_transcript_variation_allele(
            transcript, start=start, end=end, reference=reference, alternate=alternate
        )
        context = FactContext(overlap=overlap)

        assert get_fact(context, "is_snv") is expected_is_snv
        assert get_fact(context, "is_insertion") is expected_is_insertion
        assert get_fact(context, "is_deletion") is expected_is_deletion

        if expected_is_short_insertion is not None:
            assert get_fact(context, "is_short_insertion") is expected_is_short_insertion
        if expected_is_short_deletion is not None:
            assert get_fact(context, "is_short_deletion") is expected_is_short_deletion
        if expected_is_structural_insertion is not None:
            assert get_fact(context, "is_structural_insertion") is expected_is_structural_insertion
        if expected_is_structural_deletion is not None:
            assert get_fact(context, "is_structural_deletion") is expected_is_structural_deletion

        assert get_fact(context, "cds_length_change_inframe") is expected_cds_length_change_inframe

    def test_dash_alleles_treated_as_empty(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """``-`` is a common real-world spelling of an empty allele; facts must not choke on it."""
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        overlap = _create_transcript_variation_allele(
            transcript, start=1051, end=1050, reference="-", alternate="AC"
        )
        context = FactContext(overlap=overlap)
        assert get_fact(context, "is_insertion") is True
