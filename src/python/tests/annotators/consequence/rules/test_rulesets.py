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
"""Tests for fact rule-set construction, guards, pruning, validation, and subsumption."""

from __future__ import annotations

import pytest

from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.consequence.rules.catalogue import DEFAULT_TERM, EXCLUDED_TERMS
from ensembl.vep.annotators.consequence.rules.rulesets import (
    FactRuleSet,
    build_fact_ruleset,
    prune,
    subsumes_topological_order,
)
from ensembl.vep.common import ClassSOTerm, Strand
from ensembl.vep.models import Feature
from ensembl.vep.models.genic import Transcript
from ensembl.vep.models.overlap import FeatureOverlap, TranscriptVariationAllele
from ensembl.vep.models.regulatory import Motif, Regulatory
from ensembl.vep.models.variation import (
    StructuralVariation,
    StructuralVariationAllele,
    VariationAllele,
)

_FEATURE_START = 1000
_FEATURE_END = 1100


def _create_feature_overlap(feature: Feature, variation_allele: VariationAllele) -> FeatureOverlap:
    """Create the appropriate overlap type for a feature and variation allele."""
    if isinstance(feature, Transcript):
        return TranscriptVariationAllele(feature=feature, overlapping_feature=variation_allele)
    return FeatureOverlap(feature=feature, overlapping_feature=variation_allele)


def _evaluate_consequence_terms(feature: Feature, *, start: int, end: int) -> list[str]:
    """Return consequence terms for a one-base variation allele."""
    variation_allele = VariationAllele(
        start=start,
        end=end,
        strand=Strand.FORWARD,
        slice=feature.slice,
        reference_allele="A",
        alternate_allele="T",
    )
    context = FactContext(overlap=_create_feature_overlap(feature, variation_allele))
    ruleset = build_fact_ruleset()
    return [rule.term for rule in ruleset.rules if rule.when.evaluate(context)]


@pytest.fixture
def motif(synthetic_transcripts_by_stable_id: dict[str, Transcript]) -> Motif:
    """A motif feature spanning 1000-1100 on the synthetic slice."""
    return Motif(
        start=_FEATURE_START,
        end=_FEATURE_END,
        strand=Strand.FORWARD,
        slice=synthetic_transcripts_by_stable_id["ENST0001"].slice,
    )


@pytest.fixture
def regulatory(synthetic_transcripts_by_stable_id: dict[str, Transcript]) -> Regulatory:
    """A regulatory feature spanning 1000-1100 on the synthetic slice."""
    return Regulatory(
        start=_FEATURE_START,
        end=_FEATURE_END,
        strand=Strand.FORWARD,
        slice=synthetic_transcripts_by_stable_id["ENST0001"].slice,
    )


def _evaluate_allele_consequence_terms(
    feature: Feature,
    *,
    start: int,
    end: int,
    reference: str,
    alternate: str,
    structural: bool,
    class_so_term: ClassSOTerm | None = None,
) -> list[str]:
    """Return consequence terms for an explicitly configured variation allele."""
    variation_allele: VariationAllele
    if structural:
        variation = StructuralVariation(
            start=start,
            end=end,
            strand=Strand.FORWARD,
            slice=feature.slice,
            allele_string=f"{reference}/{alternate}",
        )
        if class_so_term is not None:
            variation.seed_attrs(class_so_term=class_so_term)
        variation_allele = StructuralVariationAllele(
            start=start,
            end=end,
            strand=Strand.FORWARD,
            slice=feature.slice,
            reference_allele=reference,
            alternate_allele=alternate,
            variation=variation,
        )
    else:
        variation_allele = VariationAllele(
            start=start,
            end=end,
            strand=Strand.FORWARD,
            slice=feature.slice,
            reference_allele=reference,
            alternate_allele=alternate,
        )
    context = FactContext(overlap=_create_feature_overlap(feature, variation_allele))
    ruleset = build_fact_ruleset()
    return [rule.term for rule in ruleset.rules if rule.when.evaluate(context)]


def _spanning_deletion_terms(feature: Feature, *, structural: bool) -> list[str]:
    """Return terms for a deletion spanning the feature and one flanking base."""
    length = _FEATURE_END - _FEATURE_START + 3
    return _evaluate_allele_consequence_terms(
        feature,
        start=_FEATURE_START - 1,
        end=_FEATURE_END + 1,
        reference="A" * length,
        alternate="-",
        structural=structural,
        class_so_term=ClassSOTerm.DELETION if structural else None,
    )


def _spanning_insertion_terms(feature: Feature, *, structural: bool) -> list[str]:
    """Return terms for an insertion spanning the feature and one flanking base."""
    length = _FEATURE_END - _FEATURE_START + 3
    return _evaluate_allele_consequence_terms(
        feature,
        start=_FEATURE_START - 1,
        end=_FEATURE_END + 1,
        reference="A" * length,
        alternate="A" * (length + 10),
        structural=structural,
        class_so_term=ClassSOTerm.COPY_NUMBER_GAIN if structural else None,
    )


def _contained_insertion_terms(feature: Feature, *, structural: bool) -> list[str]:
    """Return terms for an insertion contained within the feature."""
    return _evaluate_allele_consequence_terms(
        feature,
        start=_FEATURE_START + 10,
        end=_FEATURE_START + 9,
        reference="",
        alternate="AAAAA",
        structural=structural,
        class_so_term=ClassSOTerm.COPY_NUMBER_GAIN if structural else None,
    )


class TestBuildFactRuleset:
    """Compiling the catalogue into a rule set."""

    def test_build_fact_ruleset_covers_every_term_with_a_when(
        self,
        transcript_ruleset: FactRuleSet,
    ) -> None:
        """The rule set is feature-class agnostic: every catalogue term with a ``when`` is compiled.

        A term restricted to a feature class carries the class fact as a conjunct in its own formula,
        so every rule is evaluated against every overlap and settles itself.
        """
        assert len(transcript_ruleset.rules) == 39
        assert not set(transcript_ruleset.by_term) & EXCLUDED_TERMS

    def test_class_restricted_terms_carry_their_class_as_a_conjunct(
        self,
        transcript_ruleset: FactRuleSet,
    ) -> None:
        """Class-restricted terms carry their class as a conjunct in their own formula."""
        assert "is_motif" in transcript_ruleset.by_term["TF_binding_site_variant"].cheap_guard()
        assert "is_regulatory" in transcript_ruleset.by_term["regulatory_region_variant"].cheap_guard()
        assert "is_transcript" in transcript_ruleset.by_term["missense_variant"].cheap_guard()
        # ... and the classes are mutually exclusive, so neither fires on the other's feature.
        assert "is_motif" not in transcript_ruleset.by_term["regulatory_region_variant"].atoms()
        assert "is_regulatory" not in transcript_ruleset.by_term["TF_binding_site_variant"].atoms()

    def test_rules_are_rank_ordered(self, transcript_ruleset: FactRuleSet) -> None:
        """Rules are rank ordered."""
        ranks = [rule.rank for rule in transcript_ruleset.rules]
        assert ranks == sorted(ranks)

    def test_by_term_lookup_returns_the_corresponding_rule(self, transcript_ruleset: FactRuleSet) -> None:
        """Test that by-term lookup returns the corresponding rule."""
        assert transcript_ruleset.by_term["stop_gained"].accession == "SO:0001587"
        assert transcript_ruleset.by_term["missense_variant"].rank == 13


class TestCheapGuard:
    """Cheap rule-guard construction."""

    def test_stop_gained_guard_does_not_contain_expensive_atom(self, transcript_ruleset: FactRuleSet) -> None:
        """Stop-gained's guard does not contain the expensive atom ``overlap_cds``."""
        rule = transcript_ruleset.by_term["stop_gained"]
        assert frozenset({"is_transcript", "overlap_feature", "overlap_cds"}) <= rule.atoms()
        assert rule.cheap_guard() == frozenset({"is_transcript", "overlap_feature"})
        assert "overlap_cds" not in rule.cheap_guard()

    def test_guard_atoms_are_a_subset_of_the_rule_atoms(self, transcript_ruleset: FactRuleSet) -> None:
        """Guard atoms are a subset of the rule atoms."""
        for rule in transcript_ruleset.rules:
            assert rule.cheap_guard() <= rule.atoms()


class TestToOverlapConsequence:
    """Rule-to-consequence conversion."""

    def test_round_trips_catalogue_metadata(self, transcript_ruleset: FactRuleSet) -> None:
        """Catalogue metadata round-trips."""
        rule = transcript_ruleset.by_term["missense_variant"]
        consequence = rule.to_overlap_consequence()
        assert consequence.term == "missense_variant"
        assert consequence.accession == "SO:0001583"
        assert consequence.rank == 13
        assert consequence.is_default is False


class TestPrune:
    """Pruning of subsumed consequences."""

    def test_drops_subsumed_terms(self, transcript_ruleset: FactRuleSet) -> None:
        """Drops subsumed terms."""
        hits = [
            "stop_gained",
            "coding_sequence_variant",
            "protein_altering_variant",
            "coding_transcript_variant",
        ]
        assert prune(transcript_ruleset, hits) == ["stop_gained"]

    def test_keeps_terms_that_are_not_subsumed(self, transcript_ruleset: FactRuleSet) -> None:
        """Keeps terms that are not subsumed."""
        hits = ["missense_variant", "splice_region_variant"]
        assert prune(transcript_ruleset, hits, apply_default=False) == hits

    def test_apply_default_true_substitutes_default_term_when_empty(
        self, transcript_ruleset: FactRuleSet
    ) -> None:
        """Apply default true substitutes default term when empty."""
        assert prune(transcript_ruleset, [], apply_default=True) == [DEFAULT_TERM]

    def test_apply_default_false_returns_empty_list(self, transcript_ruleset: FactRuleSet) -> None:
        """Apply default false returns empty list."""
        assert prune(transcript_ruleset, [], apply_default=False) == []

    @pytest.mark.parametrize("splice_term", ["splice_acceptor_variant", "splice_donor_variant"])
    @pytest.mark.parametrize(
        "coding_term",
        [
            "stop_gained",
            "stop_lost",
            "stop_retained_variant",
            "start_lost",
            "start_retained_variant",
            "frameshift_variant",
            "missense_variant",
            "synonymous_variant",
            "inframe_deletion",
            "inframe_insertion",
        ],
    )
    def test_a_donor_or_acceptor_variant_hides_a_coding_term(
        self, transcript_ruleset: FactRuleSet, splice_term: str, coding_term: str
    ) -> None:
        """A variant over the two bases of a splice site is reported as that, not as a coding change."""
        assert prune(transcript_ruleset, [splice_term, coding_term], apply_default=False) == [splice_term]

    def test_two_terms_that_mutually_subsume_would_never_both_survive(
        self, transcript_ruleset: FactRuleSet
    ) -> None:
        """Two terms that mutually subsume would never both survive."""
        hits = ["splice_donor_region_variant", "splice_region_variant"]
        assert prune(transcript_ruleset, hits, apply_default=False) == ["splice_donor_region_variant"]


class TestTopologicalOrder:
    """The order itself."""

    def test_every_rule_appears_after_everything_that_subsumes_it(self) -> None:
        """The defining property of the order."""
        ruleset = build_fact_ruleset()
        position = {rule.term: index for index, rule in enumerate(subsumes_topological_order(ruleset))}
        for rule in ruleset.rules:
            for subsumed_term in rule.subsumes:
                if subsumed_term in position:
                    assert (
                        position[rule.term] < position[subsumed_term]
                    ), f"{rule.term} subsumes {subsumed_term} but is ordered after it"

    def test_the_order_contains_every_rule_exactly_once(self) -> None:
        """A topological sort is a permutation, not a filter."""
        ruleset = build_fact_ruleset()
        ordered = subsumes_topological_order(ruleset)
        assert len(ordered) == len(ruleset.rules)
        assert {rule.term for rule in ordered} == {rule.term for rule in ruleset.rules}

    def test_the_order_is_not_rank_order(self) -> None:
        """The catalogue's rank inversions are real, so rank order would skip the wrong rules.

        splice_donor_region_variant (rank 17) suppresses splice_region_variant (rank 16). Ordering
        by rank would visit the suppressed term first, and any skip decision taken on that basis
        would be backwards.
        """
        ruleset = build_fact_ruleset()
        ordered = subsumes_topological_order(ruleset)
        assert [rule.term for rule in ordered] != [rule.term for rule in ruleset.rules]

        inversions = {
            (rule.term, subsumed_term)
            for rule in ruleset.rules
            for subsumed_term in rule.subsumes
            if subsumed_term in ruleset.by_term and ruleset.by_term[subsumed_term].rank < rule.rank
        }
        assert ("splice_donor_region_variant", "splice_region_variant") in inversions


class TestFeatureKindGates:
    """Which terms fire for each feature kind, given the one shared rule set."""

    def test_a_transcript_overlap_never_fires_either_regulatory_term(
        self, synthetic_transcripts_by_stable_id: dict[str, Transcript]
    ) -> None:
        """Both regulatory terms require a Motif or Regulatory feature, never a Transcript.

        This is what keeps conformance on the transcript corpora byte-identical: adding a term
        restricted away from Transcript cannot change what a transcript overlap is annotated with.
        """
        transcript = synthetic_transcripts_by_stable_id["ENST0001"]
        terms = _evaluate_consequence_terms(transcript, start=1050, end=1050)
        assert "TF_binding_site_variant" not in terms
        assert "regulatory_region_variant" not in terms

    def test_the_motif_overlap_fires_the_motif_term_only(self, motif: Motif) -> None:
        """TF_binding_site_variant fires for a Motif; the regulatory-region term does not."""
        terms = _evaluate_consequence_terms(motif, start=1050, end=1050)
        assert "TF_binding_site_variant" in terms
        assert "regulatory_region_variant" not in terms

    def test_the_regulatory_overlap_fires_the_regulatory_term_only(self, regulatory: Regulatory) -> None:
        """regulatory_region_variant fires for a Regulatory feature; the motif term does not."""
        terms = _evaluate_consequence_terms(regulatory, start=1050, end=1050)
        assert "regulatory_region_variant" in terms
        assert "TF_binding_site_variant" not in terms


class TestTfBindingSiteVariant:
    """TF_binding_site_variant fires on any overlap with the motif feature."""

    def test_fires_inside_the_motif(self, motif: Motif) -> None:
        """A variant well inside the motif fires the term."""
        assert "TF_binding_site_variant" in _evaluate_consequence_terms(motif, start=1050, end=1050)

    @pytest.mark.parametrize("position", [_FEATURE_START, _FEATURE_END])
    def test_fires_on_each_boundary_base(self, motif: Motif, position: int) -> None:
        """The first and last base of the feature are inside it -- overlap() is inclusive."""
        assert "TF_binding_site_variant" in _evaluate_consequence_terms(motif, start=position, end=position)

    @pytest.mark.parametrize("position", [_FEATURE_START - 1, _FEATURE_END + 1])
    def test_does_not_fire_one_base_outside(self, motif: Motif, position: int) -> None:
        """One base beyond either end is not an overlap."""
        assert "TF_binding_site_variant" not in _evaluate_consequence_terms(
            motif, start=position, end=position
        )

    def test_fires_for_a_variant_that_merely_straddles_the_start(self, motif: Motif) -> None:
        """Partial overlap is enough: VEP's within_feature is overlap(), not containment.

        This is the case that would fail if the rule had been read as a genuine "within" test.
        """
        assert "TF_binding_site_variant" in _evaluate_consequence_terms(
            motif, start=_FEATURE_START - 10, end=_FEATURE_START
        )


class TestRegulatoryRegionVariant:
    """regulatory_region_variant behaves identically over a Regulatory feature."""

    def test_fires_inside_the_regulatory_feature(self, regulatory: Regulatory) -> None:
        """A variant inside the regulatory feature fires the term."""
        assert "regulatory_region_variant" in _evaluate_consequence_terms(regulatory, start=1050, end=1050)

    @pytest.mark.parametrize("position", [_FEATURE_START - 1, _FEATURE_END + 1])
    def test_does_not_fire_outside(self, regulatory: Regulatory, position: int) -> None:
        """One base beyond either end is not an overlap."""
        assert "regulatory_region_variant" not in _evaluate_consequence_terms(
            regulatory, start=position, end=position
        )


class TestAblationTerms:
    """feature_ablation = contain_feature and deletion."""

    def test_tfbs_ablation_fires_for_a_structural_deletion(self, motif: Motif) -> None:
        """A structural deletion covering the motif ablates it."""
        assert "TFBS_ablation" in _spanning_deletion_terms(motif, structural=True)

    def test_tfbs_ablation_does_not_fire_for_an_ordinary_deletion(self, motif: Motif) -> None:
        """A plain (non-structural) deletion does not, under the catalogue's current gate."""
        assert "TFBS_ablation" not in _spanning_deletion_terms(motif, structural=False)

    def test_regulatory_region_ablation_fires_for_a_structural_deletion(self, regulatory: Regulatory) -> None:
        """Same for a regulatory feature."""
        assert "regulatory_region_ablation" in _spanning_deletion_terms(regulatory, structural=True)

    def test_ablation_does_not_fire_when_the_feature_is_not_fully_covered(self, motif: Motif) -> None:
        """contain_feature is a containment test: a partial deletion is not an ablation."""
        assert "TFBS_ablation" not in _evaluate_allele_consequence_terms(
            motif,
            start=_FEATURE_START,
            end=_FEATURE_START + 5,
            reference="AAAAAA",
            alternate="-",
            structural=True,
            class_so_term=ClassSOTerm.DELETION,
        )


class TestAmplificationTerms:
    """feature_amplification = contain_feature and copy_number_gain -- structural only."""

    def test_tfbs_amplification_fires_for_a_structural_gain(self, motif: Motif) -> None:
        """A structural allele lengthening a span that covers the motif amplifies it."""
        assert "TFBS_amplification" in _spanning_insertion_terms(motif, structural=True)

    def test_regulatory_region_amplification_fires_for_a_structural_gain(
        self, regulatory: Regulatory
    ) -> None:
        """Same for a regulatory feature."""
        assert "regulatory_region_amplification" in _spanning_insertion_terms(regulatory, structural=True)

    def test_amplification_does_not_fire_for_an_ordinary_insertion(self, motif: Motif) -> None:
        """The is_structural gate: a plain insertion is not a copy-number gain."""
        assert "TFBS_amplification" not in _spanning_insertion_terms(motif, structural=False)

    def test_ablation_and_amplification_are_mutually_exclusive(self, motif: Motif) -> None:
        """One shortens, the other lengthens; nothing does both."""
        terms = _spanning_insertion_terms(motif, structural=True)
        assert "TFBS_ablation" not in terms


class TestStructuralGateOnFeatureLevelTerms:
    """feature_elongation/feature_truncation carry the isa gate VEP writes as an early return."""

    def test_feature_elongation_needs_a_structural_allele(self, motif: Motif) -> None:
        """VEP bails out of feature_elongation for a TranscriptVariationAllele.

        Uses an insertion entirely inside the feature (_contained_insertion_terms), not one spanning
        across it (_spanning_insertion_terms): feature_elongation's ``when`` gates on within_feature,
        the variant lying entirely inside the feature, which is the opposite shape from
        contain_feature (the feature lying entirely inside the variant) that the amplification/
        ablation terms elsewhere in this file need.
        """
        assert "feature_elongation" in _contained_insertion_terms(motif, structural=True)
        assert "feature_elongation" not in _contained_insertion_terms(motif, structural=False)

    def test_feature_truncation_needs_a_structural_allele(self, motif: Motif) -> None:
        """VEP bails out of feature_truncation for a TranscriptVariationAllele."""
        assert "feature_truncation" in _spanning_deletion_terms(motif, structural=True)
        assert "feature_truncation" not in _spanning_deletion_terms(motif, structural=False)
