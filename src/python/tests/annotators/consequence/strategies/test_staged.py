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
"""Tests for the staged consequence annotator."""

from __future__ import annotations

import pytest

from ensembl.vep.annotators.consequence.facts.base import closure_cost
from ensembl.vep.annotators.consequence.rules.rulesets import build_fact_ruleset
from ensembl.vep.annotators.consequence.strategies.staged import (
    CheapAtomAlphabet,
    StagedFactsAnnotator,
)
from ensembl.vep.caches.codon_cache import CodonEditLRUCache
from ensembl.vep.caches.zone_map import ZoneMapCache
from ensembl.vep.corpus.synthetic import SyntheticCorpus
from ensembl.vep.models.overlap import TranscriptVariationAllele


class _NotATranscriptOverlap:
    """Object that is not a transcript variation allele."""


def _annotated_terms(
    annotator: StagedFactsAnnotator, overlaps: list[TranscriptVariationAllele]
) -> list[list[str]]:
    """Return ordered consequence terms for each overlap."""
    return [
        [consequence.term for consequence in annotator.annotate_feature_overlap(overlap)]
        for overlap in overlaps
    ]


class TestCheapAtomAlphabet:
    """Tests for selecting facts evaluated in the cheap stage."""

    def test_empty_rules_create_an_empty_alphabet(self) -> None:
        """Test an empty rule collection has no cheap atoms."""
        alphabet = CheapAtomAlphabet.build(())

        assert alphabet.atom_names == ()

    def test_atoms_obey_the_cost_threshold_and_order(self) -> None:
        """Test cheap atoms are filtered and sorted by closure cost."""
        ruleset = build_fact_ruleset()

        alphabet = CheapAtomAlphabet.build(ruleset.rules, max_closure_cost=1)
        costs = [closure_cost(atom_name) for atom_name in alphabet.atom_names]

        assert alphabet.atom_names
        assert all(cost <= 1 for cost in costs)
        assert costs == sorted(costs)


class TestEarlySubsumptionSkip:
    """Tests for the optional subsumption-based early skip."""

    @pytest.fixture
    def overlaps(self, small_default_corpus: SyntheticCorpus) -> list[TranscriptVariationAllele]:
        """Create transcript variation alleles from the default corpus."""
        return [overlap for _name, overlap in small_default_corpus.named_overlaps]

    def test_enabled_skip_preserves_terms(self, overlaps: list[TranscriptVariationAllele]) -> None:
        """Test enabling the early skip does not change annotations."""
        skipped_terms = _annotated_terms(StagedFactsAnnotator(early_subsumption_skip=True), overlaps)
        baseline_terms = _annotated_terms(StagedFactsAnnotator(early_subsumption_skip=False), overlaps)

        assert skipped_terms == baseline_terms

    def test_enabled_skip_avoids_rule_evaluations(self, overlaps: list[TranscriptVariationAllele]) -> None:
        """Test enabling the early skip suppresses at least one evaluation."""
        annotator = StagedFactsAnnotator(early_subsumption_skip=True)

        _annotated_terms(annotator, overlaps)

        assert annotator.skipped_by_subsumption > 0

    def test_disabled_skip_leaves_the_counter_at_zero(
        self, overlaps: list[TranscriptVariationAllele]
    ) -> None:
        """Test disabling the early skip does not suppress evaluations."""
        annotator = StagedFactsAnnotator(early_subsumption_skip=False)

        _annotated_terms(annotator, overlaps)

        assert annotator.skipped_by_subsumption == 0


class TestInputValidation:
    """Tests for staged annotator input validation."""

    def test_non_transcript_overlap_is_rejected(self) -> None:
        """Test a non-transcript overlap raises a type error."""
        annotator = StagedFactsAnnotator()

        with pytest.raises(TypeError, match="only annotates TranscriptVariationAllele"):
            annotator.annotate_feature_overlap(_NotATranscriptOverlap())  # type: ignore[arg-type]


class TestOptionalCaches:
    """Tests for the optional strategy caches."""

    def test_zone_map_is_used_without_changing_terms(
        self, boundary_overlaps: dict[str, TranscriptVariationAllele]
    ) -> None:
        """Test the enabled zone map serves facts without changing the result."""
        overlap = boundary_overlaps["ENST0001::ENST0001_cds3_del1"]
        expected_terms = [
            consequence.term for consequence in StagedFactsAnnotator().annotate_feature_overlap(overlap)
        ]
        zone_map = ZoneMapCache()
        annotator = StagedFactsAnnotator(zone_map=zone_map)

        consequences = annotator.annotate_feature_overlap(overlap)

        assert [consequence.term for consequence in consequences] == expected_terms
        assert annotator.annotator_context.zone_map_cache is zone_map
        assert zone_map.maps_built == 1
        assert zone_map.facts_served > 0

    def test_codon_cache_is_reused_without_changing_terms(
        self, boundary_overlaps: dict[str, TranscriptVariationAllele]
    ) -> None:
        """Test the enabled codon cache serves a repeated coding edit."""
        overlap = boundary_overlaps["ENST0001::ENST0001_cds3_del1"]
        codon_cache = CodonEditLRUCache()
        annotator = StagedFactsAnnotator(codon_cache=codon_cache)

        consequences = annotator.annotate_feature_overlaps([overlap, overlap])

        assert consequences[1] == consequences[0]
        assert annotator.annotator_context.codon_cache is codon_cache
        assert codon_cache.miss_count == 1
        assert codon_cache.hit_count == 1
