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
"""Conformance tests for registered annotators."""

from __future__ import annotations

import pytest

from ensembl.vep.annotators.base import Annotator
from ensembl.vep.annotators.consequence.facts.base import verify_requires
from ensembl.vep.annotators.consequence.context.fact import FactContext
from ensembl.vep.annotators.registry import available_annotators, get_annotator
from ensembl.vep.corpus.synthetic import SyntheticCorpus
from ensembl.vep.models.overlap import OverlapConsequence


def _consequence_terms(consequences: list[OverlapConsequence]) -> list[str]:
    """Return the rank-ordered terms emitted by an annotator."""
    return [consequence.term for consequence in consequences]


def _consequence_term_set(consequences: list[OverlapConsequence]) -> frozenset[str]:
    """Return the unordered terms emitted by an annotator."""
    return frozenset(consequence.term for consequence in consequences)


def _assert_corpus_conformance(corpus: SyntheticCorpus, *, corpus_label: str) -> None:
    """Check every registered strategy against ``consequence.linear`` over a corpus."""
    reference = get_annotator("consequence.linear")
    candidate_strategy_names = [name for name in available_annotators() if name != "consequence.linear"]
    assert candidate_strategy_names, "expected at least one non-reference strategy to be registered"
    candidate_annotators: dict[str, Annotator] = {
        name: get_annotator(name) for name in candidate_strategy_names
    }

    list_mismatches: list[str] = []
    set_mismatches: list[str] = []
    for overlap_name, overlap in corpus.named_overlaps:
        reference_consequences = reference.annotate_feature_overlap(overlap)
        reference_terms = _consequence_terms(reference_consequences)
        reference_term_set = _consequence_term_set(reference_consequences)
        for strategy_name, annotator in candidate_annotators.items():
            candidate_consequences = annotator.annotate_feature_overlap(overlap)
            candidate_terms = _consequence_terms(candidate_consequences)
            if candidate_terms != reference_terms:
                list_mismatches.append(
                    f"[{corpus_label}] {strategy_name} disagrees with consequence.linear on "
                    f"{overlap_name}: {candidate_terms} != {reference_terms}"
                )
            candidate_term_set = _consequence_term_set(candidate_consequences)
            if candidate_term_set != reference_term_set:
                set_mismatches.append(
                    f"[{corpus_label}] {strategy_name} term set disagrees with consequence.linear on "
                    f"{overlap_name}: {candidate_term_set} != {reference_term_set}"
                )

    assert not list_mismatches, "\n".join(list_mismatches[:20])
    assert not set_mismatches, "\n".join(set_mismatches[:20])


class TestBoundaryAndFuzzCorpusConformance:
    """The four-transcript boundary and fuzz corpus."""

    def test_every_registered_strategy_matches_the_reference(
        self, full_default_corpus: SyntheticCorpus
    ) -> None:
        """Every registered strategy matches the reference."""
        _assert_corpus_conformance(full_default_corpus, corpus_label="boundary+fuzz")

    def test_every_requires_claim_holds(self, full_default_corpus: SyntheticCorpus) -> None:
        """Every declared ``requires`` implication holds."""
        fact_contexts = [
            FactContext(overlap=overlap) for _name, overlap in full_default_corpus.named_overlaps
        ]
        problems = verify_requires(fact_contexts)
        assert not problems


class TestHeavyLocusCorpusConformance:
    """The heavy-locus corpus."""

    def test_every_registered_strategy_matches_the_reference(
        self, full_heavy_locus_corpus: SyntheticCorpus
    ) -> None:
        """Every registered strategy matches the reference."""
        _assert_corpus_conformance(full_heavy_locus_corpus, corpus_label="heavy-locus")

    def test_every_requires_claim_holds(self, full_heavy_locus_corpus: SyntheticCorpus) -> None:
        """Every declared ``requires`` implication holds."""
        fact_contexts = [
            FactContext(overlap=overlap) for _name, overlap in full_heavy_locus_corpus.named_overlaps
        ]
        problems = verify_requires(fact_contexts)
        assert not problems


class TestBatchPathConformance:
    """Every annotator must agree with itself."""

    @pytest.fixture(scope="class")
    def batch_annotator_names(self) -> list[str]:
        """Registry names of every currently-registered strategy."""
        return list(available_annotators())

    def test_every_registered_annotator_implements_the_batch_path(self) -> None:
        """Both halves of the contract are available on every registered strategy."""
        for name in available_annotators():
            annotator = get_annotator(name)
            assert callable(annotator.annotate_feature_overlaps)
            assert callable(annotator.annotate_feature_overlap)

    def test_batch_output_matches_single_overlap_output(
        self, batch_annotator_names: list[str], small_default_corpus: SyntheticCorpus
    ) -> None:
        """Batch output matches single-overlap output."""
        overlaps = [overlap for _name, overlap in small_default_corpus.named_overlaps]
        for strategy_name in batch_annotator_names:
            annotator = get_annotator(strategy_name)
            batch_terms = [
                _consequence_terms(consequences)
                for consequences in annotator.annotate_feature_overlaps(overlaps)
            ]
            single_terms = [
                _consequence_terms(annotator.annotate_feature_overlap(overlap)) for overlap in overlaps
            ]
            assert batch_terms == single_terms, f"{strategy_name}: batch output != single-overlap output"

    def test_batch_output_matches_the_reference(
        self, batch_annotator_names: list[str], small_default_corpus: SyntheticCorpus
    ) -> None:
        """Batch output matches the reference."""
        reference = get_annotator("consequence.linear")
        overlaps = [overlap for _name, overlap in small_default_corpus.named_overlaps]
        reference_terms = [
            _consequence_terms(reference.annotate_feature_overlap(overlap)) for overlap in overlaps
        ]
        for strategy_name in batch_annotator_names:
            annotator = get_annotator(strategy_name)
            batch_terms = [
                _consequence_terms(consequences)
                for consequences in annotator.annotate_feature_overlaps(overlaps)
            ]
            assert batch_terms == reference_terms, f"{strategy_name}: batch output != consequence.linear"
