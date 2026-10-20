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
"""Tests for annotating a corpus into variations."""

from __future__ import annotations

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the annotators and the facts they read.
from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.corpus.annotate import annotate_and_attach, annotated_variations
from ensembl.vep.corpus.base import Corpus
from ensembl.vep.corpus.registry import build_corpus
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.models.variation import Variation


@pytest.fixture(name="source", scope="module")
def fixture_source() -> Corpus:
    """The boundary corpus."""
    return build_corpus(CorpusSpec(kind="synthetic.default", fuzz=0))


@pytest.fixture(name="variations", scope="module")
def fixture_variations(source: Corpus) -> list[Variation]:
    """Every variant, annotated by linear."""
    return list(annotated_variations(source, build_annotator(annotator="consequence.linear")))


class TestShape:
    """One variation per variant, carrying every overlap."""

    def test_one_variation_per_distinct_variant_in_corpus_order(
        self, source: Corpus, variations: list[Variation]
    ) -> None:
        """Names appear once each, in the order the generator produced them."""
        names = [name for name, _allele in source.named_variants]
        assert [variation.stable_id for variation in variations] == names

    def test_every_overlap_is_attached_to_its_variant(
        self, source: Corpus, variations: list[Variation]
    ) -> None:
        """The features attached are the features the corpus paired the variant with."""
        expected: dict[str, list[str]] = {}
        for name, overlap in source.named_overlaps:
            expected.setdefault(name, []).append(overlap.feature.stable_id)
        for variation in variations:
            attached = variation.variation_alleles[0].transcript_variation_alleles
            assert [overlap.feature.stable_id for overlap in attached] == expected.get(
                str(variation.stable_id), []
            )

    def test_a_variant_without_overlaps_is_still_annotated(
        self, source: Corpus, variations: list[Variation]
    ) -> None:
        """Intergenic variants remain in the annotated corpus as empty-overlap variations."""
        overlapping_names = {name for name, _overlap in source.named_overlaps}
        zero_overlap_name = next(
            name for name, _allele in source.named_variants if name not in overlapping_names
        )
        variation = next(item for item in variations if item.stable_id == zero_overlap_name)
        allele = variation.variation_alleles[0]
        assert allele.transcript_variation_alleles == []
        assert allele.overlap_consequences == []


class TestTerms:
    """The rebuilt overlaps annotate as the corpus's own do."""

    def test_terms_match_the_corpus_overlaps(self, source: Corpus, variations: list[Variation]) -> None:
        """Rebuilding against a Variation's allele changes nothing the annotator sees."""
        annotator = build_annotator(annotator="consequence.linear")
        corpus_terms = {}
        for name, overlap in source.named_overlaps:
            terms = [consequence.term for consequence in annotator.annotate_feature_overlap(overlap)]
            corpus_terms[(name, overlap.feature.stable_id)] = terms
        for variation in variations:
            for overlap in variation.variation_alleles[0].transcript_variation_alleles:
                key = (str(variation.stable_id), overlap.feature.stable_id)
                assert [consequence.term for consequence in overlap.overlap_consequences] == corpus_terms[key]

    def test_consequences_reach_the_variation(self, variations: list[Variation]) -> None:
        """The aggregate the summary reads is populated, not a stale empty cache."""
        coding = next(variation for variation in variations if variation.stable_id == "ENST0001_cds3_del3")
        most_severe = coding.most_severe_consequence
        assert most_severe is not None
        assert most_severe.term == "inframe_deletion"


class TestAttach:
    """The helper on its own."""

    def test_an_allele_with_no_overlaps_gets_an_empty_list(self, variations: list[Variation]) -> None:
        """Which is what makes it intergenic downstream."""
        allele = variations[0].variation_alleles[0]
        annotate_and_attach(build_annotator(annotator="consequence.linear"), allele, [])
        assert allele.transcript_variation_alleles == []
        assert allele.overlap_consequences == []
