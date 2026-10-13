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
"""Tests for describing a corpus."""

from __future__ import annotations

import pytest

from ensembl.vep.corpus.describe import CorpusDescription, describe_variant, variant_class
from ensembl.vep.corpus.synthetic import SyntheticCorpus
from ensembl.vep.corpus.synthetic.default import default_corpus


@pytest.fixture(name="corpus", scope="module")
def fixture_corpus() -> SyntheticCorpus:
    """The default corpus, kept small."""
    return default_corpus(fuzz_count=5)


@pytest.fixture(name="description", scope="module")
def fixture_description(corpus: SyntheticCorpus) -> CorpusDescription:
    """The whole description, built once."""
    return CorpusDescription.build(corpus, batch_size=100)


class TestShape:
    """The chromosome, transcript and class tables."""

    def test_one_chromosome_with_the_corpus_totals(
        self, corpus: SyntheticCorpus, description: CorpusDescription
    ) -> None:
        """The default corpus is one chromosome carrying four transcripts."""
        (chromosome,) = description.chromosomes
        assert chromosome.name == "chrSynthetic"
        assert chromosome.transcripts == 4
        assert chromosome.overlaps == len(corpus.named_overlaps)

    def test_transcripts_carry_strand_and_coding_status(self, description: CorpusDescription) -> None:
        """The reverse-strand and non-coding transcripts are told apart."""
        by_id = {shape.stable_id: shape for shape in description.transcripts}
        assert by_id["ENST0002"].strand == "-"
        assert by_id["ENST0003"].coding is False
        assert by_id["ENST0004"].biotype == "nonsense_mediated_decay"
        assert sum(shape.overlaps for shape in description.transcripts) == description.chromosomes[0].overlaps

    def test_variant_classes_account_for_every_variant(self, description: CorpusDescription) -> None:
        """The class counts sum to the distinct variant count."""
        assert (
            sum(count for _kind, count in description.variant_classes) == description.chromosomes[0].variants
        )
        assert {kind for kind, _count in description.variant_classes} == {"SNV", "deletion", "insertion"}

    def test_chromosome_totals_include_variants_without_overlaps(
        self, corpus: SyntheticCorpus, description: CorpusDescription
    ) -> None:
        """The description counts every generated variant, not only those in the overlap listing."""
        overlapping_names = {name for name, _overlap in corpus.named_overlaps}
        assert any(name not in overlapping_names for name, _allele in corpus.named_variants)
        assert description.chromosomes[0].variants == len(corpus.named_variants)

    def test_transcripts_per_chunk_is_bounded_by_the_transcript_count(
        self, description: CorpusDescription
    ) -> None:
        """A chunk cannot reach more distinct transcripts than the corpus has."""
        assert 1.0 <= description.transcripts_per_chunk <= 4.0
        assert description.batch_size == 100

    def test_table_rows_follow_the_description_shapes(self, description: CorpusDescription) -> None:
        """Each row adapter preserves the values held by its description shape."""
        chromosome = description.chromosomes[0]
        assert description.chromosome_rows()[0] == (
            chromosome.name,
            chromosome.length,
            chromosome.transcripts,
            chromosome.variants,
            chromosome.overlaps,
        )

        transcript = description.transcripts[0]
        assert description.transcript_rows()[0] == (
            transcript.stable_id,
            transcript.chromosome,
            transcript.strand,
            transcript.start,
            transcript.end,
            transcript.exons,
            "yes" if transcript.coding else "no",
            transcript.biotype,
            transcript.overlaps,
        )

        overlap = description.overlaps[0]
        assert description.class_rows() == list(description.variant_classes)
        assert description.overlap_rows()[0] == (
            overlap.index,
            overlap.name,
            overlap.feature,
            overlap.locus,
            overlap.alleles,
        )


class TestOverlapListing:
    """The listing that gives ``vep explain`` its selectors."""

    def test_entries_follow_corpus_order(
        self, corpus: SyntheticCorpus, description: CorpusDescription
    ) -> None:
        """The index printed is the one --index accepts."""
        names = [name for name, _overlap in corpus.named_overlaps]
        assert [entry.name for entry in description.overlaps] == names
        assert [entry.index for entry in description.overlaps] == list(range(len(names)))

    def test_select_keeps_one_variants_overlaps(self, corpus: SyntheticCorpus) -> None:
        """Narrowing by name keeps every overlap of that variant and nothing else."""
        description = CorpusDescription.build(corpus, select="ENST0001_cds3_del1")
        assert {entry.name for entry in description.overlaps} == {"ENST0001_cds3_del1"}

    def test_a_long_allele_is_abbreviated_with_its_length(self, corpus: SyntheticCorpus) -> None:
        """A whole-transcript deletion stays a table cell."""
        description = CorpusDescription.build(corpus, select="ENST0001_ablation")
        (entry,) = description.overlaps
        assert "bp)" in entry.alleles
        assert len(entry.alleles) < 40


class TestVariantSpelling:
    """How a variant is classed and printed."""

    def test_classes_come_from_the_trimmed_alleles(self, corpus: SyntheticCorpus) -> None:
        """Empty alternate is a deletion, empty reference an insertion, one-for-one an SNV."""
        by_name = {name: overlap.overlapping_feature for name, overlap in corpus.named_overlaps}
        assert variant_class(by_name["ENST0001_cds3_del1"]) == "deletion"
        assert variant_class(by_name["ENST0001_elongation"]) == "insertion"
        assert variant_class(by_name["ENST0001_intron0_donor+1"]) == "SNV"

    def test_an_empty_allele_is_printed_as_a_dash(self, corpus: SyntheticCorpus) -> None:
        """The dash is what the --variant selector reads back as empty."""
        by_name = {name: overlap.overlapping_feature for name, overlap in corpus.named_overlaps}
        _locus, alleles = describe_variant(by_name["ENST0001_cds3_del1"])
        assert alleles.endswith(">-")
