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
"""Tests for the overlap summary an explanation opens with."""

from __future__ import annotations

import pytest

from ensembl.vep.corpus import Corpus, build_corpus
from ensembl.vep.corpus.spec import CorpusSpec
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.explain.selection import FeatureSelector
from ensembl.vep.explain.summary import format_summary, summary_html, summary_rows


@pytest.fixture(name="source", scope="module")
def fixture_source() -> Corpus:
    """The smallest synthetic corpus, built once for the module."""
    return build_corpus(CorpusSpec(kind="synthetic.default", fuzz=0))


def _overlap(source: Corpus, name: str, feature: str) -> TranscriptVariationAllele:
    """Return the corpus overlap of one named variant over one transcript."""
    for overlap_name, overlap in source.named_overlaps:
        if overlap_name == name and FeatureSelector.parse(feature).matches(overlap.feature):
            return overlap
    raise LookupError(f"no overlap named {name!r} over {feature!r}")


def _rows(overlap: TranscriptVariationAllele) -> dict[str, str]:
    """Return the summary as a label -> value mapping."""
    return dict(summary_rows(overlap))


class TestCodingOverlap:
    """A frameshift in a protein-coding transcript: every row applies."""

    def test_the_change_is_described_from_codon_to_whole_peptide(self, source: Corpus) -> None:
        """Codons, residues and where the full CDS and peptide first part company."""
        rows = _rows(_overlap(source, "ENST0001_cds3_del1", "ENST0001"))

        assert rows["location"].startswith("exon ")
        assert rows["positions"].startswith("cDNA ")
        assert "CDS " in rows["positions"] and "protein " in rows["positions"]
        assert "/" in rows["codons"]
        assert "first difference at" in rows["full CDS"]
        assert "[" in rows["full CDS"] and "]" in rows["full CDS"]
        assert "first difference at" in rows["full peptide"]

    def test_the_variant_names_its_record_and_allele(self, source: Corpus) -> None:
        """The locus and alleles, then the allele's number among its record's alternates."""
        overlap = _overlap(source, "ENST0001_cds3_del1", "ENST0001")
        allele = overlap.overlapping_feature

        assert _rows(overlap)["variant"].startswith(f"chrSynthetic:{allele.start}")

    def test_the_text_block_is_aligned_under_a_heading(self, source: Corpus) -> None:
        """One ``summary:`` heading, then every label padded to the same width."""
        lines = format_summary(_overlap(source, "ENST0001_cds3_del1", "ENST0001")).splitlines()

        assert lines[0] == "summary:"
        value_columns = {len(line) - len(line[2:].split("  ", 1)[1].lstrip()) for line in lines[1:]}
        assert len(value_columns) == 1


class TestOtherOverlaps:
    """Rows that describe where a variant is when it changes no coding sequence."""

    def test_an_intronic_variant_changes_no_cds(self, source: Corpus) -> None:
        """Intronic: no codons, no positions, and the whole CDS unchanged."""
        rows = _rows(_overlap(source, "ENST0001_intron0_acceptor-17", "ENST0001"))

        assert rows["location"].startswith("intron 1/")
        assert "codons" not in rows
        assert "positions" not in rows
        assert rows["full CDS"].startswith("unchanged")

    def test_a_non_coding_transcript_says_so(self, source: Corpus) -> None:
        """No coding span, and nothing to translate."""
        rows = _rows(_overlap(source, "ENST0003_noncoding_exon", "ENST0003"))

        assert rows["coding"] == "non-coding"
        assert "codons" not in rows
        assert rows["positions"].startswith("cDNA ")

    def test_an_upstream_variant_gives_its_distance(self, source: Corpus) -> None:
        """Outside the transcript the location is a distance and a direction.

        The corpus clamps this variant to the chromosome's start, so the distance is read back off
        the overlap rather than taken from its name.
        """
        overlap = _overlap(source, "ENST0001_upstream5000", "ENST0001")
        distance = overlap.feature.start - overlap.overlapping_feature.end

        assert _rows(overlap)["location"] == f"{distance} bp upstream of the transcript"


class TestHtml:
    """The table the SVG explanation page opens with."""

    def test_every_row_is_an_escaped_table_row(self, source: Corpus) -> None:
        """Values carry ``->`` and ``[``; the table must hold them as text, not markup."""
        overlap = _overlap(source, "ENST0001_cds3_del1", "ENST0001")

        html = summary_html(overlap)

        assert html.startswith('<table class="vep-table">')
        assert html.count("<tr>") == len(summary_rows(overlap))
        assert "-&gt;" in html
