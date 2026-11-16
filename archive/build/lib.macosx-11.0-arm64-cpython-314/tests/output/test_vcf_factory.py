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
"""Tests for the VCF output factory, over annotations of the synthetic corpus."""

from __future__ import annotations

from pathlib import Path

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the annotators and the facts the catalogue compiles against.
from ensembl.vep.annotators.factory import build_annotator
from ensembl.vep.common import Strand
from ensembl.vep.corpus.store import vcf_alleles
from ensembl.vep.corpus.synthetic import SyntheticCorpus, default_corpus, heavy_locus_corpus
from ensembl.vep.models.overlap import TranscriptVariationAllele
from ensembl.vep.models.variation import Variation
from ensembl.vep.output.factory import VCF_COLS, VCFOutputFactory, format_coords
from ensembl.vep.output.parser import VCFOutputParser

TIME = "2026-01-01 00:00:00"


@pytest.fixture(name="corpus", scope="module")
def fixture_corpus() -> SyntheticCorpus:
    """The boundary corpus, no fuzz."""
    return default_corpus(fuzz_count=0)


@pytest.fixture(name="heavy", scope="module")
def fixture_heavy() -> SyntheticCorpus:
    """A sparse heavy locus: variants that each reach many exon-sharing transcripts."""
    return heavy_locus_corpus(step=600)


def _annotated(corpus: SyntheticCorpus, name: str, annotator_name: str = "consequence.linear") -> Variation:
    """Build a Variation for one named corpus variant, with every overlap annotated and attached.

    Consequences are assigned to each overlap before the overlaps are attached to the allele,
    because attaching invalidates the allele's aggregate.
    """
    overlaps = [overlap for overlap_name, overlap in corpus.named_overlaps if overlap_name == name]
    source_allele = overlaps[0].overlapping_feature
    variation = Variation(
        start=source_allele.start,
        end=source_allele.end,
        strand=Strand.FORWARD,
        slice=source_allele.slice,
        stable_id=name,
        allele_string=f"{source_allele.reference_allele or ''}/{source_allele.alternate_allele or ''}",
    )
    allele = variation.variation_alleles[0]
    annotator = build_annotator(annotator=annotator_name)
    attached = []
    for overlap in overlaps:
        rebuilt = TranscriptVariationAllele(feature=overlap.feature, overlapping_feature=allele)
        rebuilt.overlap_consequences = annotator.annotate_feature_overlap(rebuilt)
        attached.append(rebuilt)
    allele.transcript_variation_alleles = attached
    return variation


def _csq_chunks(line: str, factory: VCFOutputFactory) -> list[dict[str, str]]:
    """Split a line's CSQ into one field dict per chunk."""
    info = line.split("\t")[7]
    csq = next(item for item in info.split(";") if item.startswith(f"{factory.info_field}="))
    chunks = []
    for chunk in csq.split("=", 1)[1].split(","):
        chunks.append(dict(zip(factory.fields(), chunk.split("|"), strict=True)))
    return chunks


def _only_chunk(line: str, factory: VCFOutputFactory) -> dict[str, str]:
    """Return a line's single CSQ chunk, failing if there is not exactly one."""
    chunks = _csq_chunks(line, factory)
    assert len(chunks) == 1, f"expected one chunk, got {len(chunks)}"
    return chunks[0]


class TestFormatCoords:
    """VEP's coordinate spelling."""

    def test_every_case(self) -> None:
        """Single, range, reversed, and missing bounds."""
        assert format_coords(5, 5) == "5"
        assert format_coords(5, 9) == "5-9"
        assert format_coords(9, 8) == "8-9"
        assert format_coords(5, None) == "5-?"
        assert format_coords(None, 9) == "?-9"
        assert format_coords(None, None) == "-"


class TestHeaders:
    """The header block."""

    def test_without_input_headers(self) -> None:
        """A fileformat line first, VEP's own lines, the column heading last."""
        headers = VCFOutputFactory(time=TIME, command_line="vep run --input_file x.vcf").headers()
        assert headers[0] == "##fileformat=VCFv4.1"
        assert headers[1].startswith('##VEP="v') and f'time="{TIME}"' in headers[1]
        assert headers[2] == (
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
            f'Format: {"|".join(VCF_COLS)}">'
        )
        assert headers[3] == "##VEP-command-line='vep run --input_file x.vcf'"
        assert headers[-1].startswith("#CHROM\tPOS\tID\tREF\tALT")

    def test_input_headers_are_carried_and_an_old_csq_is_replaced(self) -> None:
        """An earlier CSQ declaration and ##VEP line go; everything else stays, heading last."""
        input_headers = [
            "##fileformat=VCFv4.2",
            "##contig=<ID=1,length=100>",
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="old">',
            '##VEP="v110"',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
        ]
        headers = VCFOutputFactory(time=TIME).headers(input_headers)
        assert headers[0] == "##fileformat=VCFv4.2"
        assert "##contig=<ID=1,length=100>" in headers
        assert not any('Description="old"' in line for line in headers)
        assert not any(line == '##VEP="v110"' for line in headers)
        assert headers[-1].endswith("FORMAT\tSAMPLE")

    def test_keep_csq_keeps_the_old_declaration(self) -> None:
        """Nothing is dropped when asked to keep it."""
        input_headers = ['##INFO=<ID=CSQ,Number=.,Type=String,Description="old">', "#CHROM\tPOS"]
        headers = VCFOutputFactory(time=TIME, keep_csq=True).headers(input_headers)
        assert any('Description="old"' in line for line in headers)

    def test_flags_extend_the_declared_format(self) -> None:
        """The column layout follows VEP's for the same flags."""
        factory = VCFOutputFactory(variant_class=True, symbol=True, allele_number=True)
        assert factory.fields()[: len(VCF_COLS)] == list(VCF_COLS)
        assert factory.fields()[len(VCF_COLS) :] == [
            "ALLELE_NUM",
            "VARIANT_CLASS",
            "SYMBOL_SOURCE",
            "HGNC_ID",
        ]


class TestInfoChunk:
    """VEP's escaping inside a chunk."""

    def test_reserved_characters_are_replaced(self) -> None:
        """Commas, semicolons, whitespace and pipes can never split a chunk."""
        factory = VCFOutputFactory()
        chunk = factory.info_chunk({"Allele": "A", "Consequence": ["a,b", "c;d"], "Gene": "x y|z"})
        fields = dict(zip(factory.fields(), chunk.split("|"), strict=True))
        assert fields["Consequence"] == "a&b&c%3Bd"
        assert fields["Gene"] == "x_y&z"

    def test_a_dash_is_empty_except_as_the_allele(self) -> None:
        """A deletion's allele stays visible; any other field's dash means nothing."""
        factory = VCFOutputFactory()
        fields = dict(
            zip(factory.fields(), factory.info_chunk({"Allele": "-", "Codons": "-"}).split("|"), strict=True)
        )
        assert fields["Allele"] == "-"
        assert fields["Codons"] == ""

    def test_zero_is_kept(self) -> None:
        """A zero is a value, not an absence."""
        factory = VCFOutputFactory()
        fields = dict(zip(factory.fields(), factory.info_chunk({"EXON": 0}).split("|"), strict=True))
        assert fields["EXON"] == "0"

    def test_the_allele_is_reverse_complemented_on_the_reverse_strand(self) -> None:
        """Only the allele, and only when the variant is on the reverse strand."""
        factory = VCFOutputFactory()
        assert factory.info_chunk({"Allele": "AC"}, strand=-1).split("|", maxsplit=1)[0] == "GT"


class TestLinesFromTheModel:
    """Lines built from a variation with no input record."""

    def test_an_snv_carries_its_terms_and_transcript_fields(self, corpus: SyntheticCorpus) -> None:
        """Consequence, impact and feature come off the annotated overlap."""
        variation = _annotated(corpus, "ENST0001_intron0_donor+1")
        factory = VCFOutputFactory()
        line = factory.line(variation)
        columns = line.split("\t")
        assert columns[:3] == ["chrSynthetic", str(variation.start), "ENST0001_intron0_donor+1"]
        (chunk,) = [c for c in _csq_chunks(line, factory) if c["Feature"] == "ENST0001"]
        assert chunk["Consequence"].split("&")[0] == "splice_donor_variant"
        assert chunk["IMPACT"] == "HIGH"
        assert chunk["Feature_type"] == "Transcript"

    def test_a_deletion_record_matches_the_corpus_writers_padding(self, corpus: SyntheticCorpus) -> None:
        """One padding rule for both writers."""
        variation = _annotated(corpus, "ENST0001_cds3_del3")
        overlap = next(o for n, o in corpus.named_overlaps if n == "ENST0001_cds3_del3")
        chromosome = overlap.overlapping_feature.slice
        assert chromosome is not None and chromosome.seq is not None
        sequence = chromosome.seq
        position, reference, alternate = vcf_alleles(overlap.overlapping_feature, sequence)
        columns = VCFOutputFactory().line(variation).split("\t")
        assert columns[1:2] + columns[3:5] == [str(position), reference, alternate]

    def test_a_coding_overlap_writes_codons_and_amino_acids(self, corpus: SyntheticCorpus) -> None:
        """What the model computes is written; what it does not stays empty."""
        variation = _annotated(corpus, "ENST0001_start_codon_pos1050")
        factory = VCFOutputFactory()
        (chunk,) = [c for c in _csq_chunks(factory.line(variation), factory) if c["Feature"] == "ENST0001"]
        assert "/" in chunk["Codons"]
        assert chunk["Amino_acids"]
        assert chunk["HGVSc"] == ""

    def test_it_round_trips_through_the_parser(self, corpus: SyntheticCorpus, tmp_path: Path) -> None:
        """Written CSQ reads back as the same terms on the same allele."""
        variation = _annotated(corpus, "ENST0001_cds3_del1")
        factory = VCFOutputFactory(time=TIME)
        path = tmp_path / "out.vcf"
        path.write_text("\n".join([*factory.headers(), factory.line(variation)]) + "\n", encoding="utf-8")
        (parsed,) = list(VCFOutputParser().parse(path))
        expected = sorted(c.term for c in variation.variation_alleles[0].overlap_consequences)
        assert sorted(c.term for c in parsed.variation_alleles[0].overlap_consequences) == expected


class TestLinesFromARecord:
    """Lines that write an input record back."""

    RECORD = ["chrSynthetic", "1201", "rs1", "A", "G", "50", "PASS", "DP=10;CSQ=old|stuff;AF=0.5"]

    def test_the_input_columns_are_kept_and_the_old_csq_replaced(self, corpus: SyntheticCorpus) -> None:
        """QUAL, FILTER and the other INFO keys survive; the previous CSQ does not."""
        variation = _annotated(corpus, "ENST0001_intron0_donor+1")
        columns = VCFOutputFactory().line(variation, self.RECORD).split("\t")
        assert columns[:7] == self.RECORD[:7]
        assert columns[7].startswith("DP=10;AF=0.5;CSQ=")
        assert "old|stuff" not in columns[7]

    def test_keep_csq_leaves_the_old_one(self, corpus: SyntheticCorpus) -> None:
        """Asked to keep it, it stays and the new one is appended."""
        variation = _annotated(corpus, "ENST0001_intron0_donor+1")
        info = VCFOutputFactory(keep_csq=True).line(variation, self.RECORD).split("\t")[7]
        assert "CSQ=old|stuff" in info
        assert info.count("CSQ=") == 2

    def test_the_uploaded_name_falls_back_to_the_record(self, corpus: SyntheticCorpus) -> None:
        """Without a name of its own, the variant is named from its input columns."""
        variation = _annotated(corpus, "ENST0001_intron0_donor+1")
        unnamed = Variation(
            start=variation.start,
            end=variation.end,
            strand=Strand.FORWARD,
            slice=variation.slice,
            allele_string=variation.allele_string,
        )
        (uploaded,) = {
            h["Uploaded_variation"] for h in VCFOutputFactory().output_hashes(unnamed, self.RECORD)
        }
        assert uploaded == "chrSynthetic_1201_A/G"


class TestSelection:
    """Squashing, picking and intergenic alleles."""

    def test_most_severe_writes_one_chunk_with_the_worst_term(self, heavy: SyntheticCorpus) -> None:
        """One record per variation, carrying only its most severe term."""
        name = next(n for n, _ in heavy.named_overlaps)
        variation = _annotated(heavy, name)
        factory = VCFOutputFactory(most_severe=True)
        chunk = _only_chunk(factory.line(variation), factory)
        most_severe = variation.most_severe_consequence
        assert most_severe is not None
        assert chunk["Consequence"] == most_severe.term
        assert chunk["Feature"] == ""

    def test_summary_writes_every_distinct_term_most_severe_first(self, heavy: SyntheticCorpus) -> None:
        """All terms, once each, in rank order."""
        name = next(n for n, _ in heavy.named_overlaps)
        variation = _annotated(heavy, name)
        factory = VCFOutputFactory(summary=True)
        chunk = _only_chunk(factory.line(variation), factory)
        terms = chunk["Consequence"].split("&")
        assert len(terms) == len(set(terms))
        most_severe = variation.most_severe_consequence
        assert most_severe is not None
        assert terms[0] == most_severe.term

    def test_chunks_are_grouped_by_feature_and_terms_ordered_by_rank(self, heavy: SyntheticCorpus) -> None:
        """One chunk per feature, each chunk's terms most severe first."""
        name = next(n for n, _ in heavy.named_overlaps)
        variation = _annotated(heavy, name)
        factory = VCFOutputFactory()
        chunks = _csq_chunks(factory.line(variation), factory)
        overlaps = variation.variation_alleles[0].transcript_variation_alleles
        assert [c["Feature"] for c in chunks] == [o.feature.stable_id for o in overlaps]
        for chunk, overlap in zip(chunks, overlaps, strict=True):
            ranks = sorted(overlap.overlap_consequences, key=lambda c: c.rank)
            assert chunk["Consequence"].split("&") == [c.term for c in ranks]

    def test_pick_keeps_the_candidate_the_pick_order_prefers(self, heavy: SyntheticCorpus) -> None:
        """With no MANE, canonical, APPRIS or TSL data, biotype then rank then length decide."""
        name = next(
            n for n, _ in heavy.named_overlaps if len([1 for m, _ in heavy.named_overlaps if m == n]) > 3
        )
        variation = _annotated(heavy, name)
        factory = VCFOutputFactory(pick=True)
        chunk = _only_chunk(factory.line(variation), factory)

        def preference(overlap: TranscriptVariationAllele) -> tuple[int, int, int]:
            transcript = overlap.feature
            biotype = transcript.biotype or (transcript.gene.biotype if transcript.gene else None)
            rank = min(c.rank for c in overlap.overlap_consequences) if overlap.overlap_consequences else 1000
            length = transcript.cds_length if transcript.translation is not None else transcript.length
            return (0 if biotype == "protein_coding" else 1, rank, -(length or 0))

        overlaps = variation.variation_alleles[0].transcript_variation_alleles
        assert chunk["Feature"] == min(overlaps, key=preference).feature.stable_id

    def test_flag_pick_marks_exactly_one(self, heavy: SyntheticCorpus) -> None:
        """Everything is written, one carries PICK."""
        name = next(
            n for n, _ in heavy.named_overlaps if len([1 for m, _ in heavy.named_overlaps if m == n]) > 3
        )
        variation = _annotated(heavy, name)
        factory = VCFOutputFactory(flag_pick=True)
        chunks = _csq_chunks(factory.line(variation), factory)
        assert len(chunks) == len(variation.variation_alleles[0].transcript_variation_alleles)
        assert [c["PICK"] for c in chunks].count("1") == 1

    def test_an_allele_with_no_overlap_is_intergenic(self, corpus: SyntheticCorpus) -> None:
        """VEP's term for it, at the catalogue's impact; nothing at all under no_intergenic."""
        source = corpus.named_overlaps[0][1].overlapping_feature
        lonely = Variation(
            start=59000,
            end=59000,
            strand=Strand.FORWARD,
            slice=source.slice,
            stable_id="far",
            allele_string="A/T",
        )
        factory = VCFOutputFactory()
        chunk = _only_chunk(factory.line(lonely), factory)
        assert chunk["Consequence"] == "intergenic_variant"
        assert chunk["IMPACT"] == "MODIFIER"
        assert VCFOutputFactory(no_intergenic=True).line(lonely).split("\t")[7] == "."

    def test_flags_add_their_values(self, corpus: SyntheticCorpus) -> None:
        """Variant class and allele number appear under their own columns."""
        variation = _annotated(corpus, "ENST0001_intron0_donor+1")
        factory = VCFOutputFactory(variant_class=True, allele_number=True, biotype=True)
        (chunk,) = [c for c in _csq_chunks(factory.line(variation), factory) if c["Feature"] == "ENST0001"]
        assert chunk["VARIANT_CLASS"] == "SNV"
        assert chunk["ALLELE_NUM"] == "1"
        assert chunk["BIOTYPE"] == "protein_coding"
