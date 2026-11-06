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
"""Tests for reading consequences back out of an annotated VCF."""

from __future__ import annotations

from pathlib import Path

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the facts the catalogue's rules compile against.
from ensembl.vep.models.overlap import OverlapConsequence
from ensembl.vep.models.variation import Variation
from ensembl.vep.output.parser import JSONOutputParser, TabOutputParser, VCFOutputParser, parser_for
from ensembl.vep.output.parser.vcf import UNKNOWN_RANK, annotation_layouts, csq_alleles

ANNOTATED_VCF = """##fileformat=VCFv4.2
##contig=<ID=1,length=248956422>
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE">
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	12345	rs1	A	T	.	.	CSQ=T|missense_variant&splice_region_variant|MODERATE|GENE|ENSG1|Transcript|ENST1|protein_coding,T|intron_variant|MODIFIER|GENE|ENSG1|Transcript|ENST2|protein_coding
1	12400	rs2	AC	A	.	.	CSQ=-|frameshift_variant|HIGH|GENE|ENSG1|Transcript|ENST1|protein_coding
1	12500	rs3	G	GTT,GT	.	.	CSQ=TT|inframe_insertion|MODERATE|GENE|ENSG1|Transcript|ENST1|protein_coding,T|frameshift_variant|HIGH|GENE|ENSG1|Transcript|ENST1|protein_coding
1	12600	rs4	G	A	.	.	ANN=A|synonymous_variant|LOW|GENE|ENSG1|transcript|ENST1|protein_coding|1/3|c.1G>A|p.X|1/3|1/3|1/1||
1	12700	rs5	G	A	.	.	CSQ=A|made_up_term|MODIFIER|GENE|ENSG1|Transcript|ENST1|protein_coding
1	12800	.	G	A	.	.	.
"""


@pytest.fixture(name="annotated", scope="module")
def fixture_annotated(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A small annotated VCF on disk."""
    path = tmp_path_factory.mktemp("output") / "annotated.vcf"
    path.write_text(ANNOTATED_VCF, encoding="utf-8")
    return path


@pytest.fixture(name="variations", scope="module")
def fixture_variations(annotated: Path) -> dict[str, Variation]:
    """Every record, keyed by ID (or position for the unnamed one)."""
    parsed = {}
    for variation in VCFOutputParser().parse(annotated):
        parsed[variation.stable_id or str(variation.start)] = variation
    return parsed


class TestRoundTripWithoutSources:
    """The claim that matters: the summary's inputs answer with zero transcripts attached."""

    def test_most_severe_and_class_answer_with_no_transcripts(self, variations: dict[str, Variation]) -> None:
        """Terms live on the allele, so the two headline charts read straight off Variation."""
        variation = variations["rs1"]
        assert variation.variation_alleles[0].transcript_variation_alleles == []
        assert variation.most_severe_consequence is not None
        assert variation.most_severe_consequence.term == "missense_variant"
        assert variation.class_so_term is not None
        assert variation.class_so_term.value == "SNV"

    def test_a_record_read_back_names_its_own_region(self, variations: dict[str, Variation]) -> None:
        """Nothing resolved this record's name, so the file's spelling is its uploaded locus."""
        variation = variations["rs1"]
        assert variation.uploaded_locus is not None
        assert variation.slice is not None
        assert variation.uploaded_locus.seq_region_name == variation.slice.seq_region_name

    def test_every_csq_entry_contributes_its_terms(self, variations: dict[str, Variation]) -> None:
        """Two entries for one allele, one carrying two terms joined by &, give three consequences."""
        terms = sorted(c.term for c in variations["rs1"].variation_alleles[0].overlap_consequences)
        assert terms == ["intron_variant", "missense_variant", "splice_region_variant"]

    def test_a_known_term_carries_the_catalogue_rank(self, variations: dict[str, Variation]) -> None:
        """Rank, accession and impact come from the rule the term is defined by."""
        (consequence,) = variations["rs2"].variation_alleles[0].overlap_consequences
        assert consequence.term == "frameshift_variant"
        assert consequence.accession.startswith("SO:")
        assert 0 < consequence.rank < UNKNOWN_RANK

    def test_an_unknown_term_is_kept_below_every_real_one(self, variations: dict[str, Variation]) -> None:
        """A term the catalogue does not model must not become the most severe by default."""
        (consequence,) = variations["rs5"].variation_alleles[0].overlap_consequences
        assert consequence.term == "made_up_term"
        assert consequence.rank == UNKNOWN_RANK

    def test_multiple_alternates_are_matched_by_the_csq_allele_spelling(
        self, variations: dict[str, Variation]
    ) -> None:
        """G>GTT,GT: VEP spells the alleles TT and T, and each entry lands on its own allele."""
        first, second = variations["rs3"].variation_alleles
        assert [c.term for c in first.overlap_consequences] == ["inframe_insertion"]
        assert [c.term for c in second.overlap_consequences] == ["frameshift_variant"]

    def test_ann_is_read_the_same_way(self, variations: dict[str, Variation]) -> None:
        """SnpEff's ANN declares its layout in quotes; the Annotation field is the term."""
        (consequence,) = variations["rs4"].variation_alleles[0].overlap_consequences
        assert consequence.term == "synonymous_variant"

    def test_an_unannotated_record_has_no_consequences(self, variations: dict[str, Variation]) -> None:
        """Nothing seeded, nothing invented."""
        assert variations["12800"].variation_alleles[0].overlap_consequences == []
        assert variations["12800"].most_severe_consequence is None


class TestSeedingOrder:
    """The trap the reader is built around, pinned rather than remembered."""

    def test_assigning_overlaps_after_seeding_wipes_the_seed(self) -> None:
        """transcript_variation_alleles invalidates the seed on assignment."""
        variation = Variation(start=1, end=1, strand=1, allele_string="A/T")  # type: ignore[arg-type]
        allele = variation.variation_alleles[0]
        allele.seed_attrs(overlap_consequences=[OverlapConsequence(term="missense_variant", accession="")])
        assert allele.overlap_consequences
        allele.transcript_variation_alleles = []
        assert allele.overlap_consequences == []


class TestLayoutsAndAlleles:
    """The header and allele conventions, on their own."""

    def test_layouts_come_from_the_header(self) -> None:
        """VEP's Format: and SnpEff's quoted list both yield field names."""
        metadata = {
            "info": {
                "CSQ": {
                    "description": "Consequence annotations from Ensembl VEP. Format: Allele|Consequence"
                },
                "ANN": {"description": "Functional annotations: 'Allele | Annotation | Impact'"},
            }
        }
        layouts = annotation_layouts(metadata)
        assert layouts["CSQ"] == ["Allele", "Consequence"]
        assert layouts["ANN"] == ["Allele", "Annotation", "Impact"]

    def test_csq_alleles_drop_shared_padding(self) -> None:
        """An indel loses its anchor base; an emptied allele becomes a dash; an SNV is untouched."""
        assert csq_alleles("AC", ["A"]) == ["-"]
        assert csq_alleles("G", ["GTT", "GT"]) == ["TT", "T"]
        assert csq_alleles("A", ["T"]) == ["T"]
        assert csq_alleles("A", ["T", "AG"]) == ["T", "AG"]


class TestParserFor:
    """Choosing a reader."""

    def test_by_extension_and_by_declared_format(self, annotated: Path) -> None:
        """The extension decides, unless the format is given."""
        assert isinstance(parser_for(annotated), VCFOutputParser)
        assert isinstance(parser_for("out.vcf.gz"), VCFOutputParser)
        assert isinstance(parser_for("out.txt"), TabOutputParser)
        assert isinstance(parser_for("anything", "json"), JSONOutputParser)

    def test_the_seams_say_they_are_not_implemented(self) -> None:
        """Declared, not hidden."""
        with pytest.raises(NotImplementedError):
            next(TabOutputParser().parse("out.tab"))
        with pytest.raises(NotImplementedError):
            next(JSONOutputParser().parse("out.json"))

    def test_an_unknown_format_is_rejected(self) -> None:
        """Named, so the caller can fix it."""
        with pytest.raises(ValueError, match="no reader"):
            parser_for("out.xyz")
