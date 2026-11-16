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
"""Tests for writing a corpus to disk."""

from __future__ import annotations

from pathlib import Path

import pytest

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ keeps the annotators-before-caches import order.
from ensembl.vep.common import CoordSystemLevel, Strand
from ensembl.vep.corpus.overlaps import SyntheticOverlapSource
from ensembl.vep.corpus.store import (
    CORPUS_FILES,
    FASTA_NAME,
    GFF3_NAME,
    MANIFEST_NAME,
    VCF_NAME,
    CorpusManifest,
    distinct_variants,
    save_corpus,
    vcf_alleles,
)
from ensembl.vep.models.sequence import CoordSystem, Slice
from ensembl.vep.models.variation import VariationAllele
from ensembl.vep.parsers.fasta import FastaParser
from ensembl.vep.parsers.gff3 import GFF3Parser
from ensembl.vep.parsers.vcf import VCFParser
from ensembl.vep.sources.file.fasta.sequence import FastaSequenceSource

READER_FEATURE_TYPES = {"gene", "transcript", "exon", "cds", "five_prime_utr", "three_prime_utr"}
"""Lower-cased feature types the GFF3 transcript reader recognises or ignores harmlessly."""


@pytest.fixture(name="source", scope="module")
def fixture_source() -> SyntheticOverlapSource:
    """Both corpora, kept small."""
    return SyntheticOverlapSource(kind="both", fuzz=5, step=600)


@pytest.fixture(name="saved", scope="module")
def fixture_saved(
    source: SyntheticOverlapSource, tmp_path_factory: pytest.TempPathFactory
) -> tuple[Path, CorpusManifest]:
    """The corpus written once, and what the writer said about it."""
    directory = tmp_path_factory.mktemp("corpus")
    return directory, save_corpus(source.corpus, directory)


def _retrim(position: int, reference: str, alternate: str) -> tuple[int, str, str]:
    """Undo VCF padding: the inverse of vcf_alleles, written independently of it."""
    if reference and alternate and len(reference) == 1 and len(alternate) == 1:
        return position, reference, alternate
    if reference[0] == alternate[0]:
        return position + 1, reference[1:], alternate[1:]
    if reference[-1] == alternate[-1]:
        return position, reference[:-1], alternate[:-1]
    return position, reference, alternate


def _slice_with(sequence: str) -> Slice:
    """Return a forward slice carrying a known sequence."""
    built = Slice(
        coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
        seq_region_name="chrTest",
        start=1,
        end=len(sequence),
        strand=Strand.FORWARD,
    )
    with built.unsafe() as unsafe_slice:
        setattr(unsafe_slice, "_seq", sequence)
    return built


class TestFasta:
    """The chromosome sequences, back through the real parser and the real source."""

    def test_one_record_per_chromosome_with_the_whole_sequence(
        self, source: SyntheticOverlapSource, saved: tuple[Path, CorpusManifest]
    ) -> None:
        """Every corpus chromosome is written in full."""
        directory, _manifest = saved
        records = {record.seqid: record.sequence for record in FastaParser.parse_file(directory / FASTA_NAME)}
        for corpus in source.corpora:
            chromosome = corpus.synthetic_slice
            assert records[chromosome.seq_region_name] == chromosome.seq

    def test_the_sequence_source_reads_it_back(
        self, source: SyntheticOverlapSource, saved: tuple[Path, CorpusManifest]
    ) -> None:
        """A Slice over the saved file yields the sequence the corpus was built from."""
        directory, _manifest = saved
        original = source.corpora[0].synthetic_slice
        read_back = Slice(
            coord_system=CoordSystem(name="chromosome", level=CoordSystemLevel.CHROMOSOME),
            seq_region_name=original.seq_region_name,
            start=1,
            end=original.end,
            strand=Strand.FORWARD,
            sequence_source=FastaSequenceSource(file_path=directory / FASTA_NAME),
        )
        assert read_back.seq == original.seq


@pytest.fixture(name="gff3_records", scope="module")
def fixture_gff3_records(saved: tuple[Path, CorpusManifest]) -> list:
    """Every GFF3 record, parsed once."""
    directory, _manifest = saved
    return GFF3Parser.parse_file(directory / GFF3_NAME, limit=None)


@pytest.fixture(name="vcf_records", scope="module")
def fixture_vcf_records(saved: tuple[Path, CorpusManifest]) -> list:
    """Every VCF record, parsed once."""
    directory, _manifest = saved
    return VCFParser.parse_file(directory / VCF_NAME, limit=None)


class TestGff3:
    """The transcripts, back through the real parser, in the shape the readers expect."""

    def test_every_feature_type_is_one_the_readers_know(self, gff3_records: list) -> None:
        """The reader lower-cases types and keys on gene/transcript/exon/cds; nothing else appears."""
        assert {record.feature_type.lower() for record in gff3_records} <= READER_FEATURE_TYPES

    def test_one_transcript_row_per_transcript_joined_by_id(
        self, source: SyntheticOverlapSource, gff3_records: list
    ) -> None:
        """Every transcript is written once, under the ID the reader joins children to."""
        rows = {
            record.attributes["ID"][0]: record
            for record in gff3_records
            if record.feature_type == "transcript"
        }
        for corpus in source.corpora:
            for transcript in corpus.transcripts:
                row = rows[f"transcript:{transcript.stable_id}"]
                assert (row.start, row.end) == (transcript.start, transcript.end)
                assert row.strand_symbol == ("+" if transcript.strand == Strand.FORWARD else "-")

    def test_exon_rows_reproduce_the_spliced_exons(
        self, source: SyntheticOverlapSource, gff3_records: list
    ) -> None:
        """Exon rows carry the exon spans and point at their transcript."""
        spans_by_parent: dict[str, set[tuple[int, int]]] = {}
        for record in gff3_records:
            if record.feature_type == "exon":
                parent = record.attributes["Parent"][0]
                spans_by_parent.setdefault(parent, set()).add((record.start, record.end))
        for corpus in source.corpora:
            for transcript in corpus.transcripts:
                expected = {(exon.start, exon.end) for exon in transcript.spliced_exons}
                assert spans_by_parent[f"transcript:{transcript.stable_id}"] == expected

    def test_cds_rows_reproduce_the_cds_spans_with_phases(
        self, source: SyntheticOverlapSource, gff3_records: list
    ) -> None:
        """CDS rows carry the coding blocks, and the first block in reading order starts in phase 0."""
        rows_by_parent: dict[str, list] = {}
        for record in gff3_records:
            if record.feature_type == "CDS":
                rows_by_parent.setdefault(record.attributes["Parent"][0], []).append(record)
        for corpus in source.corpora:
            for transcript in corpus.transcripts:
                rows = rows_by_parent.get(f"transcript:{transcript.stable_id}", [])
                spans = transcript.cds_genomic_spans or []
                assert {(row.start, row.end) for row in rows} == set(spans)
                if spans:
                    first_start = spans[0][0]
                    first_row = next(row for row in rows if row.start == first_start)
                    assert first_row.phase_value == 0

    def test_a_transcript_without_a_gene_has_no_parent(self, gff3_records: list) -> None:
        """No gene is invented for a parentless transcript; the reader accepts the absence."""
        transcript_rows = {
            record.attributes["ID"][0]: record
            for record in gff3_records
            if record.feature_type == "transcript"
        }
        assert "Parent" not in transcript_rows["transcript:ENST0001"].attributes
        assert transcript_rows["transcript:ENST0002"].attributes["Parent"] == ("gene:ENSG0002",)

    def test_a_gene_is_written_once_with_its_biotype(self, gff3_records: list) -> None:
        """A gene shared by several transcripts appears once."""
        genes = [record for record in gff3_records if record.feature_type == "gene"]
        by_id = {record.attributes["ID"][0]: record for record in genes}
        assert len(genes) == len(by_id)
        assert by_id["gene:ENSG0004"].attributes["biotype"] == ("nonsense_mediated_decay",)


class TestVcf:
    """The variants, back through the real parser, named and re-trimmable."""

    def test_one_record_per_distinct_variant_named_in_the_id_column(
        self, source: SyntheticOverlapSource, vcf_records: list
    ) -> None:
        """A variant pairing with two transcripts is still one VCF record, under its name."""
        expected = {name for corpus in source.corpora for name, _allele in distinct_variants(corpus)}
        written = [record.id[0] for record in vcf_records]
        assert len(written) == len(expected)
        assert set(written) == expected

    def test_a_variant_without_overlaps_is_written(
        self, source: SyntheticOverlapSource, vcf_records: list
    ) -> None:
        """VCF export is driven by the corpus variants, not only its transcript overlaps."""
        corpus = source.corpora[0]
        overlapping_names = {name for name, _overlap in corpus.named_overlaps}
        zero_overlap_name = next(
            name for name, _allele in corpus.named_variants if name not in overlapping_names
        )
        assert zero_overlap_name in {record.id[0] for record in vcf_records}

    def test_every_record_retrims_to_the_model_alleles(
        self, source: SyntheticOverlapSource, vcf_records: list
    ) -> None:
        """Undoing the VCF padding gives back exactly the trimmed alleles the corpus holds."""
        alleles_by_name = {
            name: allele for corpus in source.corpora for name, allele in distinct_variants(corpus)
        }
        for record in vcf_records:
            allele = alleles_by_name[record.id[0]]
            assert _retrim(record.pos, record.reference, record.alt[0]) == (
                allele.start,
                allele.reference_allele or "",
                allele.alternate_allele or "",
            )

    def test_records_are_sorted_by_position_within_a_chromosome(self, vcf_records: list) -> None:
        """Position order is what an indexer, and a reader, expects."""
        for chromosome in {record.chrom for record in vcf_records}:
            positions = [record.pos for record in vcf_records if record.chrom == chromosome]
            assert positions == sorted(positions)


class TestVcfAlleles:
    """The padding rule, case by case."""

    def test_an_snv_needs_no_padding(self) -> None:
        """A substitution is written at its own position."""
        chromosome = _slice_with("ACGTACGT")
        allele = VariationAllele(
            start=3,
            end=3,
            strand=Strand.FORWARD,
            slice=chromosome,
            reference_allele="G",
            alternate_allele="T",
        )
        assert vcf_alleles(allele, "ACGTACGT") == (3, "G", "T")

    def test_a_deletion_is_anchored_on_the_base_before(self) -> None:
        """The reference gains the anchor and the alternate becomes it."""
        chromosome = _slice_with("ACGTACGT")
        allele = VariationAllele(
            start=3,
            end=4,
            strand=Strand.FORWARD,
            slice=chromosome,
            reference_allele="GT",
            alternate_allele="",
        )
        assert vcf_alleles(allele, "ACGTACGT") == (2, "CGT", "C")

    def test_an_insertion_is_anchored_on_the_base_before_its_position(self) -> None:
        """An insertion before position 3 sits after base 2, which is the anchor."""
        chromosome = _slice_with("ACGTACGT")
        allele = VariationAllele(
            start=3,
            end=2,
            strand=Strand.FORWARD,
            slice=chromosome,
            reference_allele="",
            alternate_allele="GG",
        )
        assert vcf_alleles(allele, "ACGTACGT") == (2, "C", "CGG")

    def test_an_indel_at_position_one_is_padded_with_the_base_after(self) -> None:
        """There is no base before position 1, so the spec's rule is to pad on the right."""
        chromosome = _slice_with("ACGTACGT")
        deletion = VariationAllele(
            start=1, end=1, strand=Strand.FORWARD, slice=chromosome, reference_allele="A", alternate_allele=""
        )
        insertion = VariationAllele(
            start=1,
            end=0,
            strand=Strand.FORWARD,
            slice=chromosome,
            reference_allele="",
            alternate_allele="TT",
        )
        assert vcf_alleles(deletion, "ACGTACGT") == (1, "AC", "C")
        assert vcf_alleles(insertion, "ACGTACGT") == (1, "A", "TTA")


class TestManifestAndSafety:
    """What the manifest records, and when a save is refused."""

    def test_the_manifest_counts_what_was_written(
        self, source: SyntheticOverlapSource, saved: tuple[Path, CorpusManifest]
    ) -> None:
        """The counts are the corpus's own, and the files listed exist."""
        directory, manifest = saved
        assert manifest.kind == "both"
        assert manifest.chromosomes == ("chrSynthetic", "chrHeavyLocus")
        assert manifest.transcripts == sum(len(corpus.transcripts) for corpus in source.corpora)
        assert manifest.overlaps == sum(len(corpus.named_overlaps) for corpus in source.corpora)
        assert set(manifest.seeds) == {"default", "heavy"}
        assert manifest.parameters == source.spec.parameter_hash()
        for file_name in CORPUS_FILES:
            assert (directory / file_name).exists()
        assert MANIFEST_NAME in (directory / MANIFEST_NAME).name

    def test_saving_twice_is_byte_identical(self, source: SyntheticOverlapSource, tmp_path: Path) -> None:
        """A corpus is a pure function of its spec, and so is its identity."""
        first = save_corpus(source.corpus, tmp_path / "one")
        second = save_corpus(source.corpus, tmp_path / "two")
        assert first.identity == second.identity
        for file_name in CORPUS_FILES:
            assert (tmp_path / "one" / file_name).read_bytes() == (tmp_path / "two" / file_name).read_bytes()

    def test_an_existing_corpus_is_not_overwritten_without_force(
        self, source: SyntheticOverlapSource, tmp_path: Path
    ) -> None:
        """The refusal names the files it found."""
        save_corpus(source.corpus, tmp_path)
        with pytest.raises(FileExistsError, match=VCF_NAME):
            save_corpus(source.corpus, tmp_path)
        assert save_corpus(source.corpus, tmp_path, force=True).identity
