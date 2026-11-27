# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Tests for constructing VCF golden-corpus input from a TSV table."""

# pylint: disable=import-error,redefined-outer-name

from __future__ import annotations

import json
from pathlib import Path

import pytest

from ensembl.vep.annotators.consequence.rules.catalogue import catalogue_digest
from scripts.generate_new_golden_corpus_vcf import (
    GoldenVcfInputError,
    REQUIRED_COLUMNS,
    _read_rows,
    write_vcf,
)


@pytest.fixture
def reference_files(tmp_path: Path) -> tuple[Path, Path]:
    """Create one tiny reference with forward and reverse coding transcripts."""

    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + "ACGT" * 20 + "\n", encoding="utf-8")
    gff3 = tmp_path / "annotation.gff3"
    gff3.write_text(
        "\n".join(
            (
                "##gff-version 3",
                "chr1\ttest\tgene\t10\t30\t.\t+\t.\tID=g1",
                "chr1\ttest\tmRNA\t10\t30\t.\t+\t.\tID=tx-forward;Parent=g1",
                "chr1\ttest\texon\t10\t15\t.\t+\t.\tID=e1;Parent=tx-forward",
                "chr1\ttest\tCDS\t12\t15\t.\t+\t0\tID=p1;Parent=tx-forward",
                "chr1\ttest\texon\t20\t30\t.\t+\t.\tID=e2;Parent=tx-forward",
                "chr1\ttest\tCDS\t20\t27\t.\t+\t0\tID=p1;Parent=tx-forward",
                "chr1\ttest\tgene\t40\t60\t.\t-\t.\tID=g2",
                "chr1\ttest\tmRNA\t40\t60\t.\t-\t.\tID=tx-reverse;Parent=g2",
                "chr1\ttest\texon\t40\t45\t.\t-\t.\tID=e3;Parent=tx-reverse",
                "chr1\ttest\tCDS\t43\t45\t.\t-\t0\tID=p2;Parent=tx-reverse",
                "chr1\ttest\texon\t50\t60\t.\t-\t.\tID=e4;Parent=tx-reverse",
                "chr1\ttest\tCDS\t50\t58\t.\t-\t0\tID=p2;Parent=tx-reverse",
                "",
            )
        ),
        encoding="utf-8",
    )
    return fasta, gff3


def _row(
    variant_id: str,
    transcript_id: str,
    start: str,
    end: str,
    reference: str,
    alternate: str,
    expected: str,
    context: str,
    *,
    transcript_start: int = 10,
    transcript_end: int = 30,
    strand: str = "1",
) -> str:
    """Return one well-formed TSV test-case row."""

    return "\t".join(
        (
            variant_id,
            transcript_id,
            "gene",
            "GENE",
            "chr1",
            str(transcript_start),
            str(transcript_end),
            strand,
            start,
            end,
            reference,
            alternate,
            expected,
            context,
            "coordinates",
            "test case",
            "ENST00000000000",
        )
    )


def _write_cases(path: Path, rows: list[str]) -> None:
    """Write a test-case TSV with the generator's exact header."""

    path.write_text("\t".join(REQUIRED_COLUMNS) + "\n" + "\n".join(rows) + "\n", encoding="utf-8")


def test_writes_vcf_for_transcript_coordinates_and_indels(
    tmp_path: Path, reference_files: tuple[Path, Path]
) -> None:
    """Expressions, NA alleles, insertion anchoring, and deletion anchoring are all rendered."""

    fasta, gff3 = reference_files
    cases = tmp_path / "cases.tsv"
    _write_cases(
        cases,
        [
            _row(
                "forward-reference",
                "tx-forward",
                "t_start",
                "t_start",
                "NA",
                "NA",
                "synonymous_variant",
                "forward ref",
            ),
            _row(
                "deletion",
                "tx-forward",
                "exon_1_start",
                "exon_1_start + 1",
                "NA",
                "-",
                "feature_truncation",
                "deletion case",
            ),
            _row(
                "insertion",
                "tx-forward",
                "exon_1_end + 1",
                "exon_1_end",
                "-",
                "A",
                "insertion",
                "insert case",
            ),
            _row(
                "reverse-coding",
                "tx-reverse",
                "cds_start",
                "cds_start",
                "NA",
                "T",
                "missense_variant",
                "reverse coding",
                transcript_start=40,
                transcript_end=60,
                strand="-1",
            ),
        ],
    )
    output_vcf = tmp_path / "corpus.vcf"
    output_json = tmp_path / "terms.json"

    assert write_vcf(cases, fasta, gff3, output_vcf, output_json) == 4

    assert output_vcf.read_text(encoding="utf-8").splitlines() == [
        "##fileformat=VCFv4.2",
        "##source=ensembl_vep_golden",
        "##contig=<ID=chr1,length=80>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t9\tdeletion\tACG\tA\t.\t.\t.",
        "chr1\t10\tforward-reference\tC\tG\t.\t.\t.",
        "chr1\t15\tinsertion\tG\tGA\t.\t.\t.",
        "chr1\t43\treverse-coding\tG\tT\t.\t.\t.",
    ]
    assert json.loads(output_json.read_text(encoding="utf-8")) == {
        "format_version": 1,
        "spec": {"kind": "golden.file", "parameters": {"dataset": "default"}},
        "annotator": "consequence.linear",
        "catalogue": catalogue_digest(),
        "terms": {
            "tx-forward::deletion": ["feature_truncation"],
            "tx-forward::forward-reference": ["synonymous_variant"],
            "tx-forward::insertion": ["insertion"],
            "tx-reverse::reverse-coding": ["missense_variant"],
        },
    }


def test_uses_transcript_order_for_reverse_exons_and_introns(
    tmp_path: Path, reference_files: tuple[Path, Path]
) -> None:
    """Reverse-strand exon and intron numbering comes from the Transcript object."""

    fasta, gff3 = reference_files
    cases = tmp_path / "cases.tsv"
    _write_cases(
        cases,
        [
            _row(
                "reverse-intron",
                "tx-reverse",
                "intron_1_start",
                "intron_1_end",
                "NA",
                "A",
                "coding_sequence_variant",
                "reverse exon intron",
                transcript_start=40,
                transcript_end=60,
                strand="-1",
            )
        ],
    )
    output_vcf = tmp_path / "corpus.vcf"
    output_json = tmp_path / "terms.json"

    write_vcf(cases, fasta, gff3, output_vcf, output_json)

    assert output_vcf.read_text(encoding="utf-8").splitlines()[-1] == (
        "chr1\t46\treverse-intron\tCGTA\tA\t.\t.\t."
    )


def test_sorts_expected_consequences_by_catalogue_rank(
    tmp_path: Path, reference_files: tuple[Path, Path]
) -> None:
    """Golden terms match the annotator's catalogue-rank ordering."""

    fasta, gff3 = reference_files
    cases = tmp_path / "cases.tsv"
    _write_cases(
        cases,
        [
            _row(
                "rank-order",
                "tx-forward",
                "t_start",
                "t_start",
                "NA",
                "A",
                "coding_sequence_variant,incomplete_terminal_codon_variant",
                "rank order",
            )
        ],
    )
    output_vcf = tmp_path / "corpus.vcf"
    output_json = tmp_path / "terms.json"

    write_vcf(cases, fasta, gff3, output_vcf, output_json)

    assert json.loads(output_json.read_text(encoding="utf-8"))["terms"] == {
        "tx-forward::rank-order": [
            "incomplete_terminal_codon_variant",
            "coding_sequence_variant",
        ]
    }


def test_fails_before_writing_for_an_invalid_row(tmp_path: Path, reference_files: tuple[Path, Path]) -> None:
    """Invalid coordinates are reported with their row and do not replace an existing VCF."""

    fasta, gff3 = reference_files
    cases = tmp_path / "cases.tsv"
    _write_cases(
        cases,
        [
            _row(
                "valid-row",
                "tx-forward",
                "t_start",
                "t_start",
                "NA",
                "A",
                "synonymous_variant",
                "valid row",
            ),
            _row(
                "invalid-row",
                "tx-forward",
                "exon_3_start",
                "t_end",
                "NA",
                "A",
                "invalid row",
                "invalid row",
            ),
        ],
    )
    output_vcf = tmp_path / "corpus.vcf"
    output_json = tmp_path / "terms.json"
    output_vcf.write_text("original VCF\n", encoding="utf-8")
    output_json.write_text("original terms\n", encoding="utf-8")

    with pytest.raises(GoldenVcfInputError, match=r"row 3: AlleleStart requests exon 3"):
        write_vcf(cases, fasta, gff3, output_vcf, output_json)

    assert output_vcf.read_text(encoding="utf-8") == "original VCF\n"
    assert output_json.read_text(encoding="utf-8") == "original terms\n"


@pytest.mark.parametrize(
    ("column", "coordinate"),
    [
        ("AlleleStart", "9"),
        ("AlleleEnd", "31"),
    ],
)
def test_rejects_an_integer_allele_coordinate_outside_its_transcript_span(
    tmp_path: Path,
    reference_files: tuple[Path, Path],
    column: str,
    coordinate: str,
) -> None:
    """Literal genomic coordinates must remain inside the row's transcript span."""

    fasta, gff3 = reference_files
    cases = tmp_path / "cases.tsv"
    allele_start, allele_end = (coordinate, "10") if column == "AlleleStart" else ("10", coordinate)
    _write_cases(
        cases,
        [
            _row(
                "outside-span",
                "tx-forward",
                allele_start,
                allele_end,
                "NA",
                "A",
                "synonymous_variant",
                "integer coordinate outside transcript",
            )
        ],
    )

    with pytest.raises(
        GoldenVcfInputError,
        match=rf"row 2: {column} integer coordinate {coordinate} is outside transcript span 10-30",
    ):
        write_vcf(cases, fasta, gff3, tmp_path / "corpus.vcf", tmp_path / "terms.json")


def test_rejects_a_reordered_header(tmp_path: Path, reference_files: tuple[Path, Path]) -> None:
    """The documented table order is unambiguous and enforced."""

    fasta, gff3 = reference_files
    cases = tmp_path / "cases.tsv"
    reordered = list(REQUIRED_COLUMNS)
    reordered[0], reordered[1] = reordered[1], reordered[0]
    cases.write_text("\t".join(reordered) + "\n", encoding="utf-8")

    with pytest.raises(GoldenVcfInputError, match="unexpected TSV columns"):
        write_vcf(cases, fasta, gff3, tmp_path / "corpus.vcf", tmp_path / "terms.json")


def test_sorts_input_starts_numerically_within_each_chromosome(tmp_path: Path) -> None:
    """Transcript starts sort as integers, not their lexical representation."""

    cases = tmp_path / "cases.tsv"
    _write_cases(
        cases,
        [
            _row(
                "start-100",
                "tx-100",
                "100",
                "101",
                "NA",
                "A",
                "term",
                "case",
                transcript_start=100,
                transcript_end=101,
            ),
            _row(
                "start-20",
                "tx-20",
                "20",
                "21",
                "NA",
                "A",
                "term",
                "case",
                transcript_start=20,
                transcript_end=21,
            ),
        ],
    )

    assert [row.value("VariantID") for row in _read_rows(cases)] == ["start-20", "start-100"]
